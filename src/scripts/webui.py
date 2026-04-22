#!/usr/bin/env python3
"""AIMethFlow Web UI — port 9999, pages: /meth/  API: /methapi/"""
import json
import os
import queue
import subprocess
import sys
import threading
import uuid
from pathlib import Path
from typing import Optional

import toml
import uvicorn
from fastapi import APIRouter, FastAPI, HTTPException
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, RedirectResponse, StreamingResponse
from pydantic import BaseModel

ROOT = Path(os.environ.get("METH_WORKDIR", "/app"))
CONFIG_PATH = ROOT / "config.toml"
PYTHON = sys.executable

app = FastAPI(title="AIMethFlow")
page_router = APIRouter(prefix="/meth")   # HTML pages
api_router  = APIRouter(prefix="/methapi") # REST / SSE

STEPS = [
    ("preprocess",        "Step 4    预处理 fastp"),
    ("align",             "Step 5    比对 bwameth"),
    ("dedup",             "Step 6    去重 umi_tools"),
    ("target",            "Step 7    靶区过滤"),
    ("mbias",             "Step 7.5  M-bias 检查"),
    ("extract",           "Step 8    甲基化提取"),
    ("filter",            "Step 9    过滤"),
    ("diff",              "Step 10   差异分析"),
    ("visualize",         "Step 11   可视化"),
    ("feature_select",    "Step 14   特征选择"),
    ("outlier_detect",    "Step 14.5 异常检测"),
    ("classify",          "Step 15   诊断分类"),
    ("subtype_classify",  "Step 16   亚型分类"),
    ("annotate",          "Step 17   位点注释"),
    ("external_validate", "Step 18   公共数据验证"),
    ("report",            "Step 19   生成报告"),
]
_STEP_KEYS = {k for k, _ in STEPS}

# job_id -> {proc, queue, step, status, workdir}
_jobs: dict[str, dict] = {}

_HTML = (Path(__file__).parent / "webui.html").read_text(encoding="utf-8")


# ── helpers ───────────────────────────────────────────────────────────────────

def _cfg() -> dict:
    if not CONFIG_PATH.exists():
        raise HTTPException(404, f"Config not found: {CONFIG_PATH}")
    return toml.load(CONFIG_PATH)


def _temp_dir() -> Path:
    try:
        return ROOT / _cfg().get("directories", {}).get("temp_dir", "temp")
    except Exception:
        return ROOT / "temp"


def _result_dir() -> Path:
    try:
        return ROOT / _cfg().get("directories", {}).get("result_dir", "result")
    except Exception:
        return ROOT / "result"


def _latest_workdir() -> Optional[Path]:
    dirs = sorted(_temp_dir().glob("20[0-9]*"), reverse=True)
    return dirs[0] if dirs else None


def _latest_report() -> Optional[Path]:
    reports = sorted(_result_dir().glob("*/*_report/analysis_report.html"), reverse=True)
    return reports[0] if reports else None


# ── config API ────────────────────────────────────────────────────────────────

@api_router.get("/api/config/raw")
def get_config_raw():
    if not CONFIG_PATH.exists():
        raise HTTPException(404)
    return JSONResponse({"content": CONFIG_PATH.read_text(encoding="utf-8")})


class ConfigBody(BaseModel):
    content: str


@api_router.post("/api/config/raw")
def save_config_raw(body: ConfigBody):
    try:
        toml.loads(body.content)
    except Exception as e:
        raise HTTPException(400, f"TOML 解析错误: {e}")
    CONFIG_PATH.write_text(body.content, encoding="utf-8")
    return {"ok": True}


# ── pipeline API ──────────────────────────────────────────────────────────────

@api_router.get("/api/workdirs")
def list_workdirs():
    dirs = sorted(_temp_dir().glob("20[0-9]*"), reverse=True)
    return {"workdirs": [str(d.relative_to(ROOT)).replace("\\", "/") for d in dirs]}


@api_router.get("/api/steps")
def get_steps(workdir: Optional[str] = None):
    wd = (ROOT / workdir) if workdir else _latest_workdir()
    completed: set = set()
    if wd and wd.exists():
        p = wd / "pipeline_state.json"
        if p.exists():
            completed = set(json.loads(p.read_text()).get("completed_steps", []))

    running_step = next(
        (j["step"] for j in _jobs.values()
         if j["status"] == "running" and j.get("workdir") == (workdir or "")),
        None,
    )

    rows = []
    for key, label in STEPS:
        if key == running_step:
            status = "running"
        elif key in completed:
            status = "done"
        else:
            status = "pending"
        rows.append({"key": key, "label": label, "status": status})
    return {"steps": rows, "workdir": str(wd) if wd else None}


@api_router.post("/api/run/{step}")
def run_step(step: str, workdir: Optional[str] = None):
    if step not in _STEP_KEYS:
        raise HTTPException(400, f"未知步骤: {step}")

    wd = workdir or str(_latest_workdir() or "")
    cmd = [PYTHON, str(ROOT / "src/scripts/run_analysis.py"),
           "--plasma", "--step", step, "--config", str(CONFIG_PATH)]
    if wd:
        cmd += ["--workdir", wd]

    log_q: queue.Queue = queue.Queue()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, bufsize=1, cwd=str(ROOT))

    def _read():
        for line in proc.stdout:
            log_q.put(line)
        proc.wait()
        log_q.put(None)

    threading.Thread(target=_read, daemon=True).start()
    job_id = uuid.uuid4().hex[:8]
    _jobs[job_id] = {"proc": proc, "queue": log_q, "step": step,
                     "status": "running", "workdir": workdir or ""}
    return {"job_id": job_id}


@api_router.get("/api/logs/{job_id}")
def stream_logs(job_id: str):
    if job_id not in _jobs:
        raise HTTPException(404)
    job = _jobs[job_id]

    def _gen():
        q = job["queue"]
        while True:
            line = q.get()
            if line is None:
                rc = job["proc"].returncode
                job["status"] = "done" if rc == 0 else "error"
                yield f"data: __DONE_{rc}__\n\n"
                break
            yield f"data: {json.dumps(line.rstrip())}\n\n"

    return StreamingResponse(_gen(), media_type="text/event-stream",
                             headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"})


# ── results / reports API ─────────────────────────────────────────────────────

@api_router.get("/api/results")
def list_results():
    rdir = _result_dir()
    if not rdir.exists():
        return {"files": []}
    files = []
    for f in sorted(rdir.rglob("*"), key=lambda p: p.stat().st_mtime, reverse=True):
        if f.is_file() and not f.name.startswith("."):
            sz = f.stat().st_size
            files.append({
                "path": str(f.relative_to(ROOT)),
                "size": sz,
                "size_str": f"{sz/1024:.1f} KB" if sz < 1_048_576 else f"{sz/1_048_576:.1f} MB",
            })
    return {"files": files[:300]}


@api_router.get("/api/reports")
def list_reports():
    reports = []
    for p in sorted(_result_dir().glob("*/*_report/analysis_report.html"), reverse=True):
        reports.append({
            "name": p.parent.name,
            "report_path": str(p.relative_to(ROOT)).replace("\\", "/"),
        })
    return {"reports": reports}


# ── page routes ───────────────────────────────────────────────────────────────

@page_router.get("/report-files/{file_path:path}")
def serve_report_file(file_path: str):
    target = (ROOT / file_path).resolve()
    if not str(target).startswith(str(ROOT.resolve())):
        raise HTTPException(403)
    if not target.exists() or not target.is_file():
        raise HTTPException(404)
    return FileResponse(target)


@page_router.get("/report")
def view_report():
    r = _latest_report()
    if not r:
        return HTMLResponse("<p style='padding:2rem;font-family:sans-serif'>暂无报告，请先运行 Step 19。</p>",
                            status_code=404)
    return RedirectResponse(f"/meth/report-files/{str(r.relative_to(ROOT)).replace(chr(92), '/')}")


@page_router.get("/", response_class=HTMLResponse)
def index():
    return HTMLResponse(_HTML)


app.include_router(page_router)
app.include_router(api_router)


if __name__ == "__main__":
    uvicorn.run("webui:app", host="0.0.0.0", port=9999, reload=False,
                app_dir=str(Path(__file__).parent))
