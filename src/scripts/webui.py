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
from fastapi.responses import (
    FileResponse,
    HTMLResponse,
    JSONResponse,
    RedirectResponse,
    StreamingResponse,
)
from pydantic import BaseModel

def _load_dotenv(path: Path) -> dict:
    """Parse a .env file and return key-value pairs (no shell expansion)."""
    result = {}
    if not path.exists():
        return result
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, _, v = line.partition("=")
        result[k.strip()] = v.strip()
    return result


_DOTENV = _load_dotenv(Path(__file__).parent.parent.parent / ".env")

ROOT = Path(_DOTENV.get("METH_WORKDIR") or os.environ.get("METH_WORKDIR") or "/app").resolve()
CONFIG_PATH = ROOT / "config.toml"
PYTHON = sys.executable
print(f"Using ROOT={ROOT}, CONFIG_PATH={CONFIG_PATH}, Python={PYTHON}")

app = FastAPI(title="AIMethFlow")
page_router = APIRouter(prefix="/meth")  # HTML pages
api_router = APIRouter(prefix="/methapi")  # REST / SSE

STEPS = [
    ("preprocess",        "Step 4    Preprocessing (fastp)"),
    ("align",             "Step 5    Alignment (bwameth)"),
    ("dedup",             "Step 6    Deduplication (umi_tools)"),
    ("target",            "Step 7    Target Region Filter"),
    ("mbias",             "Step 8    M-bias Check"),
    ("extract",           "Step 9    Methylation Extraction"),
    ("filter",            "Step 10   Filtering"),
    ("diff",              "Step 11   Differential Analysis"),
    ("visualize",         "Step 12   Visualization"),
    ("feature_select",    "Step 13   Feature Selection"),
    ("outlier_detect",    "Step 14   Outlier Detection"),
    ("classify",          "Step 15   Diagnostic Classification"),
    ("subtype_classify",  "Step 16   Subtype Classification"),
    ("temporal_validate", "Step 17   Temporal Validation (plasma)"),
    ("annotate",          "Step 18   Site Annotation"),
    ("external_validate", "Step 19   External Validation"),
    ("report",            "Step 20   Report Generation"),
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
        raise HTTPException(400, f"TOML parse error: {e}")
    CONFIG_PATH.write_text(body.content, encoding="utf-8")
    return {"ok": True}


# ── pipeline API ──────────────────────────────────────────────────────────────


@api_router.get("/api/workdirs")
def list_workdirs():
    dirs = sorted(_temp_dir().glob("20[0-9]*"), reverse=True)
    return {"workdirs": [str(d) for d in dirs]}


@api_router.get("/api/steps")
def get_steps(workdir: Optional[str] = None):
    wd = Path(workdir) if workdir else _latest_workdir()
    completed: set = set()
    if wd and wd.exists():
        p = wd / "pipeline_state.json"
        if p.exists():
            completed = set(json.loads(p.read_text()).get("completed_steps", []))

    running_step = next(
        (
            j["step"]
            for j in _jobs.values()
            if j["status"] == "running" and j.get("workdir") == (workdir or "")
        ),
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
def run_step(step: str, workdir: Optional[str] = None, mode: str = "plasma"):
    if step not in _STEP_KEYS:
        raise HTTPException(400, f"Unknown step: {step}")
    if mode not in ("plasma", "tissue"):
        raise HTTPException(400, f"Invalid mode: {mode!r} — must be 'plasma' or 'tissue'")

    wd = workdir or str(_latest_workdir() or "")
    cmd = [
        PYTHON,
        str(ROOT / "src/scripts/run_analysis.py"),
        "--mode", mode,
        "--step",
        step,
        "--config",
        str(CONFIG_PATH),
    ]
    if wd:
        cmd += ["--workdir", wd]

    log_q: queue.Queue = queue.Queue()
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, cwd=str(ROOT)
    )

    def _read():
        for line in proc.stdout:
            log_q.put(line)
        proc.wait()
        log_q.put(None)

    threading.Thread(target=_read, daemon=True).start()
    job_id = uuid.uuid4().hex[:8]
    _jobs[job_id] = {
        "proc": proc,
        "queue": log_q,
        "step": step,
        "status": "running",
        "workdir": workdir or "",
    }
    return {"job_id": job_id}


@api_router.post("/api/jobs/{job_id}/kill")
def kill_job(job_id: str):
    if job_id not in _jobs:
        raise HTTPException(404)
    job = _jobs[job_id]
    proc = job["proc"]
    if proc.poll() is None:
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except Exception:
            proc.kill()
    job["status"] = "killed"
    return {"ok": True}


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

    return StreamingResponse(
        _gen(),
        media_type="text/event-stream",
        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
    )


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
            files.append(
                {
                    "path": str(f.relative_to(rdir)),
                    "size": sz,
                    "size_str": f"{sz/1024:.1f} KB" if sz < 1_048_576 else f"{sz/1_048_576:.1f} MB",
                }
            )
    return {"files": files[:300]}


@api_router.get("/api/reports")
def list_reports():
    reports = []
    for p in sorted(_result_dir().glob("*/*_report/analysis_report.html"), reverse=True):
        # Strip leading "/" so the path can be embedded in a URL segment
        reports.append(
            {
                "name": p.parent.name,
                "report_path": str(p).lstrip("/"),
            }
        )
    return {"reports": reports}


# ── page routes ───────────────────────────────────────────────────────────────


@page_router.get("/report-files/{file_path:path}")
def serve_report_file(file_path: str):
    # Restore the leading "/" that was stripped for URL embedding
    target = Path("/" + file_path).resolve()
    safe = {str(_result_dir().resolve()), str(_temp_dir().resolve()), str(ROOT.resolve())}
    if not any(str(target).startswith(r) for r in safe):
        raise HTTPException(403)
    if not target.exists() or not target.is_file():
        raise HTTPException(404)
    return FileResponse(target)


@page_router.get("/report")
def view_report():
    r = _latest_report()
    if not r:
        return HTMLResponse(
            "<p style='padding:2rem;font-family:sans-serif'>No report found. Run Step 19 first.</p>",
            status_code=404,
        )
    return RedirectResponse(f"/meth/report-files/{str(r).lstrip('/')}")


@page_router.get("/", response_class=HTMLResponse)
def index():
    return HTMLResponse(_HTML)


app.include_router(page_router)
app.include_router(api_router)


if __name__ == "__main__":
    uvicorn.run(
        "webui:app", host="0.0.0.0", port=9999, reload=True, app_dir=str(Path(__file__).parent)
    )
