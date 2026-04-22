#!/usr/bin/env python3
"""AIMethFlow Web UI — port 9999"""
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
from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse, JSONResponse, StreamingResponse
from pydantic import BaseModel

ROOT = Path(os.environ.get("METH_WORKDIR", "/app"))
CONFIG_PATH = ROOT / "config.toml"
PYTHON = sys.executable

app = FastAPI(title="AIMethFlow")

STEPS = [
    ("preprocess",        "Step 4   预处理 fastp"),
    ("align",             "Step 5   比对 bwameth"),
    ("dedup",             "Step 6   去重 umi_tools"),
    ("target",            "Step 7   靶区过滤"),
    ("mbias",             "Step 7.5 M-bias 检查"),
    ("extract",           "Step 8   甲基化提取"),
    ("filter",            "Step 9   过滤"),
    ("diff",              "Step 10  差异分析"),
    ("visualize",         "Step 11  可视化"),
    ("feature_select",    "Step 14  特征选择"),
    ("outlier_detect",    "Step 14.5 异常检测"),
    ("classify",          "Step 15  诊断分类"),
    ("subtype_classify",  "Step 16  亚型分类"),
    ("annotate",          "Step 17  位点注释"),
    ("external_validate", "Step 18  公共数据验证"),
    ("report",            "Step 19  生成报告"),
]

# job_id -> {proc, queue, step, status}
_jobs: dict[str, dict] = {}


# ── helpers ───────────────────────────────────────────────────────────────────

def _cfg() -> dict:
    if not CONFIG_PATH.exists():
        raise HTTPException(404, f"Config not found: {CONFIG_PATH}")
    return toml.load(CONFIG_PATH)


def _latest_workdir() -> Optional[Path]:
    try:
        temp = ROOT / _cfg().get("directories", {}).get("temp_dir", "temp")
    except Exception:
        temp = ROOT / "temp"
    dirs = sorted(temp.glob("20[0-9]*"), reverse=True)
    return dirs[0] if dirs else None


def _pipeline_state() -> dict:
    wd = _latest_workdir()
    if not wd:
        return {}
    p = wd / "pipeline_state.json"
    return json.loads(p.read_text()) if p.exists() else {}


def _latest_report() -> Optional[Path]:
    try:
        rdir = ROOT / _cfg().get("directories", {}).get("result_dir", "result")
    except Exception:
        rdir = ROOT / "result"
    reports = sorted(rdir.glob("*/*_report/analysis_report.html"), reverse=True)
    return reports[0] if reports else None


# ── API ───────────────────────────────────────────────────────────────────────

@app.get("/api/config/raw")
def get_config_raw():
    if not CONFIG_PATH.exists():
        raise HTTPException(404)
    return JSONResponse({"content": CONFIG_PATH.read_text(encoding="utf-8")})


class ConfigBody(BaseModel):
    content: str


@app.post("/api/config/raw")
def save_config_raw(body: ConfigBody):
    try:
        toml.loads(body.content)
    except Exception as e:
        raise HTTPException(400, f"TOML 解析错误: {e}")
    CONFIG_PATH.write_text(body.content, encoding="utf-8")
    return {"ok": True}


@app.get("/api/steps")
def get_steps():
    state = _pipeline_state()
    wd = _latest_workdir()
    rows = []
    for key, label in STEPS:
        s = state.get(key, {})
        rows.append({
            "key": key,
            "label": label,
            "status": s.get("status", "pending"),
            "completed_at": s.get("completed_at", ""),
        })
    return {"steps": rows, "workdir": str(wd) if wd else None}


@app.post("/api/run/{step}")
def run_step(step: str, workdir: Optional[str] = None):
    if step not in {k for k, _ in STEPS}:
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
    _jobs[job_id] = {"proc": proc, "queue": log_q, "step": step, "status": "running"}
    return {"job_id": job_id}


@app.get("/api/logs/{job_id}")
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
                             headers={"Cache-Control": "no-cache",
                                      "X-Accel-Buffering": "no"})


@app.get("/api/results")
def list_results():
    try:
        rdir = ROOT / _cfg().get("directories", {}).get("result_dir", "result")
    except Exception:
        rdir = ROOT / "result"
    if not rdir.exists():
        return {"files": []}
    files = []
    for f in sorted(rdir.rglob("*"), key=lambda p: p.stat().st_mtime, reverse=True):
        if f.is_file() and not f.name.startswith("."):
            sz = f.stat().st_size
            files.append({"path": str(f.relative_to(ROOT)),
                          "size": sz,
                          "size_str": f"{sz/1024:.1f} KB" if sz < 1048576 else f"{sz/1048576:.1f} MB"})
    return {"files": files[:300]}


@app.get("/report", response_class=HTMLResponse)
def view_report():
    r = _latest_report()
    if not r:
        return HTMLResponse("<p style='padding:2rem;font-family:sans-serif'>暂无报告，请先运行 Step 19。</p>",
                            status_code=404)
    return HTMLResponse(r.read_text(encoding="utf-8"))


# ── Web UI ────────────────────────────────────────────────────────────────────

@app.get("/", response_class=HTMLResponse)
def index():
    return HTMLResponse(HTML)


HTML = """<!DOCTYPE html>
<html lang="zh">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>AIMethFlow</title>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css">
<style>
body{background:#f8f9fa}
.navbar-brand{font-weight:700;letter-spacing:.05em}
#log-box{background:#1e1e1e;color:#d4d4d4;font-family:monospace;font-size:.8rem;
         height:320px;overflow-y:auto;padding:.75rem;border-radius:4px;white-space:pre-wrap}
.step-row td{vertical-align:middle}
.badge-pending{background:#6c757d}
.badge-done{background:#198754}
.badge-running{background:#0d6efd}
.badge-error{background:#dc3545}
.badge-skipped{background:#adb5bd}
textarea.toml{font-family:monospace;font-size:.82rem;height:70vh}
#results-list{font-size:.85rem}
.iframe-wrap{height:calc(100vh - 120px);border:none}
.tab-content{padding-top:1rem}
</style>
</head>
<body>
<nav class="navbar navbar-dark bg-dark px-3">
  <span class="navbar-brand">⚗ AIMethFlow</span>
  <span class="text-secondary small" id="nav-workdir"></span>
</nav>

<div class="container-fluid mt-2">
<ul class="nav nav-tabs" id="tabs">
  <li class="nav-item"><a class="nav-link active" data-tab="pipeline" href="#">流程</a></li>
  <li class="nav-item"><a class="nav-link" data-tab="config" href="#">配置</a></li>
  <li class="nav-item"><a class="nav-link" data-tab="results" href="#">结果</a></li>
  <li class="nav-item"><a class="nav-link" data-tab="report" href="#">报告</a></li>
</ul>

<!-- 流程 -->
<div id="tab-pipeline" class="tab-content">
  <div class="row">
    <div class="col-lg-5">
      <table class="table table-sm table-hover" id="step-table">
        <thead><tr><th>步骤</th><th>状态</th><th>完成时间</th><th></th></tr></thead>
        <tbody id="step-tbody"></tbody>
      </table>
    </div>
    <div class="col-lg-7">
      <div class="d-flex justify-content-between align-items-center mb-1">
        <small class="text-muted" id="log-title">— 日志 —</small>
        <button class="btn btn-sm btn-outline-secondary" onclick="clearLog()">清空</button>
      </div>
      <div id="log-box"></div>
    </div>
  </div>
</div>

<!-- 配置 -->
<div id="tab-config" class="tab-content d-none">
  <div class="d-flex gap-2 mb-2">
    <button class="btn btn-primary btn-sm" onclick="saveConfig()">保存</button>
    <button class="btn btn-outline-secondary btn-sm" onclick="loadConfig()">重新加载</button>
    <span class="text-muted small align-self-center" id="cfg-status"></span>
  </div>
  <textarea class="form-control toml" id="cfg-editor" spellcheck="false"></textarea>
</div>

<!-- 结果 -->
<div id="tab-results" class="tab-content d-none">
  <button class="btn btn-sm btn-outline-secondary mb-2" onclick="loadResults()">刷新</button>
  <div id="results-list"></div>
</div>

<!-- 报告 -->
<div id="tab-report" class="tab-content d-none">
  <div class="d-flex gap-2 mb-2">
    <button class="btn btn-sm btn-outline-secondary" onclick="reloadReport()">刷新报告</button>
  </div>
  <iframe id="report-frame" class="iframe-wrap w-100" src="/report"></iframe>
</div>
</div>

<script>
const STATUS_CLASS = {pending:'badge-pending',done:'badge-done',running:'badge-running',
                      error:'badge-error',skipped:'badge-skipped'};
const STATUS_LABEL = {pending:'待执行',done:'完成',running:'执行中',error:'失败',skipped:'跳过'};
let currentJob = null;

// ── tabs ─────────────────────────────────────────────────────────────────────
document.querySelectorAll('[data-tab]').forEach(a => {
  a.addEventListener('click', e => {
    e.preventDefault();
    document.querySelectorAll('[data-tab]').forEach(x => x.classList.remove('active'));
    a.classList.add('active');
    document.querySelectorAll('.tab-content').forEach(x => x.classList.add('d-none'));
    const t = document.getElementById('tab-' + a.dataset.tab);
    t.classList.remove('d-none');
    if (a.dataset.tab === 'config') loadConfig();
    if (a.dataset.tab === 'results') loadResults();
  });
});

// ── pipeline ──────────────────────────────────────────────────────────────────
async function loadSteps() {
  const r = await fetch('/api/steps').then(r => r.json());
  if (r.workdir) document.getElementById('nav-workdir').textContent = r.workdir;
  const tbody = document.getElementById('step-tbody');
  tbody.innerHTML = '';
  r.steps.forEach(s => {
    const cls = STATUS_CLASS[s.status] || 'badge-pending';
    const lbl = STATUS_LABEL[s.status] || s.status;
    const ts = s.completed_at ? `<small class="text-muted">${s.completed_at.slice(0,19)}</small>` : '';
    tbody.innerHTML += `<tr class="step-row">
      <td>${s.label}</td>
      <td><span class="badge ${cls}">${lbl}</span></td>
      <td>${ts}</td>
      <td><button class="btn btn-sm btn-outline-primary py-0" onclick="runStep('${s.key}','${s.label}')">▶</button></td>
    </tr>`;
  });
}

async function runStep(key, label) {
  if (currentJob) {
    if (!confirm('当前有任务正在执行，确定要新建任务吗？')) return;
  }
  document.getElementById('log-title').textContent = `▶ ${label}`;
  const logBox = document.getElementById('log-box');
  logBox.innerHTML = '';
  const r = await fetch(`/api/run/${key}`, {method:'POST'}).then(r => r.json());
  if (r.job_id) streamLogs(r.job_id);
}

function streamLogs(jobId) {
  currentJob = jobId;
  const logBox = document.getElementById('log-box');
  const es = new EventSource(`/api/logs/${jobId}`);
  es.onmessage = e => {
    if (e.data.startsWith('__DONE_')) {
      const rc = parseInt(e.data.replace('__DONE_',''));
      logBox.innerHTML += `\\n[退出码: ${rc}]\\n`;
      es.close();
      currentJob = null;
      loadSteps();
      return;
    }
    const line = JSON.parse(e.data);
    logBox.innerHTML += escHtml(line) + '\\n';
    logBox.scrollTop = logBox.scrollHeight;
  };
  es.onerror = () => { es.close(); currentJob = null; };
}

function clearLog() { document.getElementById('log-box').innerHTML = ''; }
function escHtml(s) { return s.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;'); }

// ── config ────────────────────────────────────────────────────────────────────
async function loadConfig() {
  const r = await fetch('/api/config/raw').then(r => r.json());
  document.getElementById('cfg-editor').value = r.content;
  document.getElementById('cfg-status').textContent = '';
}

async function saveConfig() {
  const content = document.getElementById('cfg-editor').value;
  const r = await fetch('/api/config/raw', {
    method:'POST', headers:{'Content-Type':'application/json'},
    body: JSON.stringify({content})
  });
  const j = await r.json();
  const el = document.getElementById('cfg-status');
  if (r.ok) { el.style.color='green'; el.textContent = '✓ 已保存'; }
  else       { el.style.color='red';   el.textContent = '✗ ' + (j.detail || '保存失败'); }
}

// ── results ───────────────────────────────────────────────────────────────────
async function loadResults() {
  const r = await fetch('/api/results').then(r => r.json());
  const el = document.getElementById('results-list');
  if (!r.files.length) { el.innerHTML = '<p class="text-muted">暂无结果文件</p>'; return; }
  el.innerHTML = '<table class="table table-sm"><thead><tr><th>文件</th><th>大小</th></tr></thead><tbody>'
    + r.files.map(f => `<tr><td><code>${f.path}</code></td><td class="text-muted">${f.size_str}</td></tr>`).join('')
    + '</tbody></table>';
}

// ── report ────────────────────────────────────────────────────────────────────
function reloadReport() {
  const f = document.getElementById('report-frame');
  f.src = '/report?' + Date.now();
}

// ── init ──────────────────────────────────────────────────────────────────────
loadSteps();
setInterval(loadSteps, 8000);
</script>
</body>
</html>"""


if __name__ == "__main__":
    uvicorn.run("webui:app", host="0.0.0.0", port=9999, reload=False,
                app_dir=str(Path(__file__).parent))
