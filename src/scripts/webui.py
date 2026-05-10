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
from fastapi import APIRouter, FastAPI, File, Form, HTTPException, UploadFile
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

# Pipeline step list — must stay in sync with PIPELINE_STEPS in run_analysis.py.
# (We don't import run_analysis here because it has top-level side effects:
#  load_config() reads config.toml and setup_logger() creates a log file.)
# Each entry is (key, "Step <num>  <description>", "<module>").
STEPS = [
    ("data_check",        "Step 1   Validate config",                 "A · Setup"),
    ("ref_index",         "Step 2   bwameth hg38 index (one-time)",   "A · Setup"),
    ("panel_derive",      "Step 3   Derive panel BED",                "A · Setup"),
    ("preprocess",        "Step 4   Preprocessing (fastp + UMI)",     "B · Preprocessing"),
    ("align",             "Step 5   Alignment (bwameth)",             "B · Preprocessing"),
    ("dedup",             "Step 6   Deduplication (umi_tools)",       "B · Preprocessing"),
    ("target",            "Step 7   On-target filter",                "B · Preprocessing"),
    ("mbias",             "Step 8   M-bias check",                    "B · Preprocessing"),
    ("extract",           "Step 9   Methylation extraction",          "C · Methylation"),
    ("filter",            "Step 10  Depth + repeat filter",           "C · Methylation"),
    ("diff",              "Step 11  Differential methylation",        "C · Methylation"),
    ("visualize",         "Step 12  Density / scatter / heatmap",     "C · Methylation"),
    ("outlier_detect",    "Step 13  Sample outlier detection",        "D · QC + Features"),
    ("feature_select",    "Step 14  Feature selection (LASSO + RF)",  "D · QC + Features"),
    ("fs_benchmark",      "Step 14b Feature selection ablation",      "D · QC + Features"),
    ("classify",          "Step 15  Diagnostic classification",       "E · Classification"),
    ("final_model",       "Step 15c Train production models (full cohort)", "E · Classification"),
    ("subtype_classify",  "Step 16  Subtype classification",          "E · Classification"),
    ("temporal_validate", "Step 17  Temporal split validation",       "E · Classification"),
    ("annotate",          "Step 18  Site annotation",                 "E · Classification"),
    ("external_validate", "Step 19  External validation",             "F · Output"),
    ("report",            "Step 20  HTML + JSON report",              "F · Output"),
]
_STEP_KEYS = {s[0] for s in STEPS}

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
    for key, label, module in STEPS:
        if key == running_step:
            status = "running"
        elif key in completed:
            status = "done"
        else:
            status = "pending"
        rows.append({"key": key, "label": label, "module": module, "status": status})
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


# ── Aided Diagnosis API ───────────────────────────────────────────────────────
#
# Endpoints support the "Aided Diagnosis" tab: pick a pre-trained workdir
# (result/<wd>/model/), upload a bedGraph (or fq.gz, future), enter optional
# CA-125, and get per-model probabilities + per-CpG contribution against the
# training cohort medians.

def _result_root() -> Path:
    """Return absolute path to result/, regardless of cwd."""
    return ROOT / "result"


def _list_model_workdirs():
    """Walk result/*/model/ and return entries with manifest summary."""
    out = []
    root = _result_root()
    if not root.exists():
        return out
    for wd_dir in sorted(root.iterdir(), reverse=True):
        if not wd_dir.is_dir():
            continue
        model_dir = wd_dir / "model"
        if not model_dir.is_dir():
            continue
        manifest_path = model_dir / "manifest.json"
        if not manifest_path.exists():
            continue
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except Exception:
            continue
        out.append({
            "workdir":     wd_dir.name,
            "n_models":    len(manifest.get("trained_models", [])),
            "n_features":  manifest.get("panel", {}).get("n_features", 0),
            "ca125":       bool(manifest.get("ca125_available", False)),
            "n_train":     manifest.get("training_cohort", {}).get("n_total", 0),
        })
    return out


@api_router.get("/api/diagnosis/model_workdirs")
def get_diagnosis_model_workdirs():
    """List workdirs that have a trained model/ directory (output of Step 15c)."""
    return {"model_workdirs": _list_model_workdirs()}


def _parse_bedgraph_to_feature_vector(file_bytes: bytes, feature_names, filename: str):
    """
    Parse a bedGraph file (optionally gzipped) and extract β-values for the
    panel CpG sites in `feature_names` (formatted as 'chrom:pos').

    Returns (feature_dict, n_resolved) where feature_dict[name] = β (0-1) or
    np.nan when the site isn't present in the bedGraph.
    """
    import gzip, io
    import numpy as np

    if filename.endswith(".gz"):
        try:
            text = gzip.decompress(file_bytes).decode("utf-8", errors="replace")
        except Exception as exc:
            raise HTTPException(400, f"Failed to gunzip bedGraph: {exc}")
    else:
        text = file_bytes.decode("utf-8", errors="replace")

    # Build a {(chrom, pos): beta} map; bedGraph schema: chrom start end meth_pct ...
    # MethylDackel writes meth_pct as 0-100 (β × 100). Convert to 0-1 here.
    site_map = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "track", "chrom\t")):
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        try:
            chrom = parts[0]
            pos   = int(parts[1])
            val   = float(parts[3])
        except (ValueError, TypeError):
            continue
        # Heuristic: values >1 → 0-100 percent → divide by 100; values ≤1 → already 0-1
        if val > 1.0:
            val = val / 100.0
        site_map[(chrom, pos)] = val

    feat_dict = {}
    n_resolved = 0
    for name in feature_names:
        if ":" not in name:
            feat_dict[name] = float("nan")
            continue
        chrom, _, pos_str = name.rpartition(":")
        try:
            pos = int(pos_str)
        except ValueError:
            feat_dict[name] = float("nan")
            continue
        v = site_map.get((chrom, pos), float("nan"))
        feat_dict[name] = v
        if v == v:  # not NaN
            n_resolved += 1
    return feat_dict, n_resolved


@api_router.post("/api/diagnosis/predict")
async def predict_diagnosis(
    workdir: str = Form(...),
    file: UploadFile = File(...),
    ca125: Optional[str] = Form(None),
):
    """
    Run a diagnosis using pre-trained models in result/<workdir>/model/.
    Currently supports bedGraph input. fq.gz support is stubbed (returns 400
    pending Step 9 background-run integration).
    """
    import joblib
    import numpy as np

    # 1. Locate model directory
    model_dir = _result_root() / workdir / "model"
    if not model_dir.is_dir():
        raise HTTPException(404, f"Model directory not found: result/{workdir}/model")
    manifest_path = model_dir / "manifest.json"
    pre_path      = model_dir / "preprocessing.joblib"
    if not (manifest_path.exists() and pre_path.exists()):
        raise HTTPException(404, "manifest.json or preprocessing.joblib missing in model dir")

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    pre      = joblib.load(pre_path)
    feature_names = list(pre["feature_names"])
    col_medians   = np.asarray(pre["col_medians"], dtype=float)

    # 2. Validate / parse uploaded file
    fname = (file.filename or "").lower()
    file_bytes = await file.read()
    if fname.endswith((".bedgraph", ".bedgraph.gz", ".bg", ".bg.gz")) or \
       (fname.endswith(".gz") and "bedgraph" in fname):
        feat_dict, n_resolved = _parse_bedgraph_to_feature_vector(file_bytes, feature_names, fname)
    elif fname.endswith((".fq", ".fq.gz", ".fastq", ".fastq.gz")):
        raise HTTPException(
            400,
            "FASTQ input is not yet wired into the diagnosis path. Workaround: "
            "run --step preprocess/align/dedup/target/extract on the FASTQ to "
            "produce a *.bedGraph, then upload that file."
        )
    else:
        raise HTTPException(400, f"Unsupported file type: {file.filename!r}. "
                                  "Expected .bedGraph or .bedGraph.gz.")

    # 3. Build feature vector + per-feature contribution table
    x_raw  = np.array([feat_dict[name] for name in feature_names], dtype=float)
    x_imp  = np.where(np.isnan(x_raw), col_medians, x_raw)
    x_sc   = pre["scaler"].transform(x_imp.reshape(1, -1))

    features_payload = []
    for i, name in enumerate(feature_names):
        beta_raw   = x_raw[i]
        train_med  = float(col_medians[i])
        is_missing = bool(np.isnan(beta_raw))
        beta_used  = train_med if is_missing else float(beta_raw)
        features_payload.append({
            "name":            name,
            "beta":            None if is_missing else float(beta_raw),
            "training_median": train_med,
            "delta":           beta_used - train_med,   # positive → above cohort median
            "missing":         is_missing,
        })

    # 4. Optional CA-125 augmented vector for fusion models
    ca_val = None
    x_aug_sc = None
    if ca125 not in (None, ""):
        try:
            ca_val = float(ca125)
        except ValueError:
            raise HTTPException(400, f"Invalid ca125 value: {ca125!r}")
    if pre.get("ca125_scaler") is not None:
        ca_for_scale = ca_val if ca_val is not None else pre.get("ca125_median", 0.0)
        ca_sc = pre["ca125_scaler"].transform(np.array([[ca_for_scale]], dtype=float))
        x_aug_sc = np.hstack([x_sc, ca_sc])

    # 5. Run each available model
    import warnings as _warnings
    models_payload = []
    for name in manifest.get("trained_models", []):
        clf_path = model_dir / f"{name}.joblib"
        if not clf_path.exists():
            continue
        try:
            clf = joblib.load(clf_path)
        except Exception:
            continue
        is_fusion = name.endswith("+CA125")
        if is_fusion:
            if x_aug_sc is None or ca_val is None:
                # User didn't provide CA-125 → skip fusion models silently
                continue
            x_in = x_aug_sc
        else:
            x_in = x_sc
        try:
            with _warnings.catch_warnings():
                _warnings.filterwarnings(
                    "ignore", message="X does not have valid feature names",
                    category=UserWarning,
                )
                p = float(clf.predict_proba(x_in)[0, 1])
        except Exception as exc:
            print(f"[diagnosis] {name} predict failed: {exc}")
            continue
        models_payload.append({
            "name":            name,
            "uses_ca125":      is_fusion,
            "prob_malignant":  p,
            "predicted_label": "malignant" if p >= 0.5 else "benign",
        })

    # 6. Consensus (majority vote across all available models)
    n_voters = len(models_payload)
    mal_votes = sum(1 for m in models_payload if m["predicted_label"] == "malignant")
    mean_prob = (sum(m["prob_malignant"] for m in models_payload) / n_voters) if n_voters else 0.0
    if n_voters == 0:
        consensus = {"label": "unknown", "n_voters": 0, "malignant_votes": 0, "mean_prob": 0.0}
    elif mal_votes > n_voters / 2:
        consensus = {"label": "malignant", "n_voters": n_voters, "malignant_votes": mal_votes, "mean_prob": mean_prob}
    elif mal_votes < n_voters / 2:
        consensus = {"label": "benign",    "n_voters": n_voters, "malignant_votes": mal_votes, "mean_prob": mean_prob}
    else:
        consensus = {"label": "borderline", "n_voters": n_voters, "malignant_votes": mal_votes, "mean_prob": mean_prob}

    return {
        "workdir":            workdir,
        "n_panel_features":   len(feature_names),
        "n_features_resolved": int(sum(1 for f in features_payload if not f["missing"])),
        "ca125_supplied":     ca_val,
        "models":             models_payload,
        "consensus":          consensus,
        "features":           features_payload,
    }


app.include_router(page_router)
app.include_router(api_router)


if __name__ == "__main__":
    uvicorn.run(
        "webui:app", host="0.0.0.0", port=9999, reload=True, app_dir=str(Path(__file__).parent)
    )
