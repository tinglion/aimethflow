#!/usr/bin/env python3
"""report.helpers — shared low-level utilities for the report package.

Holds parsers (fastp/dedup/bedGraph/diff), numeric formatters, image
copy/embed helpers, the CSV reader, and the JS table-pagination widget.

All public symbols (no leading underscore) and the underscore-prefixed
helpers are re-exported by `scripts.modules.report.__init__` so existing
call sites that do ``from scripts.modules.report_utils import X`` keep
working through the shim.
"""

import os
import re
import json
import glob
import shutil
import logging
from datetime import datetime

try:
    import pandas as pd
    _HAS_PANDAS = True
except ImportError:
    _HAS_PANDAS = False

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────
#  Data-reading helpers
# ─────────────────────────────────────────────

def _parse_fastp_json(path):
    try:
        with open(path, encoding="utf-8") as f:
            d = json.load(f)
        sm = d.get("summary", {})
        bf = sm.get("before_filtering", {})
        af = sm.get("after_filtering", {})
        dup_rate = d.get("duplication", {}).get("rate", None)
        return {
            "total_reads_before": bf.get("total_reads", 0),
            "total_bases_before": bf.get("total_bases", 0),
            "q30_before":         bf.get("q30_rate", 0),
            "total_reads_after":  af.get("total_reads", 0),
            "total_bases_after":  af.get("total_bases", 0),
            "q30_after":          af.get("q30_rate", 0),
            "passed_filter_reads": d.get("filtering_result", {}).get("passed_filter_reads", 0),
            "low_quality_reads":   d.get("filtering_result", {}).get("low_quality_reads", 0),
            "too_short_reads":     d.get("filtering_result", {}).get("too_short_reads", 0),
            "dup_rate": dup_rate,
        }
    except Exception as e:
        logger.warning(f"Cannot parse fastp JSON {path}: {e}")
        return {}


def _parse_dedup_stats(path):
    result = {}
    try:
        with open(path, encoding="utf-8") as f:
            content = f.read()
        m = re.search(r"Input Reads:\s*([\d,]+)", content)
        if m:
            result["input_reads"] = int(m.group(1).replace(",", ""))
        m = re.search(r"Read pairs:\s*([\d,]+)", content)
        if m:
            result["input_pairs"] = int(m.group(1).replace(",", ""))
        m = re.search(r"Number of reads out:\s*([\d,]+)", content)
        if m:
            result["output_reads"] = int(m.group(1).replace(",", ""))
        if "input_pairs" in result and "output_reads" in result and result["input_pairs"] > 0:
            result["dup_rate"] = 1.0 - result["output_reads"] / result["input_pairs"]
    except Exception as e:
        logger.warning(f"Cannot parse dedup stats {path}: {e}")
    return result


def _count_bedgraph_sites(path):
    try:
        count = 0
        meth_sum = 0.0
        with open(path, encoding="utf-8") as f:
            for line in f:
                if line.startswith("#") or line.startswith("track"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    count += 1
                    try:
                        meth_sum += float(parts[3])
                    except ValueError:
                        pass
        avg_meth = meth_sum / count if count > 0 else 0
        return count, avg_meth
    except Exception:
        return 0, 0.0


def _parse_diff_results(path):
    rows = []
    try:
        with open(path, encoding="utf-8") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 5:
                    rows.append({
                        "chrom": parts[0],
                        "pos":   int(parts[1]),
                        "meth1": float(parts[2]),
                        "meth2": float(parts[3]),
                        "diff":  float(parts[4]),
                    })
    except Exception as e:
        logger.warning(f"Cannot parse diff results {path}: {e}")
    return rows


def _fmt_num(n):
    if isinstance(n, int):
        return f"{n:,}"
    return str(n)

def _fmt_pct(r, digits=1):
    if r is None:
        return "N/A"
    return f"{r * 100:.{digits}f}%"

def _fmt_mb(b):
    return f"{b / 1024**2:.0f} MB" if b else "—"

def _fmt_gb(b):
    return f"{b / 1024**3:.2f} GB" if b else "—"

def _safe_pct(num, denom, digits=1):
    """Return percentage string; returns '—' when denominator is 0."""
    if not denom:
        return "—"
    return f"{num / denom * 100:.{digits}f}%"


def _copy_img(src_path, imgs_dir):
    """Copy *src_path* to *imgs_dir* and return a relative URL ``imgs/<fname>``.
    Returns None if the source file does not exist."""
    if not src_path or not os.path.exists(src_path):
        return None
    os.makedirs(imgs_dir, exist_ok=True)
    fname = os.path.basename(src_path)
    shutil.copy2(src_path, os.path.join(imgs_dir, fname))
    return f"imgs/{fname}"


def _img_tag(src_rel, alt):
    if not src_rel:
        return f"<p style='color:#aaa;font-size:.85em'>{alt} not yet generated — run --step visualize first</p>"
    return f"<img src='{src_rel}' alt='{alt}' style='max-width:100%;border-radius:6px;margin-top:6px'/>"


def _read_csv_safe(path):
    """Read CSV and return list[dict]; returns None if file missing or pandas unavailable."""
    if not _HAS_PANDAS or not path or not os.path.exists(path):
        return None
    try:
        import pandas as pd
        df = pd.read_csv(path)
        return df.to_dict(orient="records")
    except Exception:
        return None


def _pagination_widget(tbl_id: str, ctrl_id: str, n_total: int, page_size: int = 20):
    """Return (controls_html, js_html) for a paginated <table id=tbl_id>.

    Each <tr> must carry a data-row='index' attribute (0-based).
    The table must have a <tbody>.
    """
    n_pages = max(1, (n_total + page_size - 1) // page_size)
    controls = f"""
<div id="{ctrl_id}_bar" style="display:flex;align-items:center;gap:8px;flex-wrap:wrap;margin:8px 0;font-size:.83em">
  <button id="{ctrl_id}_prev"
          style="padding:2px 10px;cursor:pointer;border:1px solid #ccc;background:#fff;border-radius:3px"
          onclick="window.{ctrl_id}_show(window.{ctrl_id}_cur - 1)">&#8592; Prev</button>
  <span id="{ctrl_id}_pages" style="display:inline-flex;flex-wrap:wrap;gap:2px"></span>
  <button id="{ctrl_id}_next"
          style="padding:2px 10px;cursor:pointer;border:1px solid #ccc;background:#fff;border-radius:3px"
          onclick="window.{ctrl_id}_show(window.{ctrl_id}_cur + 1)">Next &#8594;</button>
  <span id="{ctrl_id}_info" style="color:#666;margin-left:6px"></span>
</div>"""
    js = f"""
<script>
(function() {{
  var PAGE  = {page_size};
  var TOTAL = {n_total};
  var PAGES = {n_pages};
  var cur   = 0;
  function show(p) {{
    p = Math.max(0, Math.min(p, PAGES - 1));
    cur = p;
    window.{ctrl_id}_cur = cur;
    var rows = document.querySelectorAll('#{tbl_id} tbody tr');
    var start = cur * PAGE;
    rows.forEach(function(tr) {{
      var i = parseInt(tr.getAttribute('data-row'));
      tr.style.display = (i >= start && i < start + PAGE) ? '' : 'none';
    }});
    var info = document.getElementById('{ctrl_id}_info');
    if (info) info.textContent = 'Showing ' + (start+1) + '\u2013' + Math.min(start+PAGE, TOTAL) + ' of ' + TOTAL;
    var prev = document.getElementById('{ctrl_id}_prev');
    var next = document.getElementById('{ctrl_id}_next');
    if (prev) prev.disabled = (cur === 0);
    if (next) next.disabled = (cur === PAGES - 1);
    document.querySelectorAll('.{ctrl_id}_pg').forEach(function(btn) {{
      var active = parseInt(btn.getAttribute('data-p')) === cur;
      btn.style.fontWeight  = active ? 'bold'    : 'normal';
      btn.style.background  = active ? '#e3f0ff' : '#fff';
      btn.style.borderColor = active ? '#1a6fbb' : '#ccc';
    }});
  }}
  window.{ctrl_id}_show = show;
  window.{ctrl_id}_cur  = 0;
  document.addEventListener('DOMContentLoaded', function() {{
    var pgArea = document.getElementById('{ctrl_id}_pages');
    if (pgArea) {{
      for (var p = 0; p < PAGES; p++) {{
        (function(pp) {{
          var b = document.createElement('button');
          b.textContent = pp + 1;
          b.setAttribute('data-p', pp);
          b.className = '{ctrl_id}_pg';
          b.style.cssText = 'margin:0 1px;padding:2px 7px;cursor:pointer;border:1px solid #ccc;background:#fff;border-radius:3px';
          b.onclick = function() {{ show(pp); }};
          pgArea.appendChild(b);
        }})(p);
      }}
    }}
    show(0);
  }});
}})();
</script>"""
    return controls, js
