"""Interaction and accessibility checks for the portfolio web page (portfolio/web/index.html).

Usage: python3 portfolio/src/web_check.py [path-to-chrome]
Wraps the page in a minimal document, appends the checks in web_interaction_test.html.part,
runs headless Chrome at 1280, 820 and 500 px, and prints PASS/FAIL per check.
Exits non-zero if any check fails or no results are produced.
"""
import html, re, subprocess, sys, tempfile
from pathlib import Path

P = Path(__file__).resolve().parents[1]
CHROME = sys.argv[1] if len(sys.argv) > 1 else "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
HEAD = ('<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">'
        '<style>:root{padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px)}body{margin:0;font-size:14px}[hidden]{display:none!important}</style></head><body>')
page = HEAD + (P / "web" / "index.html").read_text(encoding="utf-8") + (P / "src" / "web_interaction_test.html.part").read_text(encoding="utf-8") + "</body></html>"
failed = False
with tempfile.TemporaryDirectory() as d:
    f = Path(d) / "test.html"
    f.write_text(page, encoding="utf-8")
    for w in (1280, 820, 500):
        out = subprocess.run([CHROME, "--headless=new", "--disable-gpu", "--no-first-run", f"--user-data-dir={d}/c{w}", f"--window-size={w},900",
                              "--virtual-time-budget=6000", "--dump-dom", f.as_uri()], capture_output=True, text=True, timeout=120).stdout
        m = re.search(r'<pre id="qa-results"[^>]*>(.*?)</pre>', out, re.S)
        res = html.unescape(m.group(1)) if m else ""
        npass, nfail = res.count("PASS "), res.count("FAIL ")
        print(f"width {w}px: {npass} PASS, {nfail} FAIL")
        for line in res.splitlines():
            if line.startswith("FAIL") or line.startswith("errors"):
                print("  " + line)
        failed |= (not res) or nfail > 0 or "errors: []" not in res
sys.exit(1 if failed else 0)
