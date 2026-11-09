"""Build the public website (GitHub Pages) from the portfolio web page source.

Usage: python3 portfolio/src/build_site.py [out_dir]   (default: site/)
portfolio/web/index.html is a page body (title, styles, markup, script) without a document skeleton.
This wraps it in a complete HTML document and copies the slide deck PDF next to it. Output is not committed;
.github/workflows/pages.yml builds and deploys it.
"""
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "site"
DESC = ("What cancer pathway scores really measure in bulk tumours: a tested audit tool, a replicated stromal reading of EMT "
        "confirmed in single cells, and a rigorous recurrence-prediction test in 452 colon cancers.")
HEAD = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">\n'
        f'<meta name="description" content="{DESC}">\n<meta name="author" content="Prasad Akula">\n'
        '<style>:root{box-sizing:border-box;padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px)}'
        'html{scroll-padding-top:env(safe-area-inset-top,0px)}body{margin:0;padding:0}img{max-width:100%}'
        '[hidden]:not([hidden=until-found i]){display:none!important}</style>\n</head>\n<body>\n')


def main():
    body = (ROOT / "portfolio" / "web" / "index.html").read_text(encoding="utf-8")
    if "<html" in body.lower()[:200]:
        raise SystemExit("web source already has a document skeleton; expected a page body")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "index.html").write_text(HEAD + body.rstrip() + "\n</body>\n</html>\n", encoding="utf-8")
    shutil.copy2(ROOT / "portfolio" / "presentation" / "GSE39582_recurrence_portfolio.pdf", OUT / "GSE39582_recurrence_portfolio.pdf")
    (OUT / ".nojekyll").write_text("")
    print("built", OUT, sorted(p.name for p in OUT.iterdir()))


if __name__ == "__main__":
    main()
