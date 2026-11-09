"""Build the slide visuals for the 21-slide presentation (SVG with editable text, plus PNG).

Usage: python portfolio/src/make_slide_visuals.py
Writes portfolio/diagrams/slides/V*.{svg,png}. All numbers come from values.py (frozen outputs,
read-only). Biological mechanisms are textbook background and are labelled as such; no visual
shows measured pathway activation, mutation status, protein function or cell composition.
"""
import math
import textwrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle, Ellipse, FancyBboxPatch, Rectangle  # noqa: E402

import values  # noqa: E402

OUT = values.ROOT / "portfolio" / "diagrams" / "slides"
W, H = 12.2, 6.0
INK, MUTED, FAINT, GRID, PANEL, WHITE = "#1f2933", "#52606d", "#8a94a3", "#dde3ea", "#f4f6f9", "#ffffff"
PRIMARY, S7, FIVE, APPARENT = "#2a78d6", "#eb6834", "#1baf7a", "#9aa5b1"
PW = {"EMT": "#c2366c", "E2F targets": "#6650c4", "Wnt/β-catenin": "#00907e", "p53 pathway": "#a86b12"}
TAGS = {"est": ("ESTABLISHED MECHANISM", "#52606d", WHITE, "-"), "lit": ("CRC LITERATURE", "#7b5ea7", WHITE, "-"),
        "meas": ("MEASURED IN THIS STUDY", INK, WHITE, "-"), "nt": ("NOT TESTED / HYPOTHESIS", WHITE, MUTED, "--"),
        "concept": ("CONCEPTUAL ILLUSTRATION", WHITE, MUTED, "--")}
CELL = {"epi": "#e7a6bd", "epin": "#4b3a8c", "fib": "#d4ad86", "fibn": "#6b4a2b", "imm": "#5d6db8", "endo": "#9a7fc7", "rbc": "#d9534f", "matrix": "#f7e9ef"}

plt.rcParams.update({"font.family": ["Arial", "DejaVu Sans"], "font.size": 12, "svg.fonttype": "none",
                     "svg.hashsalt": "gse39582-slides", "text.color": INK, "axes.edgecolor": FAINT,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.labelcolor": MUTED})


def canvas(w=W, h=H):
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.axis("off")
    return fig, ax


def wrap(t, w, fs):
    return "\n".join(textwrap.fill(p, max(8, int(w * 72 / (fs * 0.52)))) for p in t.split("\n"))


def box(ax, x, y, w, h, text=None, fc=WHITE, ec=GRID, fs=12, bold=False, color=INK, ls="-", lw=1.3, r=0.08, ha="center", pad=0.1):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw, ls=ls))
    if text:
        tx = x + w / 2 if ha == "center" else x + pad
        ax.text(tx, y + h / 2, wrap(text, w - 2 * pad, fs), ha=ha, va="center", fontsize=fs, fontweight="bold" if bold else "normal",
                color=color, linespacing=1.2)


def act(ax, x1, y1, x2, y2, color=INK, lw=1.8, rad=0.0, label=None, lfs=10.5, loff=(0, -0.12)):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, mutation_scale=16,
                connectionstyle=f"arc3,rad={rad}", shrinkA=2, shrinkB=2))
    if label:
        ax.text((x1 + x2) / 2 + loff[0], (y1 + y2) / 2 + loff[1], label, ha="center", va="bottom", fontsize=lfs, color=color, style="italic")


def inh(ax, x1, y1, x2, y2, color="#b42318", lw=1.8, rad=0.0, label=None, lfs=10.5, loff=(0, -0.12)):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="-[, widthB=0.45, lengthB=0", color=color, lw=lw,
                connectionstyle=f"arc3,rad={rad}", shrinkA=2, shrinkB=3))
    if label:
        ax.text((x1 + x2) / 2 + loff[0], (y1 + y2) / 2 + loff[1], label, ha="center", va="bottom", fontsize=lfs, color=color, style="italic")


def tag(ax, x, y, kind, fs=9):
    label, fc, tc, ls = TAGS[kind]
    t = ax.text(x, y, label, fontsize=fs, fontweight="bold", color=tc, va="center", ha="left",
                bbox=dict(boxstyle="round,pad=0.28,rounding_size=0.12", fc=fc, ec=MUTED if fc == WHITE else fc, ls=ls, lw=1))
    return t


def evidence(ax, x, y, w, items, fs=11.2, gap=0.14):
    """Stacked evidence panel: items = [(kind, text), ...]; returns bottom y."""
    for kind, text in items:
        tag(ax, x, y + 0.13, kind, fs=8.6)
        t = wrap(text, w, fs)
        n = t.count("\n") + 1
        ax.text(x, y + 0.36, t, fontsize=fs, va="top", linespacing=1.22, color=INK)
        y += 0.36 + n * fs * 1.22 / 72 + gap
    return y


def legend_arrows(ax, x, y):
    act(ax, x, y, x + 0.45, y)
    ax.text(x + 0.55, y, "activates / leads to", fontsize=10, va="center", color=MUTED)
    inh(ax, x + 2.15, y, x + 2.6, y)
    ax.text(x + 2.7, y, "inhibits / represses", fontsize=10, va="center", color=MUTED)


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.svg", metadata={"Date": None})
    fig.savefig(OUT / f"{name}.png", dpi=200, metadata={"Software": None})
    plt.close(fig)


def fmt_hr(t):
    return f"{t[0]:.3f} ({t[1]:.3f}–{t[2]:.3f})"


# ------------------------------------------------------------------ V01 cover illustration
def v01(V):
    fig, ax = canvas(6.0, 6.0)
    ax.add_patch(FancyBboxPatch((0.15, 0.15), 5.7, 5.7, boxstyle="round,pad=0,rounding_size=0.25", fc=CELL["matrix"], ec="none"))
    # normal crypt (left): a test-tube of epithelium with lumen
    def crypt(cx, top, depth, n, cell_fc, wobble=0.0, seed=1):
        s = seed
        def rnd():
            nonlocal s
            s = (s * 16807) % 2147483647
            return s / 2147483647
        for side in (-1, 1):
            for k in range(n):
                yy = top + k * depth / n
                xx = cx + side * (0.42 + wobble * (rnd() - 0.5))
                ax.add_patch(FancyBboxPatch((xx - 0.13, yy), 0.26, depth / n * 0.92, boxstyle="round,pad=0,rounding_size=0.05",
                                            fc=cell_fc, ec=WHITE, lw=0.8))
                ax.add_patch(Ellipse((xx + side * 0.04, yy + depth / n * 0.62), 0.11, 0.13, fc=CELL["epin"], ec="none"))
        for k in range(5):
            a = math.pi * (k + 0.5) / 5
            xx, yy = cx + 0.42 * math.cos(a) * -1, top + depth + 0.42 * math.sin(a)
            ax.add_patch(Circle((xx, yy), 0.13, fc=cell_fc, ec=WHITE, lw=0.8))
            ax.add_patch(Circle((xx, yy), 0.055, fc=CELL["epin"], ec="none"))
    crypt(1.25, 1.1, 2.9, 9, CELL["epi"], 0.0, 3)
    ax.text(1.25, 0.75, "normal crypt", ha="center", fontsize=11, color=MUTED)
    # tumour glands (right): crowded, irregular
    import random
    rng = random.Random(39582)
    for (gx, gy, r) in [(3.5, 2.0, 0.75), (4.75, 3.0, 0.62), (3.6, 3.85, 0.66), (4.6, 1.35, 0.5)]:
        ax.add_patch(Ellipse((gx, gy), r * 0.9, r * 0.65, fc=WHITE, ec="none", alpha=0.9))
        n = int(r * 22)
        for k in range(n):
            a = 2 * math.pi * k / n
            rr = r * (0.62 + 0.12 * rng.random())
            xx, yy = gx + math.cos(a) * rr, gy + math.sin(a) * rr * 0.78
            ax.add_patch(Ellipse((xx, yy), 0.2, 0.3, angle=math.degrees(a) + 90, fc=CELL["epi"], ec=WHITE, lw=0.6))
            ax.add_patch(Circle((xx + 0.02 * rng.random(), yy), 0.05 + 0.02 * rng.random(), fc=CELL["epin"], ec="none"))
    for _ in range(26):
        x, y = 2.6 + rng.random() * 3.0, 0.5 + rng.random() * 4.9
        ax.add_patch(Ellipse((x, y), 0.38, 0.08, angle=rng.uniform(-35, 35), fc=CELL["fib"], ec="none"))
    for _ in range(22):
        x, y = 2.6 + rng.random() * 3.0, 0.5 + rng.random() * 4.9
        ax.add_patch(Circle((x, y), 0.06, fc=CELL["imm"], ec="none"))
    ax.text(4.2, 5.35, "tumour glands in stroma", ha="center", fontsize=11, color=MUTED)
    ax.text(3.0, 5.72, "Conceptual illustration, not patient data", ha="center", fontsize=9.5, color=MUTED, style="italic")
    save(fig, "V01_cover_illustration")


# ------------------------------------------------------------------ V02 problem
def v02(V):
    fig, ax = canvas()
    steps = [("The clinical problem", "After surgery for stage II–III colon cancer, some patients recur and others do not. Stage alone sorts that risk imperfectly."),
             ("The question", "Do four predefined tumour pathway-expression scores add recurrence-discrimination information beyond clinical variables?"),
             ("The trap", "A larger model usually fits its own patients better. That optimism makes signatures look more predictive than they are."),
             ("The approach", "Fix the features before outcomes, compare clinical vs clinical + pathways, and measure the gain after removing optimism, with its uncertainty.")]
    bw, gap = 2.72, 0.33
    for i, (h, t) in enumerate(steps):
        x = 0.15 + i * (bw + gap)
        box(ax, x, 0.35, bw, 2.9, fc=PANEL if i != 3 else "#e8f1fc", ec=GRID if i != 3 else PRIMARY)
        ax.text(x + 0.2, 0.7, f"{i + 1}", fontsize=26, fontweight="bold", color=PRIMARY if i == 3 else FAINT, va="center")
        ax.text(x + 0.2, 1.25, h, fontsize=15, fontweight="bold", va="center")
        ax.text(x + 0.2, 1.6, wrap(t, bw - 0.4, 12.5), fontsize=12.5, va="top", linespacing=1.3)
        if i < 3:
            act(ax, x + bw + 0.03, 2.0, x + bw + gap - 0.03, 2.0, color=FAINT)
    c = V["cohort"]
    facts = [(f"{c['n']}", "patients, stage II–III"), (f"{c['events']}", "recurrences"), ("4", "Hallmark pathway scores"),
             (f"{V['probes']['total']:,}", "probe sets, GEO GSE39582")]
    for i, (n_, l) in enumerate(facts):
        x = 0.15 + i * (bw + gap)
        ax.text(x + 0.2, 4.35, n_, fontsize=28, fontweight="bold", color=INK, va="center")
        ax.text(x + 0.2, 4.85, l, fontsize=12.5, color=MUTED, va="center")
    ax.plot([0.15, W - 0.15], [3.95, 3.95], color=GRID, lw=1)
    ax.text(0.15, 5.55, "Endpoint: cause-specific hazard of recurrence (not absolute risk). Validation: internal only.", fontsize=12, color=MUTED)
    save(fig, "V02_problem_question")


# ------------------------------------------------------------------ V03 relationship map
def v03(V):
    fig, ax = canvas()
    # compartments
    box(ax, 0.15, 0.25, 8.9, 5.45, fc="#fbfaf7", ec=GRID, r=0.3)
    ax.text(0.35, 0.5, "Tumour epithelial cell", fontsize=11, color=MUTED, fontweight="bold")
    ax.add_patch(Ellipse((5.4, 3.85), 6.6, 2.5, fc="#f1eefb", ec="#c9c0ee", lw=1.3))
    ax.text(5.4, 4.95, "nucleus", fontsize=11, color="#6650c4", style="italic", ha="center", va="center")
    # Wnt axis
    box(ax, 0.35, 0.8, 1.35, 0.55, "Wnt signal", fc=WHITE, ec=PW["Wnt/β-catenin"], fs=11.5, bold=True)
    box(ax, 0.35, 1.75, 1.75, 0.75, "APC · AXIN · GSK3β destruction complex", fc=WHITE, ec=PW["Wnt/β-catenin"], fs=10.5)
    box(ax, 2.6, 1.85, 1.2, 0.55, "β-catenin", fc=PW["Wnt/β-catenin"], ec=PW["Wnt/β-catenin"], fs=11.5, bold=True, color=WHITE)
    inh(ax, 1.0, 1.35, 1.0, 1.75, label=None)
    inh(ax, 2.1, 2.12, 2.6, 2.12)
    box(ax, 2.55, 3.0, 1.5, 0.55, "TCF/LEF", fc=WHITE, ec=PW["Wnt/β-catenin"], fs=11.5)
    act(ax, 3.2, 2.4, 3.25, 3.0)
    box(ax, 2.35, 4.0, 1.9, 0.6, "MYC, cyclin D (CCND1)", fc=WHITE, ec=PW["Wnt/β-catenin"], fs=10.5)
    act(ax, 3.3, 3.55, 3.3, 4.0)
    # E2F axis
    box(ax, 4.85, 4.0, 1.55, 0.6, "cyclin D–CDK4/6", fc=WHITE, ec=PW["E2F targets"], fs=10.5)
    act(ax, 4.25, 4.3, 4.85, 4.3)
    box(ax, 5.6, 2.9, 1.35, 0.55, "RB", fc=WHITE, ec=PW["E2F targets"], fs=11.5, bold=True)
    inh(ax, 6.2, 4.0, 6.2, 3.45, label="phosphorylates", lfs=9.5, loff=(-0.62, 0.05))
    box(ax, 7.45, 2.9, 1.25, 0.55, "E2F", fc=PW["E2F targets"], ec=PW["E2F targets"], fs=11.5, bold=True, color=WHITE)
    inh(ax, 6.95, 3.17, 7.45, 3.17)
    box(ax, 6.6, 4.0, 2.25, 0.6, "G1/S genes: cyclin E, MCM, PCNA", fc=WHITE, ec=PW["E2F targets"], fs=10)
    act(ax, 8.05, 3.45, 7.85, 4.0)
    box(ax, 6.75, 5.0, 2.0, 0.5, "proliferation", fc=PANEL, ec=PW["E2F targets"], fs=11.5, bold=True)
    act(ax, 7.7, 4.6, 7.75, 5.0)
    # p53 axis
    box(ax, 6.0, 0.8, 1.65, 0.55, "DNA damage", fc=WHITE, ec=PW["p53 pathway"], fs=11.5)
    box(ax, 6.2, 1.75, 1.2, 0.55, "p53", fc=PW["p53 pathway"], ec=PW["p53 pathway"], fs=11.5, bold=True, color=WHITE)
    act(ax, 6.8, 1.35, 6.8, 1.75)
    box(ax, 8.0, 1.0, 0.9, 0.5, "MDM2", fc=WHITE, ec=PW["p53 pathway"], fs=11)
    act(ax, 7.4, 1.9, 8.0, 1.35, color=PW["p53 pathway"], rad=-0.2)
    inh(ax, 8.1, 1.5, 7.4, 2.1, rad=-0.2)
    box(ax, 4.3, 2.9, 1.0, 0.5, "p21", fc=WHITE, ec=PW["p53 pathway"], fs=11)
    act(ax, 6.2, 2.15, 4.95, 2.9, color=PW["p53 pathway"])
    inh(ax, 4.85, 3.4, 5.25, 4.0, label="blocks CDK", lfs=9.5, loff=(-0.55, 0.0))
    # EMT axis (right panel, outside the cell)
    box(ax, 9.3, 0.25, 2.75, 5.45, fc="#fbf5f8", ec=GRID, r=0.3)
    ax.text(9.5, 0.55, "EMT and tumour stroma", fontsize=11, color=MUTED, fontweight="bold")
    box(ax, 9.55, 0.85, 2.25, 0.5, "TGF-β", fc=WHITE, ec=PW["EMT"], fs=11.5)
    box(ax, 9.55, 1.75, 2.25, 0.55, "SNAIL · ZEB · TWIST", fc=PW["EMT"], ec=PW["EMT"], fs=11, bold=True, color=WHITE)
    act(ax, 10.67, 1.35, 10.67, 1.75)
    box(ax, 9.55, 2.75, 1.05, 0.62, "E-cadherin (CDH1)", fc=WHITE, ec=PW["EMT"], fs=9.5)
    box(ax, 10.75, 2.75, 1.05, 0.62, "vimentin, fibronectin", fc=WHITE, ec=PW["EMT"], fs=9.5)
    inh(ax, 10.2, 2.3, 10.07, 2.75)
    act(ax, 11.15, 2.3, 11.27, 2.75)
    box(ax, 9.55, 3.8, 2.25, 0.6, "mesenchymal, migratory state", fc=PANEL, ec=PW["EMT"], fs=10.5, bold=True)
    act(ax, 11.27, 3.37, 10.9, 3.8)
    box(ax, 9.55, 4.75, 2.25, 0.7, "fibroblasts (CAFs) also express mesenchymal genes", fc=WHITE, ec=CELL["fib"], fs=9.8, ls="--")
    legend_arrows(ax, 0.4, 5.35)
    tag(ax, 4.6, 0.52, "est", fs=8.6)
    save(fig, "V03_biology_map")


# ------------------------------------------------------------------ pathway mechanism slides
def pathway_slide(V, name, draw_mech, items):
    fig, ax = canvas()
    box(ax, 0.1, 0.1, 6.9, 5.8, fc="#fbfaf7", ec=GRID, r=0.25)
    tag(ax, 0.3, 0.38, "est", fs=8.6)
    draw_mech(ax)
    evidence(ax, 7.3, 0.2, 4.8, items)
    save(fig, name)


def mech_wnt(ax):
    col = PW["Wnt/β-catenin"]
    for (y0, title) in ((0.75, "APC intact, no Wnt signal"), (3.25, "APC lost or CTNNB1 mutated (most colorectal cancers)")):
        ax.text(0.35, y0, title, fontsize=12, fontweight="bold", va="center")
    # state 1
    box(ax, 0.35, 1.05, 2.4, 0.9, "destruction complex\nAPC · AXIN · GSK3β · CK1", fc=WHITE, ec=col, fs=10.5)
    ax.text(3.0, 1.0, "phosphorylates", fontsize=9.5, color=INK, style="italic", ha="center")
    box(ax, 3.25, 1.2, 1.25, 0.6, "β-catenin", fc=col, ec=col, fs=11, bold=True, color=WHITE)
    act(ax, 2.75, 1.5, 3.25, 1.5)
    box(ax, 5.0, 1.2, 1.7, 0.6, "proteasomal degradation", fc=PANEL, ec=GRID, fs=10)
    act(ax, 4.5, 1.5, 5.0, 1.5)
    ax.text(0.35, 2.35, "Target genes stay off", fontsize=11, color=MUTED, style="italic")
    ax.plot([0.35, 6.75], [2.8, 2.8], color=GRID, lw=1)
    # state 2 with nucleus
    box(ax, 0.35, 3.55, 2.4, 0.9, "complex disabled", fc=WHITE, ec=col, fs=10.5, ls="--")
    box(ax, 3.0, 3.7, 1.25, 0.6, "β-catenin", fc=col, ec=col, fs=11, bold=True, color=WHITE)
    ax.text(2.6, 4.75, "β-catenin escapes degradation", fontsize=9.5, color=MUTED, style="italic", ha="center")
    ax.text(3.62, 4.55, "accumulates", fontsize=9.5, ha="center", color=MUTED, style="italic")
    ax.add_patch(Ellipse((5.55, 4.6), 2.3, 1.95, fc="#f1eefb", ec="#c9c0ee", lw=1.3))
    ax.text(5.55, 3.82, "nucleus", fontsize=10, color="#6650c4", style="italic", ha="center")
    box(ax, 4.75, 4.05, 1.6, 0.5, "β-catenin + TCF/LEF", fc=WHITE, ec=col, fs=10)
    act(ax, 4.25, 4.0, 4.75, 4.25)
    box(ax, 4.6, 4.85, 1.9, 0.6, "MYC, AXIN2, cyclin D: transcribed", fc=WHITE, ec=col, fs=9.8)
    act(ax, 5.55, 4.55, 5.55, 4.85)


def mech_e2f(ax):
    col = PW["E2F targets"]
    ax.add_patch(Ellipse((3.55, 3.9), 6.5, 3.2, fc="#f1eefb", ec="#c9c0ee", lw=1.3))
    ax.text(0.35, 0.8, "G1 → S transition", fontsize=12, fontweight="bold", va="center")
    box(ax, 0.35, 1.15, 1.6, 0.6, "growth signals (e.g. MYC)", fc=WHITE, ec=GRID, fs=10)
    box(ax, 2.35, 1.15, 1.85, 0.6, "cyclin D–CDK4/6", fc=WHITE, ec=col, fs=10.5, bold=True)
    act(ax, 1.95, 1.45, 2.35, 1.45)
    # RB-E2F complex
    box(ax, 0.6, 2.6, 2.6, 0.75, "RB bound to E2F: E2F target genes repressed", fc=PANEL, ec=col, fs=10)
    inh(ax, 3.27, 1.75, 2.6, 2.6, label="phosphorylates RB", lfs=9.5, loff=(0.75, 0.05))
    box(ax, 3.85, 2.65, 1.1, 0.65, "RB-P (inactive)", fc=WHITE, ec=GRID, fs=9.8)
    box(ax, 5.3, 2.65, 1.35, 0.65, "E2F free", fc=col, ec=col, fs=11, bold=True, color=WHITE)
    act(ax, 3.2, 2.97, 3.85, 2.97)
    act(ax, 4.95, 2.97, 5.3, 2.97)
    ax.text(0.55, 4.3, "nucleus", fontsize=10, color="#6650c4", style="italic")
    box(ax, 2.9, 4.1, 2.9, 0.65, "transcribes cyclin E, MCM2–7, PCNA, CDC6", fc=WHITE, ec=col, fs=10)
    act(ax, 5.97, 3.3, 5.5, 4.1)
    box(ax, 2.3, 5.2, 2.2, 0.5, "DNA replication", fc=PANEL, ec=col, fs=10.5, bold=True)
    box(ax, 4.75, 5.2, 1.85, 0.5, "proliferation", fc=PANEL, ec=col, fs=10.5, bold=True)
    act(ax, 3.8, 4.75, 3.4, 5.2)
    act(ax, 5.0, 4.75, 5.65, 5.2)
    box(ax, 4.75, 1.15, 1.95, 0.6, "p21 (from p53)", fc=WHITE, ec=PW["p53 pathway"], fs=10)
    inh(ax, 4.75, 1.45, 4.2, 1.45)


def mech_p53(ax):
    col = PW["p53 pathway"]
    ax.text(0.35, 0.8, "Stress response", fontsize=12, fontweight="bold", va="center")
    box(ax, 0.35, 1.2, 1.9, 0.6, "DNA damage, oncogenic stress", fc=WHITE, ec=GRID, fs=10)
    box(ax, 2.7, 1.2, 1.3, 0.6, "ATM / ATR", fc=WHITE, ec=col, fs=10.5)
    act(ax, 2.25, 1.5, 2.7, 1.5)
    box(ax, 2.6, 2.4, 1.5, 0.7, "p53 stabilised", fc=col, ec=col, fs=11, bold=True, color=WHITE)
    act(ax, 3.35, 1.8, 3.35, 2.4)
    box(ax, 5.0, 2.45, 1.3, 0.6, "MDM2", fc=WHITE, ec=col, fs=11)
    act(ax, 4.1, 2.62, 5.0, 2.62, color=col, label="induces", lfs=9.5)
    inh(ax, 5.0, 2.92, 4.1, 2.92, label="degrades p53", lfs=9.5, loff=(0, 0.3))
    ax.add_patch(Ellipse((3.6, 4.55), 6.0, 1.75, fc="#f1eefb", ec="#c9c0ee", lw=1.3))
    ax.text(6.55, 3.85, "nucleus: target genes", fontsize=10, color="#6650c4", style="italic", ha="right", va="bottom")
    for x, t, out in ((0.9, "CDKN1A (p21)", "cell-cycle arrest"), (2.85, "BAX, PUMA", "apoptosis"), (4.8, "GADD45", "DNA repair")):
        box(ax, x, 4.2, 1.6, 0.55, t, fc=WHITE, ec=col, fs=10)
        act(ax, 3.35, 3.1, x + 0.8, 4.2, color=col)
        box(ax, x, 5.15, 1.6, 0.5, out, fc=PANEL, ec=GRID, fs=10, bold=True)
        act(ax, x + 0.8, 4.75, x + 0.8, 5.15)
    box(ax, 4.6, 1.2, 2.2, 0.6, "TP53 missense mutation impairs target activation", fc=WHITE, ec="#b42318", fs=9.2, ls="--")


def mech_emt(ax):
    col = PW["EMT"]
    ax.text(0.35, 0.8, "Two biological sources of the same mesenchymal genes", fontsize=12, fontweight="bold", va="center")
    # possibility 1: tumour cells
    ax.text(0.35, 1.3, "1 · Tumour cells undergoing EMT", fontsize=11, color=col, fontweight="bold")
    for i in range(4):
        ax.add_patch(FancyBboxPatch((0.4 + i * 0.36, 1.6), 0.32, 0.6, boxstyle="round,pad=0,rounding_size=0.04", fc=CELL["epi"], ec=WHITE, lw=1))
    ax.text(1.1, 2.45, "E-cadherin junctions", fontsize=9.5, ha="center", color=MUTED)
    act(ax, 1.95, 1.9, 2.9, 1.9, color=col, label="TGF-β →\nSNAIL/ZEB/TWIST", lfs=9.2, loff=(0, -0.1))
    for k, (dx, dy, ang) in enumerate(((3.55, 1.75, -10), (4.0, 2.1, 15))):
        ax.add_patch(Ellipse((dx, dy), 0.85, 0.2, angle=ang, fc=col, ec=WHITE, lw=1))
    ax.text(3.75, 2.4, "CDH1 down;\nvimentin, fibronectin up", fontsize=9.5, ha="center", va="top", color=MUTED)
    # possibility 2: stroma
    ax.text(0.35, 3.0, "2 · Fibroblasts in the tumour stroma", fontsize=11, color=col, fontweight="bold")
    import random
    rng = random.Random(7)
    for _ in range(9):
        x, y = 0.5 + rng.random() * 3.6, 3.35 + rng.random() * 0.75
        ax.add_patch(Ellipse((x, y), 0.6, 0.13, angle=rng.uniform(-25, 25), fc=col, ec="none", alpha=0.85))
    ax.text(0.4, 4.3, "CAFs express collagens, vimentin and fibronectin\nwithout any change in the tumour cells", fontsize=9.5, ha="left", va="top", color=MUTED)
    # both feed the score
    box(ax, 4.75, 2.6, 2.0, 1.1, "one bulk EMT score (198 genes)", fc=WHITE, ec=col, fs=11, bold=True, lw=2)
    act(ax, 4.55, 2.1, 4.95, 2.6, color=col)
    act(ax, 4.2, 3.7, 4.95, 3.55, color=col)
    box(ax, 0.35, 4.8, 6.4, 0.85, "An elevated bulk EMT-associated score does not establish tumour-cell EMT.", fc="#fbe9f0", ec=col, fs=12.5, bold=True)


def v04(V):
    p = next(x for x in V["pathways"] if x["label"] == "Wnt/β-catenin")
    items = [("lit", "Wnt-pathway genes are genomically altered in 93 % of colorectal cancers; APC or CTNNB1 in about 80 % (TCGA 2012). Marked WNT/MYC activation characterises CMS2 (Guinney 2015). Genomic alteration is not measured activation."),
             ("meas", f"Mean standardised mRNA of {p['genes']} target genes. Primary HR {fmt_hr(p['hr'])}; S-7 {fmt_hr(p['hr_s7'])}. Within-set coherence {p['coherence']:.3f} (loosely co-expressed). {p['lobo_sign_changes']} sign change when single batches are omitted. Post-hoc (same data): lower in BRAF-mutant ({V['X']['braf_wnt_diff']:+.2f} SD) and CIMP+ ({V['X']['cimp_wnt_diff']:+.2f} SD) tumours, one serrated axis."),
             ("nt", "Pathway activation (no nuclear β-catenin or APC data). No explanation is offered for the near-zero estimate.")]
    pathway_slide(V, "V04_wnt", mech_wnt, items)


def v05(V):
    p = next(x for x in V["pathways"] if x["label"] == "E2F targets")
    f5 = V["five_set_e2f"]
    items = [("lit", "Evidence on proliferation and outcome is mixed: lower proliferation-signature expression in relapsed primaries (Anjomshoaa 2009), but high Ki-67 associated with poorer survival in a meta-analysis (Luo 2019)."),
             ("meas", f"Mean standardised mRNA of {p['genes']} E2F target genes; coherence {p['coherence']:.3f}. Primary HR {fmt_hr(p['hr'])}; S-7 {fmt_hr(p['hr_s7'])}. Five-set model {fmt_hr(f5)}: sensitivity only, coincided with added collinearity, not promoted."),
             ("nt", "Which cells are dividing, and whether proliferation causes or protects against recurrence. The direction is not established biology.")]
    pathway_slide(V, "V05_e2f", mech_e2f, items)


def v06(V):
    p = next(x for x in V["pathways"] if x["label"] == "p53 pathway")
    items = [("lit", "TP53 is mutated in about 60 % of non-hypermutated and about 20 % of hypermutated colorectal cancers (TCGA 2012)."),
             ("meas", f"Mean standardised mRNA of {p['genes']} p53 target genes. Primary HR {fmt_hr(p['hr'])}; S-7 {fmt_hr(p['hr_s7'])}. Strongest MMR association of the four scores (ε² {p['mmr_e2']:.4f}, exploratory); {p['lobo_sign_changes']} sign changes under batch omission. Post-hoc: lower in TP53-mutant tumours ({V['X']['tp53_p53_diff']:+.2f} SD; also within pMMR), replicated in GSE41258 and TCGA-COAD (rank-biserial {V['E']['r1_41258']:+.2f}, {V['E']['r1_tcga']:+.2f})."),
             ("nt", "p53 protein function and activity, and the mechanism of the MMR association. A mutation–expression link is not a function assay.")]
    pathway_slide(V, "V06_p53", mech_p53, items)


def v07(V):
    p = next(x for x in V["pathways"] if x["label"] == "EMT")
    items = [("lit", "The mesenchymal CMS4 subtype had worse relapse-free survival (Guinney 2015; GSE39582 contributed, so not independent). Mesenchymal-subtype transcripts in bulk tumours are mostly stromal (Calon 2015; Isella 2015)."),
             ("meas", f"Mean standardised mRNA of {p['genes']} genes. Primary HR {fmt_hr(p['hr'])}, above 1 in every fit but never distinguishable from 1; S-7 {fmt_hr(p['hr_s7'])}. MMR association ε² {p['mmr_e2']:.4f} (exploratory). Post-hoc: tracks fibroblast markers with shared genes removed (ρ {min(v[0] for v in V['E']['disjoint'].values()):.2f}–{max(v[0] for v in V['E']['disjoint'].values()):.2f}) in GSE39582 and three independent cohorts."),
             ("nt", "The cellular source. Composition was estimated computationally, not measured; bulk data cannot separate tumour-cell EMT from stroma.")]
    pathway_slide(V, "V07_emt", mech_emt, items)


# ------------------------------------------------------------------ V08 bulk tissue
def v08(V):
    fig, ax = canvas()
    ax.add_patch(FancyBboxPatch((0.1, 0.1), 6.6, 5.8, boxstyle="round,pad=0,rounding_size=0.25", fc=CELL["matrix"], ec="none"))
    import random
    rng = random.Random(20250)
    glands = [(1.6, 1.6, 0.8), (4.3, 1.4, 0.65), (2.6, 3.9, 0.75), (5.2, 3.9, 0.7)]
    for gx, gy, r in glands:
        ax.add_patch(Ellipse((gx, gy), r * 0.95, r * 0.7, fc=WHITE, ec="none"))
        n = int(r * 24)
        for k in range(n):
            a = 2 * math.pi * k / n
            xx, yy = gx + math.cos(a) * r * 0.72, gy + math.sin(a) * r * 0.55
            ax.add_patch(Ellipse((xx, yy), 0.17, 0.3, angle=math.degrees(a) + 90, fc=CELL["epi"], ec=WHITE, lw=0.6))
            ax.add_patch(Circle((xx, yy), 0.045, fc=CELL["epin"], ec="none"))
    def free(x, y):
        return all(((x - gx) / (r * 0.95)) ** 2 + ((y - gy) / (r * 0.8)) ** 2 > 1 for gx, gy, r in glands) and all(((x - vx) / 0.6) ** 2 + ((y - vy) / 0.4) ** 2 > 1 for vx, vy in ((0.9, 5.0), (6.0, 0.8)))
    for vx, vy in ((0.9, 5.0), (6.0, 0.8)):
        ax.add_patch(Ellipse((vx, vy), 0.8, 0.42, fc=WHITE, ec=CELL["endo"], lw=4))
        for k in range(4):
            ax.add_patch(Ellipse((vx - 0.24 + k * 0.16, vy + (0.04 if k % 2 else -0.04)), 0.11, 0.08, fc=CELL["rbc"], ec="none"))
    placed = 0
    while placed < 38:
        x, y = 0.3 + rng.random() * 6.2, 0.3 + rng.random() * 5.4
        if free(x, y):
            ax.add_patch(Ellipse((x, y), 0.42, 0.09, angle=rng.uniform(-35, 35), fc=CELL["fib"], ec="none"))
            placed += 1
    placed = 0
    while placed < 32:
        x, y = 0.3 + rng.random() * 6.2, 0.3 + rng.random() * 5.4
        if free(x, y):
            ax.add_patch(Circle((x, y), 0.07, fc=CELL["imm"], ec="none"))
            placed += 1
    ax.text(0.3, 5.75, "Conceptual illustration · not a GSE39582 sample · cell counts arbitrary", fontsize=9.5, color=MUTED, style="italic")
    # legend with statuses
    rows = [(CELL["epi"], "Tumour epithelial cells", "NOT SEPARATED", "purity estimated computationally (ESTIMATE)"),
            (CELL["fib"], "Fibroblasts and stroma", "NOT MEASURED", "abundance estimated (ESTIMATE, MCP-counter)"),
            (CELL["imm"], "Immune cells", "NOT MEASURED", "abundance estimated (ESTIMATE, MCP-counter)"),
            (CELL["endo"], "Endothelial cells and vessels", "NOT MEASURED", "abundance estimated (MCP-counter)")]
    y = 0.35
    for c, name, st, note in rows:
        ax.add_patch(Circle((7.1, y + 0.2), 0.13, fc=c, ec=MUTED, lw=0.6))
        ax.text(7.35, y + 0.2, name, fontsize=12, fontweight="bold", va="center")
        ax.text(7.35, y + 0.55, f"{st} · {note}", fontsize=10.5, color=MUTED, va="center")
        y += 0.88
    # flow
    yy = 4.05
    for i, (h, t) in enumerate((("Mixed tissue", "all cell types"), ("One array", "pooled mRNA, one value per gene"), ("One score", "mean of member-gene z-scores"))):
        x = 7.0 + i * 1.75
        box(ax, x, yy, 1.55, 0.95, None, fc=WHITE, ec=GRID)
        ax.text(x + 0.12, yy + 0.3, h, fontsize=11.5, fontweight="bold", va="center")
        ax.text(x + 0.12, yy + 0.62, wrap(t, 1.35, 9.5), fontsize=9.5, color=MUTED, va="center")
        if i < 2:
            act(ax, x + 1.55, yy + 0.47, x + 1.75, yy + 0.47, color=FAINT)
    box(ax, 7.0, 5.2, 5.05, 0.65, f"Post-hoc estimates from the same arrays: EMT score ρ {V['X']['emt_fib']:.2f} with estimated fibroblast content. Estimates, not measured fractions.", fc=PANEL, ec=GRID, fs=10.2, ha="left", pad=0.15)
    save(fig, "V08_bulk_tissue")


# ------------------------------------------------------------------ V10 preprocessing + panel
def v10(V):
    P, pn = V["probes"], V["panel"]
    fig, ax = canvas()
    stages = [(P["total"], "Affymetrix probe sets"), (P["single_gene"], "map to one gene"), (P["genes"], "unique genes"), (pn["union4"], "pathway-member genes"), (4, "scores per patient")]
    notes = ["HG-U133 Plus 2.0; deposited RMA + ComBat values", f"{P['unmapped']:,} unmapped and {P['ambiguous']:,} ambiguous probes excluded",
             f"median of probes per gene ({P['multi_probe_genes']:,} genes have several)", "EMT 198 · E2F 195 · Wnt 41 · p53 194", "standardised inside each training fold"]
    maxw, x0 = 4.9, 0.2
    for i, ((n_, lab), note) in enumerate(zip(stages, notes)):
        w = max(2.7, maxw * (math.log10(n_) / math.log10(P["total"])))
        y = 0.25 + i * 1.08
        last = i == 4
        ax.add_patch(FancyBboxPatch((x0 + (maxw - w) / 2, y), w, 0.78, boxstyle="round,pad=0,rounding_size=0.06", fc=PRIMARY if last else "#e8f1fc", ec=PRIMARY, lw=1.3))
        ax.text(x0 + maxw / 2, y + 0.39, f"{n_:,}  {lab}", ha="center", va="center", fontsize=13, fontweight="bold", color=WHITE if last else INK)
        ax.text(x0 + maxw + 0.2, y + 0.39, wrap(note, 2.7, 10.5), fontsize=10.5, color=MUTED, va="center", linespacing=1.2)
    ax.text(x0, 5.75, "Bar widths follow log feature count (last bar enlarged for legibility).", fontsize=9.5, color=MUTED, style="italic")
    ax.plot([8.25, 8.25], [0.3, 5.7], color=GRID, lw=1)
    ax.text(8.5, 0.45, "Outcome-blind panel selection", fontsize=13, fontweight="bold")
    ax.text(8.5, 0.75, wrap("Ten candidate Hallmark sets formed correlation clusters (six on one mesenchymal axis, r 0.60–0.91 with EMT). "
                            "One representative per cluster was chosen without outcome data.", 3.55, 10.5), fontsize=10.5, va="top", color=INK, linespacing=1.25)
    a = fig.add_axes([8.5 / W, 1 - 5.2 / H, 3.4 / W, 2.75 / H])
    rows = [("10 candidate sets", 77, APPARENT), ("5-set (sensitivity)", pn["cond5"], S7), ("4-set (primary)", pn["cond4"], PRIMARY)]
    for i, (l, v, c) in enumerate(rows):
        a.barh(2 - i, v, color=c, height=0.42)
        a.text(1.05, 2 - i + 0.3, l, fontsize=10.5, va="bottom", color=INK)
        a.text(v * 1.08, 2 - i, f"{v:.2f}" if v < 50 else f"{v:.0f}", va="center", fontsize=11)
    a.set_xscale("log")
    a.set_xlim(1, 250)
    a.set_ylim(-0.4, 2.85)
    a.set_yticks([])
    a.set_xticks([1, 10, 100], ["1", "10", "100"])
    a.set_xlabel("condition number (descriptive, log scale)", fontsize=10)
    for sp in ("top", "right", "left"):
        a.spines[sp].set_visible(False)
    save(fig, "V10_preprocessing_panel")


# ------------------------------------------------------------------ V11 model architecture
def v11(V):
    fig, ax = canvas()
    def model(x, title, sub, blocks, ec):
        box(ax, x, 0.3, 3.7, 3.55, fc=PANEL, ec=ec, lw=1.6)
        ax.text(x + 0.2, 0.65, title, fontsize=15, fontweight="bold", color=ec)
        ax.text(x + 0.2, 0.98, sub, fontsize=10.5, color=MUTED)
        y = 1.35
        for t, fc, ecb in blocks:
            box(ax, x + 0.2, y, 3.3, 1.0, t, fc=fc, ec=ecb, fs=12)
            y += 1.15
    clin = "Clinical block, unpenalised: age (spline), sex, stage II/III, tumour location, adjuvant chemotherapy"
    model(0.15, "M1 · clinical", "Cox model, cause-specific hazard", [(clin, WHITE, GRID)], INK)
    model(4.25, "M2 · primary", f"M1 + four pathway scores (ridge, λ = {V['lambda']['primary']:.3f})",
          [(clin, WHITE, GRID), ("Molecular block, ridge-penalised: EMT · E2F · Wnt · p53 scores", "#e8f1fc", PRIMARY)], PRIMARY)
    model(8.35, "Five-set model", "sensitivity analysis only, never promoted",
          [(clin, WHITE, GRID), ("+ PI3K/AKT/mTOR score; collinear with p53 (r 0.636)", "#e6f7f0", FIVE)], FIVE)
    act(ax, 3.85, 2.4, 4.25, 2.4, color=FAINT)
    for i, t in enumerate(["M1 is exactly M2 with the pathway coefficients set to zero: the comparison isolates the molecular block.",
                           "The clinical block is unpenalised, which keeps the comparison scale-equivariant; shrinking only the molecular block makes the gain in C conservative.",
                           "λ is tuned by five-fold, event-stratified cross-validation inside each training set; M1 is never re-tuned."]):
        ax.text(0.15, 4.2 + i * 0.55, "•  " + wrap(t, 11.6, 12.5), fontsize=12.5, va="top", color=INK)
    save(fig, "V11_model_architecture")


# ------------------------------------------------------------------ V12 nested bootstrap
def v12(V):
    B = V["B"]
    fig, ax = canvas()
    box(ax, 0.15, 0.3, 2.2, 1.2, f"Cohort\n{V['cohort']['n']} patients\n{V['cohort']['events']} recurrences", fc=WHITE, ec=GRID, fs=11.5, bold=True)
    act(ax, 2.35, 0.9, 2.85, 0.9)
    box(ax, 2.85, 0.15, 6.4, 4.25, fc="#f4f8fe", ec=PRIMARY, lw=1.6)
    ax.text(3.05, 0.48, f"Outer loop · {B['out']} patient-level, event-stratified bootstrap resamples", fontsize=12.5, fontweight="bold", color=PRIMARY)
    box(ax, 3.05, 0.8, 2.6, 0.85, "refit M1 and M2 on the resample (scores standardised in-fold)", fc=WHITE, ec=PRIMARY, fs=10.2)
    box(ax, 6.4, 0.8, 2.6, 0.85, "apparent ΔC on the resample", fc=WHITE, ec=PRIMARY, fs=10.5)
    act(ax, 5.65, 1.22, 6.4, 1.22)
    box(ax, 3.05, 1.95, 5.95, 1.45, fc=WHITE, ec="#7aa8e6", ls="--")
    ax.text(3.2, 2.22, f"Inner loop · {B['in']} bootstrap resamples of the outer resample", fontsize=11.5, fontweight="bold", color="#2a5fa8")
    ax.text(3.2, 2.55, wrap("Refit M1 and M2 on each inner resample (paired: the same resample for both). Optimism = C on the inner resample − C on the outer resample. ΔOptimism(k) = mean optimism of M2 − mean optimism of M1.", 5.6, 10.8), fontsize=10.8, va="top")
    act(ax, 7.7, 1.65, 7.7, 1.95)
    box(ax, 3.05, 3.6, 5.95, 0.6, "corrected ΔC(k) = apparent ΔC(k) − ΔOptimism(k)", fc="#e8f1fc", ec=PRIMARY, fs=11.5, bold=True)
    act(ax, 6.0, 3.4, 6.0, 3.6)
    act(ax, 9.25, 2.2, 9.7, 2.2)
    box(ax, 9.7, 0.3, 2.35, 1.75, f"95 % interval = 2.5th and 97.5th percentiles of {B['out']} corrected ΔC values", fc=WHITE, ec=PRIMARY, fs=10.5)
    box(ax, 9.7, 2.35, 2.35, 1.95, "repeated in 5 runs that differ only in random seed: Monte Carlo replicates, not independent experiments", fc=WHITE, ec=GRID, fs=10.3)
    mc = V["mc"]
    ax.text(0.15, 4.75, wrap(f"Precision rule (pre-specified): largest endpoint Monte Carlo SE ≤ {mc['threshold']*100:.0f} % of interval width; "
                        f"achieved {mc['rho_max']*100:.1f} % after topping up from 500 to {B['out']} outer resamples.", 11.9, 11.5), fontsize=11.5, color=INK, va="top")
    ax.text(0.15, 5.2, wrap(f"Point estimate: full-data corrected ΔC = +{V['theta']:.6f}. Reporting was gated on a frozen verifier "
                        f"(PASS {V['verification']['status_counts']['PASS']}/60); a separate re-implementation recomputed the reported quantities from the checkpoints.", 11.9, 11.5), fontsize=11.5, color=INK, va="top", linespacing=1.25)
    save(fig, "V12_nested_bootstrap")


# ------------------------------------------------------------------ V13 central result
def v13(V):
    fig, ax = canvas()
    a = fig.add_axes([0.16, 0.33, 0.31, 0.55])
    rows = [("Clinical (M1)", "M1", INK), ("+ 4 pathways (M2)", "M2_primary", PRIMARY)]
    for i, (l, k, c) in enumerate(rows):
        y = 1 - i
        ap, co = V["C"][k]["apparent"], V["C"][k]["corrected"]
        a.plot([co, ap], [y, y], color=GRID, lw=3, zorder=1)
        a.scatter([ap], [y], s=90, facecolor=WHITE, edgecolor=APPARENT, lw=2, zorder=2)
        a.scatter([co], [y], s=90, color=c, zorder=3)
        a.text(ap + 0.0012, y + 0.16, f"{ap:.4f}", fontsize=11, color=MUTED)
        a.text(co - 0.0012, y + 0.16, f"{co:.4f}", fontsize=11, color=c, ha="right")
    a.set_yticks([1, 0], [r[0] for r in rows], fontsize=12)
    a.set_xlim(0.594, 0.646)
    a.set_ylim(-0.6, 1.6)
    a.set_xlabel("Harrell's C · open = apparent, filled = optimism-corrected", fontsize=10.5)
    for s in ("top", "right", "left"):
        a.spines[s].set_visible(False)
    a.tick_params(axis="y", length=0)
    a.grid(axis="x", color=GRID)
    a.set_axisbelow(True)
    ax.text(0.2, 0.35, "1 · Most of the apparent gain was optimism", fontsize=13.5, fontweight="bold")
    b = fig.add_axes([0.57, 0.33, 0.4, 0.55])
    for r, (l, u) in enumerate(V["runs"]):
        b.plot([l, u], [4 - r, 4 - r], color=PRIMARY, lw=7, solid_capstyle="round")
    b.axvline(0, color=INK, lw=1.2, ls=(0, (3, 3)))
    b.axvline(V["theta"], color=INK, lw=2)
    b.text(V["theta"] + 0.0012, 4.6, r"$\hat{\theta}$" + f" = +{V['theta']:.6f}", fontsize=11)
    b.set_yticks(range(5), [f"Run {r}" for r in (4, 3, 2, 1, 0)], fontsize=11.5)
    b.set_xlim(-0.02, 0.06)
    b.set_ylim(-0.6, 5.0)
    b.set_xlabel("Corrected ΔC, 95 % percentile interval per run", fontsize=10.5)
    for s in ("top", "right", "left"):
        b.spines[s].set_visible(False)
    b.tick_params(axis="y", length=0)
    b.grid(axis="x", color=GRID)
    b.set_axisbelow(True)
    ax.text(6.6, 0.35, "2 · All five corrected intervals include zero", fontsize=13.5, fontweight="bold")
    d_app = V["C"]["M2_primary"]["apparent"] - V["C"]["M1"]["apparent"]
    box(ax, 0.2, 4.55, 11.8, 1.3, None, fc="#e8f1fc", ec=PRIMARY)
    ax.text(0.4, 4.9, f"Apparent ΔC +{d_app:.4f}  →  corrected ΔC +{V['theta']:.6f}  ·  intervals ≈ −0.009 to +0.05  ·  calibration slopes "
                      f"{V['slopes']['M1']:.3f} (M1), {V['slopes']['M2_primary']:.3f} (M2): residual overfitting", fontsize=11.5)
    ax.text(0.4, 5.4, "The study did not demonstrate a reproducible improvement in recurrence discrimination from adding the four pathway scores.",
            fontsize=12.5, fontweight="bold")
    save(fig, "V13_central_result")


# ------------------------------------------------------------------ V14 forest
def v14(V):
    fig, ax = canvas()
    a = fig.add_axes([0.15, 0.2, 0.46, 0.72])
    lo, hi = 0.7, 1.35
    for i, p in enumerate(V["pathways"]):
        y0 = 3 - i
        for dy, (h, l, u), c, m, lab in ((0.14, p["hr"], PRIMARY, "o", "Primary"), (-0.14, p["hr_s7"], S7, "s", "S-7")):
            a.plot([l, u], [y0 + dy] * 2, color=c, lw=2.5, solid_capstyle="round")
            a.plot([h], [y0 + dy], marker=m, color=c, ms=8, mec=WHITE, mew=1.2)
        a.text(1.02, y0 + 0.14, fmt_hr(p["hr"]), fontsize=10.5, color=INK, va="center", transform=a.get_yaxis_transform())
        a.text(1.02, y0 - 0.14, fmt_hr(p["hr_s7"]), fontsize=10.5, color=MUTED, va="center", transform=a.get_yaxis_transform())
    a.axvline(1, color=INK, lw=1.2, ls=(0, (3, 3)))
    a.set_xscale("log")
    a.set_xlim(lo, hi)
    a.set_xticks([0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3], ["0.7", "0.8", "0.9", "1.0", "1.1", "1.2", "1.3"])
    a.minorticks_off()
    a.set_yticks([3, 2, 1, 0], [p["label"] for p in V["pathways"]], fontsize=12.5)
    for t, p in zip(a.get_yticklabels(), V["pathways"]):
        t.set_color(PW[p["label"]])
        t.set_fontweight("bold")
    a.set_ylim(-0.6, 3.6)
    a.set_xlabel("Hazard ratio per training-fold SD (log scale)", fontsize=11)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
    a.grid(axis="x", color=GRID)
    a.set_axisbelow(True)
    ax.plot([0.25, 0.55], [0.3, 0.3], color=PRIMARY, lw=2.5)
    ax.text(0.65, 0.3, "Primary (median representation)", fontsize=11, va="center")
    ax.plot([3.4, 3.7], [0.3, 0.3], color=S7, lw=2.5)
    ax.text(3.8, 0.3, "S-7 (representative probe, same arrays: not replication)", fontsize=11, va="center")
    f5 = V["five_set_e2f"]
    box(ax, 9.7, 1.0, 2.35, 2.6, None, fc="#e6f7f0", ec=FIVE, ls="--")
    ax.text(9.85, 1.3, "Sensitivity only", fontsize=12, fontweight="bold", color=FIVE)
    ax.text(9.85, 1.62, wrap(f"Five-set model E2F HR {fmt_hr(f5)}. Excludes 1 only after PI3K was added (collinear with p53); fixed in advance as sensitivity, not promoted.", 2.1, 10.5), fontsize=10.5, va="top")
    ax.text(9.7, 4.0, wrap("No primary-model interval excludes 1. Wald intervals are approximate under ridge penalisation.", 2.4, 11), fontsize=11, va="top", fontweight="bold")
    save(fig, "V14_forest")


# ------------------------------------------------------------------ V15 robustness
def v15(V):
    fig, ax = canvas()
    labels = [p["label"] for p in V["pathways"]]
    allc = {p["label"]: p["coef_all"] for p in V["pathways"]}
    for i, lab in enumerate(labels):
        a = fig.add_axes([0.065 + (i % 2) * 0.235, 0.6 - (i // 2) * 0.4, 0.19, 0.25])
        xs = [b["n"] for b in V["lobo_batches"]]
        ys = [b["coef"][lab] for b in V["lobo_batches"]]
        flips = [(y > 0) != (allc[lab] > 0) for y in ys]
        a.axhline(0, color=INK, lw=0.9, ls=(0, (3, 3)))
        a.axhline(allc[lab], color=PW[lab], lw=1.6)
        a.scatter([x for x, f in zip(xs, flips) if not f], [y for y, f in zip(ys, flips) if not f], s=40, color=PW[lab], zorder=3)
        a.scatter([x for x, f in zip(xs, flips) if f], [y for y, f in zip(ys, flips) if f], s=40, facecolor=WHITE, edgecolor=PW[lab], lw=1.6, zorder=3)
        a.set_title(f"{lab} · {sum(flips)} sign change{'s' if sum(flips) != 1 else ''}", fontsize=11, color=PW[lab], loc="left", fontweight="bold")
        a.set_xlim(0, 110)
        a.tick_params(labelsize=9)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        if i // 2 == 1:
            a.set_xlabel("patients in omitted batch", fontsize=9.5)
        if i % 2 == 0:
            a.set_ylabel("log-HR", fontsize=9.5)
    ax.text(0.3, 0.3, "Leave one processing batch out (9 refits)", fontsize=12.5, fontweight="bold")
    m = fig.add_axes([0.66, 0.6, 0.3, 0.27])
    for i, p in enumerate(V["pathways"]):
        m.barh(3 - i, p["mmr_e2"], color=PW[p["label"]], height=0.6)
        m.text(p["mmr_e2"] + 0.004, 3 - i, f"{p['mmr_e2']:.4f}", va="center", fontsize=10)
    m.set_yticks([3, 2, 1, 0], [p["label"] for p in V["pathways"]], fontsize=10.5)
    m.set_xlim(0, 0.22)
    m.set_xlabel("ε², association with MMR status (n = 410, exploratory)", fontsize=9.5)
    m.tick_params(labelsize=9)
    for s in ("top", "right"):
        m.spines[s].set_visible(False)
    ax.text(7.5, 0.3, "Scores track molecular subtype", fontsize=12.5, fontweight="bold")
    c = fig.add_axes([0.66, 0.22, 0.3, 0.2])
    sl = [("M1 clinical", V["slopes"]["M1"], INK), ("M2 four-set", V["slopes"]["M2_primary"], PRIMARY), ("Five-set (sens.)", V["slopes"]["M2_sensitivity"], FIVE)]
    for i, (l, v, col) in enumerate(sl):
        c.barh(2 - i, v, color=col, height=0.55)
        c.text(v + 0.01, 2 - i, f"{v:.3f}", va="center", fontsize=10)
    c.axvline(1, color=INK, lw=1, ls=(0, (3, 3)))
    c.set_yticks([2, 1, 0], [s[0] for s in sl], fontsize=10.5)
    c.set_xlim(0.6, 1.05)
    c.set_xlabel("corrected calibration slope (1 = ideal)", fontsize=9.5)
    c.tick_params(labelsize=9)
    for s in ("top", "right"):
        c.spines[s].set_visible(False)
    ax.text(7.5, 3.15, "Residual overfitting", fontsize=12.5, fontweight="bold")
    ax.text(0.3, 5.45, wrap(f"Omitted batches hold {V['lobo_sizes'][0]}–{V['lobo_sizes'][1]} patients, so the spread mixes batch identity with sample size. "
                       "Open circle = sign differs from the all-patient estimate (solid line). Batch effects were characterised, not resolved.", 11.6, 10.5), fontsize=10.5, color=INK, va="top")
    save(fig, "V15_robustness")


# ------------------------------------------------------------------ V16 what this project established
def _axes(fig, x, y, w, h):
    """Inset axes from canvas inches (top-left origin)."""
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])


def _block(ax, x, y, w, h, num, title, col, bg, status):
    box(ax, x, y, w, h, None, fc=bg, ec=col, lw=1.5)
    ax.add_patch(Rectangle((x, y + 0.04), 0.07, h - 0.08, fc=col, ec="none"))
    ax.add_patch(Circle((x + 0.36, y + 0.3), 0.17, fc=col, ec="none"))
    ax.text(x + 0.36, y + 0.31, str(num), ha="center", va="center", fontsize=12, color=WHITE, fontweight="bold")
    ax.text(x + 0.62, y + 0.3, title, fontsize=13, fontweight="bold", color=col, va="center")
    ax.text(x + w - 0.14, y + 0.3, status, fontsize=8.4, fontweight="bold", color=col, va="center", ha="right",
            bbox=dict(boxstyle="round,pad=0.25,rounding_size=0.1", fc=WHITE, ec=col, lw=0.9))


def _lines(ax, x, y, w, items, fs=10.6, gap=0.08, color=INK):
    for t, c, b in items:
        tw = wrap(t, w, fs)
        ax.text(x, y, tw, fontsize=fs, va="top", linespacing=1.2, color=c or color, fontweight="bold" if b else "normal")
        y += (tw.count("\n") + 1) * fs * 1.2 / 72 + gap
    return y


def _mini(a):
    for sp in ("top", "right", "left"):
        a.spines[sp].set_visible(False)
    a.tick_params(axis="y", length=0, labelsize=9.6)
    a.tick_params(axis="x", labelsize=9)
    a.patch.set_alpha(0)


def v16(V):
    E, X = V["E"], V["X"]
    fig, ax = canvas()
    TOOL, BIO, P53, PRED = ("#6650c4", "#f1eefb"), ("#1f7a4d", "#e9f5ee"), ("#a86b12", "#fbf4e7"), (PRIMARY, "#eaf2fc")
    # 1 · tool
    x0, y0, w0, h0 = 0.08, 0.06, 5.98, 2.92
    _block(ax, x0, y0, w0, h0, 1, "A reusable pathway-audit framework", *TOOL, "TESTED SOFTWARE")
    _lines(ax, x0 + 0.24, y0 + 0.6, w0 - 0.4, [("Python package and command-line tool that checks what a bulk gene-set score measures before it is interpreted.", None, False)])
    steps = ["Data-quality\ngate", "Gene drivers,\ncoherence", "Tissue\ncomposition", "Molecular\nassociations", "Random-set\nnull", "Survival\nscreen"]
    cw, cg = 0.88, 0.075
    for i, s in enumerate(steps):
        cx = x0 + 0.24 + i * (cw + cg)
        box(ax, cx, y0 + 1.2, cw, 0.6, None, fc=WHITE, ec=TOOL[0], lw=1.1)
        ax.text(cx + cw / 2, y0 + 1.5, s, ha="center", va="center", fontsize=8.8, linespacing=1.15)
        if i < len(steps) - 1:
            act(ax, cx + cw, y0 + 1.52, cx + cw + cg, y0 + 1.52, color=TOOL[0], lw=1.2)
    facts = [("11/11", "package tests pass"), ("4/4", "planted truths recovered in the demo"), ("5", "cohorts in the case study (1 failed the QC gate)")]
    for i, (n_, t) in enumerate(facts):
        fx = x0 + 0.24 + i * 1.92
        ax.text(fx, y0 + 2.04, n_, fontsize=17, fontweight="bold", color=TOOL[0], va="center")
        ax.text(fx, y0 + 2.24, wrap(t, 1.78, 9.6), fontsize=9.6, color=INK, va="top", linespacing=1.15)
    ax.text(x0 + 0.24, y0 + h0 - 0.13, "Docker image and GitHub Actions workflow written, not yet run.", fontsize=9.2, color=MUTED, va="bottom", style="italic")
    # 2 · EMT
    x1, y1, w1, h1 = 6.14, 0.06, 5.98, 2.92
    _block(ax, x1, y1, w1, h1, 2, "Replicated: bulk EMT score tracks stroma", *BIO, "REPLICATED · 4 COHORTS")
    yy = _lines(ax, x1 + 0.24, y1 + 0.6, 2.3, [("Gene-disjoint ρ 0.94–0.98 with fibroblast markers in all four usable cohorts.", None, True),
                                              ("Each above its own random-gene-set baseline (grey band, 95 % range).", None, False)], fs=10.3)
    # bulk-tissue icon: tumour cells + fibroblasts pooled into one score
    ix, iy = x1 + 0.3, y1 + 1.86
    ax.add_patch(Circle((ix + 0.42, iy + 0.32), 0.36, fc=CELL["matrix"], ec=FAINT, lw=1))
    for dx, dy in ((0.22, 0.2), (0.4, 0.14), (0.3, 0.4), (0.52, 0.34)):
        ax.add_patch(Circle((ix + dx, iy + dy), 0.065, fc=CELL["epi"], ec=CELL["epin"], lw=0.6))
    for dx, dy, a_ in ((0.62, 0.22, 30), (0.58, 0.5, -20), (0.2, 0.52, 10)):
        ax.add_patch(Ellipse((ix + dx, iy + dy), 0.24, 0.07, angle=a_, fc=CELL["fib"], ec=CELL["fibn"], lw=0.6))
    act(ax, ix + 0.82, iy + 0.32, ix + 1.12, iy + 0.32, color=MUTED, lw=1.2)
    ax.text(ix + 1.16, iy + 0.32, "one bulk\nscore", fontsize=9, color=MUTED, va="center")
    ax.text(x1 + 0.24, y1 + h0 - 0.13, "Bulk-tissue association; not tumour-cell EMT.", fontsize=9.6, color=BIO[0], va="bottom", fontweight="bold")
    a = _axes(fig, x1 + 3.42, y1 + 0.62, 2.4, 1.85)
    order = ["GSE39582", "GSE41258", "GSE143985", "TCGA-COAD"]
    for i, c in enumerate(order):
        yv = len(order) - 1 - i
        b = E["baseline"][c]
        a.plot([b["p2_5"], b["p97_5"]], [yv, yv], color=FAINT, lw=6, solid_capstyle="butt", alpha=0.75)
        est, lo, hi = E["disjoint"][c]
        a.plot([lo, hi], [yv, yv], color=INK, lw=1.4)
        a.plot(est, yv, "o", ms=7, color=BIO[0], mec=WHITE, mew=1)
        a.text(1.04, yv, f"{est:.2f}", va="center", fontsize=9.4, color=INK)
    a.axvline(0, color=MUTED, lw=0.8)
    a.set_yticks(range(4), order[::-1])
    a.set_xlim(-0.8, 1.3)
    a.set_xticks([-0.5, 0, 0.5, 1.0])
    a.set_ylim(-0.6, 3.6)
    a.set_xlabel("Spearman ρ", fontsize=9, labelpad=1)
    _mini(a)
    # 3 · p53 anchoring
    x2, y2, w2, h2 = 0.08, 3.06, 7.22, 1.98
    _block(ax, x2, y2, w2, h2, 3, "Molecular anchoring: p53 score tracks TP53", *P53, "REPLICATED · 2 OF 3 EXTERNAL")
    _lines(ax, x2 + 0.24, y2 + 0.58, 3.45, [("Lower in TP53-mutant tumours; replicated in GSE41258 and TCGA-COAD, and after excluding hypermutated tumours.", None, True),
                                          ("An mRNA-score association, not a measure of p53 protein activity or pathway function.", None, False),
                                          ("Not established: E2F–TP53 and p53–MSI (not replicated); TCGA BRAF links (explained by hypermutation).", MUTED, False)], fs=9.6, gap=0.05)
    a = _axes(fig, x2 + 5.0, y2 + 0.58, 1.6, 1.0)
    r1 = {r["cohort"]: r for r in values.csv.DictReader(open(values.ROOT / "results/research_extensions/e3_replication_tests.csv")) if r["hypothesis"] == "R1"}
    hx = E["tcga_hyper"]["R1"]["excluding_hypermutated"]["rank_biserial"]
    rows = [("GSE41258", (float(r1["GSE41258"]["effect"]), float(r1["GSE41258"]["lo"]), float(r1["GSE41258"]["hi"])), True),
            ("GSE143985", (float(r1["GSE143985"]["effect"]), float(r1["GSE143985"]["lo"]), float(r1["GSE143985"]["hi"])), False),
            ("TCGA-COAD", (E["r1_tcga"], *E["r1_tcga_ci"]), True),
            ("TCGA, non-hypermut.", tuple(hx), True)]
    for i, (lab, (e_, lo, hi), ok) in enumerate(rows):
        yv = len(rows) - 1 - i
        a.plot([lo, hi], [yv, yv], color=INK if ok else FAINT, lw=1.4)
        a.plot(e_, yv, "o", ms=6.5, color=P53[0] if ok else WHITE, mec=P53[0] if ok else FAINT, mew=1.2)
        a.text(0.25, yv, f"{e_:+.2f}" + ("" if ok else " n.s."), va="center", fontsize=9, color=INK if ok else MUTED)
    a.axvline(0, color=MUTED, lw=0.8)
    a.set_yticks(range(4), [r[0] for r in rows][::-1])
    a.set_xlim(-0.65, 0.62)
    a.set_xticks([-0.5, 0])
    a.set_ylim(-0.6, 3.6)
    a.set_xlabel("rank-biserial, mutant vs wild-type", fontsize=8.6, labelpad=1)
    _mini(a)
    # 4 · prediction
    x3, y3, w3, h3 = 7.38, 3.06, 4.74, 1.98
    _block(ax, x3, y3, w3, h3, 4, "Rigorous predictive evaluation", *PRED, "PRE-SPECIFIED")
    _lines(ax, x3 + 0.24, y3 + 0.58, 2.35, [(f"Corrected ΔC +{V['theta']:.6f}; all five 95 % intervals include zero.", None, True),
                                          ("No reproducible gain in recurrence discrimination. Does not rule out a biological role.", None, False)], fs=9.8, gap=0.06)
    a = _axes(fig, x3 + 2.75, y3 + 0.6, 1.8, 0.95)
    for i, (lo, hi) in enumerate(V["runs"]):
        a.plot([lo, hi], [4 - i, 4 - i], color=PRIMARY, lw=1.6)
        a.plot(V["theta"], 4 - i, "|", ms=7, color=INK)
    a.axvline(0, color=S7, lw=1.1, ls="--")
    a.set_xlim(-0.02, 0.06)
    a.set_xticks([0, 0.05])
    a.set_yticks([])
    a.set_ylim(-0.5, 4.5)
    a.set_xlabel("corrected ΔC, 5 runs", fontsize=8.6, labelpad=1)
    _mini(a)
    # what remains unresolved
    yb = 5.12
    ax.text(0.1, yb + 0.02, "What remains unresolved", fontsize=11.2, fontweight="bold", color=INK, va="top")
    gaps = [("Tumour-cell-specific EMT", "needs spatial, single-cell or cell-resolved protein evidence"),
            ("Direct pathway activation", "needs functional or protein-level measurement"),
            ("External predictive validation", f"no eligible cohort: {X['ext_eligible']} of {X['ext_assessed']} assessed; max {X['ext_max_events']} vs ≈{X['ext_needed']} recurrences needed"),
            ("Recurrence-time unit", "undocumented; provisionally months")]
    gw = 3.0
    for i, (h_, t) in enumerate(gaps):
        gx = 0.1 + i * (gw + 0.02)
        ax.plot([gx, gx], [yb + 0.28, yb + 0.86], color=FAINT, lw=2)
        ax.text(gx + 0.1, yb + 0.3, h_, fontsize=10, fontweight="bold", color=INK, va="top")
        ax.text(gx + 0.1, yb + 0.52, wrap(t, gw - 0.25, 9.4), fontsize=9.4, color=MUTED, va="top", linespacing=1.15)
    save(fig, "V16_project_established")


# ------------------------------------------------------------------ V17 decisions timeline
def v17(V):
    fig, ax = canvas()
    stages = ["Foundation", "Cohort and endpoint", "Gene mapping", "Feasibility", "Protocol and panel", "QC and exploration", "Model construction", "Nested bootstrap", "Verification"]
    dec = {1: "DC3 CIT partition internal only", 4: "DC1 fix every choice before outcomes\nDC2 clinical block unpenalised",
           7: "DC4 seed deviation adjudicated, full re-run\nDC5 pre-specified 500 → 750 top-up", 8: "DC6 misleading results excluded\nDC7 audits commissioned and reviewed"}
    x0, x1, y = 0.6, 11.6, 1.2
    ax.plot([x0, x1], [y, y], color=GRID, lw=3)
    for i, s in enumerate(stages):
        x = x0 + i * (x1 - x0) / 8
        hasd = i in dec
        ax.add_patch(Circle((x, y), 0.2, fc="#6650c4" if hasd else WHITE, ec="#6650c4" if hasd else FAINT, lw=2, zorder=3))
        ax.text(x, y + 0.04, str(i + 1), ha="center", va="center", fontsize=9.5, color=WHITE if hasd else MUTED, fontweight="bold", zorder=4)
        ax.text(x, y - 0.42, wrap(s, 1.25, 10), ha="center", va="bottom", fontsize=10, color=INK)
    bx = {1: 0.5, 4: 3.75, 7: 7.0, 8: 9.6}
    for i, t in dec.items():
        x = x0 + i * (x1 - x0) / 8
        ax.plot([x, bx[i] + 1.225], [y + 0.2, 1.8], color="#6650c4", lw=1.2)
        box(ax, bx[i], 1.8, 2.45, 1.05, t, fc="#f1eefb", ec="#6650c4", fs=10)
    facts = [("60 / 60", "frozen verifier checks passed before reporting"), ("2.2×10⁻¹⁶", "maximum difference in the independent re-computation (from checkpoints; not a model refit)"),
             ("138", "reported quantities recomputed independently"), ("1,097", "automated tests passing")]
    for i, (n_, t) in enumerate(facts):
        x = 0.2 + i * 3.0
        ax.text(x, 3.5, n_, fontsize=22, fontweight="bold", color=INK)
        ax.text(x, 3.9, wrap(t, 2.75, 11), fontsize=11, color=MUTED, va="top")
    ax.text(0.2, 5.3, "My contribution: scientific direction and methodological decisions, each recorded in the project's decision records (DC1–DC7).", fontsize=11.5, color=INK)
    save(fig, "V17_decisions_audit")


# ------------------------------------------------------------------ V18 overall conclusion + next questions
def v18(V):
    X = V["X"]
    fig, ax = canvas()
    ax.text(0.15, 0.3, "How the project answered its question", fontsize=13.2, fontweight="bold", color=INK, va="center")
    steps = [("Asked", "Do four pathway scores improve recurrence prediction after surgery?", PRIMARY, "#eaf2fc"),
             ("Found", f"No reproducible gain: corrected ΔC +{V['theta']:.6f}, every interval includes zero.", PRIMARY, "#eaf2fc"),
             ("Asked next", "What do these scores measure in bulk tumour tissue?", "#1f7a4d", "#e9f5ee"),
             ("Found", "Composition- and genotype-linked expression: EMT score tracks stroma; p53 score tracks TP53 status.", "#1f7a4d", "#e9f5ee"),
             ("Delivered", "pathway-audit, so these checks can be run before a score is interpreted.", "#6650c4", "#f1eefb")]
    y = 0.62
    for i, (lab, t, col, bg) in enumerate(steps):
        box(ax, 0.15, y, 5.25, 0.62, None, fc=bg, ec=col, lw=1.2)
        ax.text(0.3, y + 0.31, lab, fontsize=10.4, fontweight="bold", color=col, va="center")
        ax.text(1.42, y + 0.31, wrap(t, 3.85, 10.6), fontsize=10.6, va="center", linespacing=1.15)
        if i < len(steps) - 1:
            act(ax, 0.75, y + 0.62, 0.75, y + 0.74, color=MUTED, lw=1.2)
        y += 0.74
    box(ax, 0.15, 4.42, 5.25, 1.22, None, fc=INK, ec=INK)
    ax.text(0.35, 4.66, "Take-home", fontsize=10.6, fontweight="bold", color="#9fd3ff", va="center")
    ax.text(0.35, 4.9, wrap("Audit bulk pathway scores for tissue composition, genotype confounding and data quality before reading them as biology or adding them to prognostic models.", 4.85, 11.4),
            fontsize=11.4, color=WHITE, va="top", linespacing=1.2)
    ax.text(5.85, 0.3, "Next scientific questions", fontsize=13.2, fontweight="bold", color=INK, va="center")
    tag(ax, 8.65, 0.3, "nt", fs=8.4)
    qs = [("Is the EMT-associated expression in tumour cells or in stroma?",
           "Spatial transcriptomics, single-cell RNA-seq or cell-resolved protein staining (e.g. E-cadherin, vimentin) on matched tumours.",
           "Separates tumour-cell EMT from fibroblast content, which bulk data cannot do."),
          ("Are the Wnt, E2F, p53 and EMT pathways actually active?",
           "Functional or protein-level readouts: reporter assays, phospho-protein and target-protein measurement.",
           "Turns expression associations into evidence about pathway function."),
          ("Does any score add recurrence information in a new cohort?",
           f"A pre-specified test in an independent stage II–III cohort with ≈{X['ext_needed']} recurrences and documented follow-up time.",
           f"None of the {X['ext_assessed']} public cohorts assessed was eligible (largest: {X['ext_max_events']} recurrences).")]
    y = 0.62
    for i, (q, how, why) in enumerate(qs):
        box(ax, 5.85, y, 6.2, 1.58, None, fc=WHITE, ec=GRID, lw=1.2)
        ax.add_patch(Circle((6.17, y + 0.32), 0.17, fc=PANEL, ec=MUTED, lw=1))
        ax.text(6.17, y + 0.33, f"Q{i + 1}", ha="center", va="center", fontsize=9, fontweight="bold", color=MUTED)
        ax.text(6.46, y + 0.32, wrap(q, 5.45, 11.6), fontsize=11.6, fontweight="bold", va="center")
        ax.text(6.05, y + 0.62, "How:", fontsize=10.2, fontweight="bold", color=MUTED, va="top")
        ax.text(6.62, y + 0.62, wrap(how, 5.3, 10.2), fontsize=10.2, va="top", linespacing=1.15)
        ax.text(6.05, y + 1.14, "Why:", fontsize=10.2, fontweight="bold", color=MUTED, va="top")
        ax.text(6.62, y + 1.14, wrap(why, 5.3, 10.2), fontsize=10.2, va="top", linespacing=1.15, color=MUTED)
        y += 1.73
    save(fig, "V18_conclusion_questions")


def v19(V):
    E = V["E"]
    fig, ax = canvas()
    ax.text(0.15, 0.35, "1 · Replicated finding: the EMT score reads stroma", fontsize=13, fontweight="bold")
    a = fig.add_axes([0.13, 0.43, 0.33, 0.42])
    order = ["GSE39582", "GSE41258", "GSE143985", "TCGA-COAD"]
    for i, c in enumerate(order):
        est, lo, hi = E["disjoint"][c]
        y = len(order) - i
        a.barh(y, est, color=PW["EMT"], height=0.55)
        a.plot([lo, hi], [y, y], color=INK, lw=1.6)
        a.text(hi + 0.02, y, f"{est:.2f}", va="center", fontsize=11)
    a.barh(0, 0, height=0.55)
    a.text(0.18, 0, f"failed quality gate (control r {E['qc']['GSE103479']:.2f})", va="center", fontsize=10.5, color=MUTED)
    for i, c in enumerate(order):
        b = E["baseline"][c]
        y = len(order) - i
        a.plot([max(b["p2_5"], 0), max(b["p97_5"], 0.005)], [y + 0.36, y + 0.36], color=MUTED, lw=3, alpha=0.6, solid_capstyle="butt")
    a.text(0.0, 4.62, "grey: random gene sets, 95 % range (cohort-specific; ≤ 0 not shown)", fontsize=9, color=MUTED)
    a.set_yticks([4, 3, 2, 1, 0], order + ["GSE103479"], fontsize=11)
    a.set_xlim(0, 1.12)
    a.set_ylim(-0.6, 5.0)
    a.set_xlabel(f"Spearman ρ, EMT score vs fibroblast markers ({E['n_shared_fib']} shared genes removed)", fontsize=10)
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)
    ax.text(6.35, 0.35, "2 · Delivered as a tool: pathway-audit", fontsize=13, fontweight="bold")
    steps = [("Data-quality gate", "control gene pairs must co-express"), ("Gene drivers", "coherence, leave-one-out, stable modules"),
             ("Composition", "stromal / immune estimates, gene-disjoint test"), ("Genotype effects", "robust to, or explained by, MMR and composition"),
             ("Survival screen", "optimism-corrected gain over clinical model")]
    for i, (h, t) in enumerate(steps):
        y = 0.65 + i * 0.6
        box(ax, 6.35, y, 5.7, 0.5, None, fc=PANEL, ec=GRID)
        ax.text(6.5, y + 0.25, f"{i + 1}  {h}", fontsize=11.5, fontweight="bold", va="center")
        ax.text(8.75, y + 0.25, t, fontsize=10.5, color=MUTED, va="center")
    ax.text(6.35, 3.85, "Case study, five colorectal cohorts", fontsize=11.5, fontweight="bold")
    lines = [f"EMT composition-dominated in all 4 usable cohorts", f"GSE103479 stopped at the quality gate",
             f"BRAF → p53 score explained by MMR and composition ({E['braf_p53'][0]:+.2f} → {E['braf_p53'][1]:+.2f} SD)",
             f"TCGA: BRAF effects explained by hypermutation; TP53 → p53 score robust",
             f"Survival screen on GSE39582: corrected ΔC +{E['survival_screen']:.4f} (frozen +{V['theta']:.4f})"]
    for i, t in enumerate(lines):
        ax.text(6.45, 4.17 + i * 0.3, "•  " + t, fontsize=10.4, va="center")
    ax.text(0.15, 5.65, "Post-hoc registered extensions, outcome-free; tool tested (11/11) incl. ESTIMATE match to 1e-14. Composition values are estimates, not measured fractions.",
            fontsize=10, color=MUTED)
    save(fig, "V19_solution_tool")


ALL = [v01, v02, v03, v04, v05, v06, v07, v08, v10, v11, v12, v13, v14, v15, v16, v17, v18, v19]

if __name__ == "__main__":
    V = values.load()
    V["X"] = values.load_exploratory()
    V["E"] = values.load_extensions()
    for f in ALL:
        f(V)
    print("wrote", len(list(OUT.glob("*.svg"))), "SVG and", len(list(OUT.glob("*.png"))), "PNG")
