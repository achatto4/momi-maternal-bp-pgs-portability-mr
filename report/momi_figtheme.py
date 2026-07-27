#!/usr/bin/env python3
"""
momi_figtheme.py -- shared publication theme for all MOMI figures.

One place for palette, typography, cohort naming and the save helper, so every figure
in the manuscript and deck looks like it came from the same study. Colourblind-safe
(Okabe-Ito) throughout. Imported by build_figures.py.
"""
import os
import matplotlib as mpl
import matplotlib.pyplot as plt

# ---- Okabe-Ito colourblind-safe palette -------------------------------------
OI = dict(black="#000000", orange="#E69F00", skyblue="#56B4E9", green="#009E73",
          yellow="#F0E442", blue="#0072B2", vermillion="#D55E00", purple="#CC79A7",
          grey="#999999")

# Evidence lines in the triangulation (distinct hue + marker, not colour alone)
ARM = {
    "MR_external":   dict(color=OI["green"],      marker="o", label="MR, external (European)"),
    "observational": dict(color=OI["blue"],       marker="s", label="Observational (adjusted)"),
    "MR_ours_all":   dict(color=OI["black"],       marker="D", label="MR, this study (all cohorts)"),
    "MR_ours_SAS":   dict(color=OI["orange"],      marker="^", label="MR, this study (South Asian)"),
    "MR_ours_AFR":   dict(color=OI["vermillion"], marker="v", label="MR, this study (African)"),
}

# Ancestry accent (score training ancestry)
ANC = {"EUR": OI["blue"], "SAS": OI["orange"], "AFR": OI["vermillion"],
       "EAS": OI["green"], "DIVERSE": OI["purple"], "MVP": OI["purple"]}

# Sequential map for the transferability heatmap
HEATMAP_CMAP = "cividis"          # perceptually uniform + colourblind-safe

# ---- cohort naming and ordering ---------------------------------------------
# canonical short two-line labels, grouped South Asian then African
COHORT = {
    "AMANHI-Bangladesh":      dict(short="Sylhet",  sub="AMANHI · Bangladesh", anc="SAS"),
    "AMANHI (Sylhet, Bangladesh)": dict(short="Sylhet", sub="AMANHI · Bangladesh", anc="SAS"),
    "AMANHI-Karachi":         dict(short="Karachi", sub="AMANHI · Pakistan",   anc="SAS"),
    "AMANHI (Karachi, Pakistan)": dict(short="Karachi", sub="AMANHI · Pakistan", anc="SAS"),
    "GAPPS-Bangladesh":       dict(short="Matlab",  sub="PreSSMat · Bangladesh", anc="SAS"),
    "PreSSMat (Matlab, Bangladesh)": dict(short="Matlab", sub="PreSSMat · Bangladesh", anc="SAS"),
    "AMANHI-Pemba":           dict(short="Pemba",   sub="AMANHI · Tanzania",   anc="AFR"),
    "AMANHI (Pemba, Tanzania)": dict(short="Pemba", sub="AMANHI · Tanzania",   anc="AFR"),
    "ZAPPS-Lusaka":           dict(short="Lusaka",  sub="ZAPPS · Zambia",      anc="AFR"),
    "ZAPPS (Lusaka, Zambia)": dict(short="Lusaka",  sub="ZAPPS · Zambia",      anc="AFR"),
}
COHORT_ORDER = ["Sylhet", "Karachi", "Matlab", "Pemba", "Lusaka"]  # SAS then AFR

def cohort_short(name):
    return COHORT.get(str(name).strip(), dict(short=str(name))).get("short", str(name))

def cohort_label(name):
    d = COHORT.get(str(name).strip())
    return f"{d['short']}\n{d['sub']}" if d else str(name)

# ---- global rcParams --------------------------------------------------------
def use_theme():
    mpl.rcParams.update({
        "figure.dpi": 120,
        "savefig.dpi": 300,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.labelsize": 10,
        "axes.edgecolor": "#444444",
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": "#DDDDDD",
        "grid.linewidth": 0.6,
        "xtick.color": "#333333",
        "ytick.color": "#333333",
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.frameon": False,
        "legend.fontsize": 9,
        "figure.titlesize": 12,
        "figure.titleweight": "bold",
    })

def save(fig, name, outdir, pad=0.25):
    """Write both PDF (vector, for the manuscript) and PNG (for the deck/preview)."""
    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"{name}.{ext}"), bbox_inches="tight", pad_inches=pad)
    plt.close(fig)
