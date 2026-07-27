#!/usr/bin/env python3
"""
build_figures.py -- regenerate the manuscript's publication figures from the results
tables, using one shared theme (momi_figtheme). Runs locally from the .tsv tables; no
cluster or R needed. Run AFTER the analysis pipeline and BEFORE build_paper.py:

  python3 report/build_figures.py --results DIR --out DIR/figures

Writes PDF (manuscript) + PNG (deck/preview) for each figure.
"""
import argparse, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import momi_figtheme as T

ap = argparse.ArgumentParser()
ap.add_argument("--results", required=True)
ap.add_argument("--out", default=None)
A = ap.parse_args()
TBL = os.path.join(A.results, "tables")
OUT = A.out or os.path.join(A.results, "figures")
T.use_theme()

def tsv(name):
    return pd.read_csv(os.path.join(TBL, name + ".tsv"), sep="\t", dtype=str)

def num(s):
    return pd.to_numeric(s, errors="coerce")

# canonical cohort ordering key
def coh_key(c):
    s = T.cohort_short(c)
    return T.COHORT_ORDER.index(s) if s in T.COHORT_ORDER else 99

CANON = {"AMANHI-Pakistan": "AMANHI-Karachi", "GAPPS-Zambia": "ZAPPS-Lusaka"}

# =====================================================================
# F2 -- transferability heatmap (score ancestry x cohort, per trait)
# =====================================================================
def fig_F2():
    d = tsv("S4_transfer_grid")
    d = d[d["definition"] == "mean"].copy()
    d["cohort"] = d["cohort"].replace(CANON)
    d["short"] = d["cohort"].map(T.cohort_short)
    d["R2"] = num(d["R2pct"])
    anc_order = ["EUR", "SAS", "DIVERSE", "MVP", "EAS"]
    anc_lab = {"EUR": "European", "SAS": "South Asian", "DIVERSE": "Multi-ancestry",
               "MVP": "MVP (multi)", "EAS": "East Asian"}
    cohorts = [c for c in T.COHORT_ORDER if c in set(d["short"])]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    vmax = np.nanpercentile(d["R2"], 98)
    for ax, trait in zip(axes, ["SBP", "DBP"]):
        sub = d[d["trait"] == trait]
        M = np.full((len(anc_order), len(cohorts)), np.nan)
        for i, a in enumerate(anc_order):
            for j, c in enumerate(cohorts):
                v = sub[(sub["score_anc"] == a) & (sub["short"] == c)]["R2"]
                if len(v): M[i, j] = np.nanmax(v.values)
        im = ax.imshow(M, cmap=T.HEATMAP_CMAP, vmin=0, vmax=vmax, aspect="auto")
        ax.set_xticks(range(len(cohorts)))
        ax.set_xticklabels([T.cohort_label(c if c not in T.COHORT_ORDER else
                            [k for k in T.COHORT.values() if k["short"] == c][0] and c)
                            for c in cohorts], fontsize=8.5)
        # use short two-line labels
        ax.set_xticklabels([f"{c}\n{'South Asian' if c in ['Sylhet','Karachi','Matlab'] else 'African'}"
                            for c in cohorts], fontsize=8.5)
        ax.set_yticks(range(len(anc_order)))
        ax.set_yticklabels([anc_lab[a] for a in anc_order] if ax is axes[0] else [])
        ax.set_title(f"{trait}", fontsize=11)
        ax.grid(False)
        for i in range(len(anc_order)):
            for j in range(len(cohorts)):
                if not np.isnan(M[i, j]):
                    frac = M[i, j] / vmax if vmax else 0   # cividis: low=dark navy, high=yellow
                    ax.text(j, i, f"{M[i,j]:.1f}", ha="center", va="center",
                            fontsize=8.5, color="white" if frac < 0.55 else "#222222")
        # divider between SAS and AFR cohorts
        nsas = sum(1 for c in cohorts if c in ["Sylhet", "Karachi", "Matlab"])
        ax.axvline(nsas - 0.5, color="white", lw=2.5)
        ax.set_xticks(np.arange(-.5, len(cohorts), 1), minor=True)
        ax.set_yticks(np.arange(-.5, len(anc_order), 1), minor=True)
        ax.grid(which="minor", color="white", lw=1.2)
        ax.tick_params(which="minor", length=0)
    axes[0].set_ylabel("Polygenic score training ancestry")
    cb = fig.colorbar(im, ax=axes, shrink=0.8, pad=0.02)
    cb.set_label("Incremental $R^2$ (%)")
    T.save(fig, "F2_transfer", OUT)
    print("  F2_transfer")

# =====================================================================
# F4 / F4b -- triangulation small multiples
# =====================================================================
def _arm_key(row):
    if row["arm"] == "MR_ours":
        return f"MR_ours_{row['stratum']}"
    return row["arm"]

def _triangulation(outcomes, value_cols, xlabel, logx, null, fname, title):
    d = tsv("T4_triangulation")
    d["armkey"] = d.apply(_arm_key, axis=1)
    est, lo, hi = value_cols
    d["est"], d["lo"], d["hi"] = num(d[est]), num(d[lo]), num(d[hi])
    d["Fmin"], d["Fmax"] = num(d["F_min"]), num(d["F_max"])
    row_order = ["MR_external", "observational", "MR_ours_all", "MR_ours_SAS", "MR_ours_AFR"]
    traits = ["SBP", "DBP"]
    nR, nC = len(outcomes), len(traits)
    fig, axes = plt.subplots(nR, nC, figsize=(10.5, 1.35 * len(row_order) * nR / 2 + 1.2),
                             sharex="col", squeeze=False)
    for r, oc in enumerate(outcomes):
        # x-limits driven by the ESTIMATES (not the huge underpowered CIs), so the
        # informative points stay readable; CIs beyond the frame get an arrow cap.
        est_row = num(d[d["outcome"] == oc]["est"]).dropna()
        if logx:
            lo_lim = min(est_row.min(), null) * 0.55
            hi_lim = max(est_row.max(), null) * 1.8
        else:
            span = max(abs(est_row.min() - null), abs(est_row.max() - null))
            lo_lim, hi_lim = null - span * 1.5, null + span * 1.5
        for c, tr in enumerate(traits):
            ax = axes[r][c]
            sub = d[(d["outcome"] == oc) & (d["trait"] == tr)]
            for y, ak in enumerate(row_order):
                rr = sub[sub["armkey"] == ak]
                yy = len(row_order) - 1 - y
                if not len(rr):
                    continue
                rr = rr.iloc[0]
                if np.isnan(rr["est"]):
                    continue
                st = T.ARM[ak]
                hollow = ak == "MR_ours_AFR"
                # clip whisker to frame, arrow-cap where it exceeds
                clo, chi = max(rr["lo"], lo_lim), min(rr["hi"], hi_lim)
                ax.plot([clo, chi], [yy, yy], color=st["color"], lw=1.6, zorder=2,
                        solid_capstyle="butt")
                if rr["lo"] < lo_lim:
                    ax.annotate("", xy=(lo_lim, yy), xytext=(lo_lim * 1.06 if logx else lo_lim + (hi_lim-lo_lim)*0.02, yy),
                                arrowprops=dict(arrowstyle="-|>", color=st["color"], lw=1.4))
                if rr["hi"] > hi_lim:
                    ax.annotate("", xy=(hi_lim, yy), xytext=(hi_lim * 0.94 if logx else hi_lim - (hi_lim-lo_lim)*0.02, yy),
                                arrowprops=dict(arrowstyle="-|>", color=st["color"], lw=1.4))
                if lo_lim <= rr["est"] <= hi_lim:
                    ax.scatter([rr["est"]], [yy], s=46, marker=st["marker"],
                               facecolor="white" if hollow else st["color"],
                               edgecolor=st["color"], linewidth=1.6, zorder=3)
                if ak in ("MR_ours_SAS", "MR_ours_AFR") and not np.isnan(rr["Fmin"]):
                    ax.text(0.985, yy, f"F {rr['Fmin']:.0f}–{rr['Fmax']:.0f}",
                            transform=ax.get_yaxis_transform(), ha="right", va="center",
                            fontsize=6.5, color="#888888",
                            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85))
            ax.axvline(null, color="#555555", ls="--", lw=1)
            if logx:
                ax.set_xscale("log")
            ax.set_xlim(lo_lim, hi_lim)
            ax.set_yticks(range(len(row_order)))
            ax.set_yticklabels([T.ARM[a]["label"] for a in reversed(row_order)]
                               if c == 0 else [])
            ax.set_ylim(-0.6, len(row_order) - 0.4)
            ax.grid(axis="y", visible=False)
            ax.grid(axis="x", color="#EEEEEE")
            if r == 0:
                ax.set_title(tr)
            if c == nC - 1:
                ax.text(1.02, 0.5, oc, transform=ax.transAxes, rotation=-90,
                        va="center", ha="left", fontweight="bold", fontsize=10)
        # pad xlim
        for c in range(nC):
            pass
    for c in range(nC):
        axes[-1][c].set_xlabel(xlabel)
    handles = [Line2D([0], [0], color=T.ARM[a]["color"], marker=T.ARM[a]["marker"],
               lw=1.6, markersize=7,
               markerfacecolor="white" if a == "MR_ours_AFR" else T.ARM[a]["color"],
               markeredgecolor=T.ARM[a]["color"], label=T.ARM[a]["label"])
               for a in row_order]
    fig.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout(rect=[0, 0.02, 1, 0.97])
    T.save(fig, fname, OUT)
    print(f"  {fname}")

def fig_F4():
    _triangulation(["PTB", "LBW", "SGA"],
                   ("effect", "eff_lo", "eff_hi"),   # OR scale (est/lo/hi are log-odds)
                   "Odds ratio per 10 mmHg (log scale)", True, 1.0,
                   "F4_triangulation",
                   "Triangulation: maternal blood pressure and binary perinatal outcomes")

def fig_F4b():
    _triangulation(["BWT"],
                   ("effect_sd", "eff_lo_sd", "eff_hi_sd"),
                   "Change in birth weight (SD) per 10 mmHg", False, 0.0,
                   "F4b_triangulation_bwt",
                   "Triangulation: maternal blood pressure and birth weight")

# =====================================================================
# SF2 -- genetic distance from EUR training vs transferability
# =====================================================================
def fig_SF2():
    d = tsv("SF2_distance")
    if "pc_dist_eur" not in d.columns:
        print("  SF2_distance_vs_r2 SKIPPED (no pc_dist_eur column; needs PCs) — keeping existing PDF")
        return
    d["cohort"] = d["cohort"].replace(CANON)
    d["short"] = d["cohort"].map(T.cohort_short)
    d["dist"], d["R2"] = num(d["pc_dist_eur"]), num(d["R2pct"])
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.3), sharey=True)
    for ax, tr in zip(axes, ["SBP", "DBP"]):
        sub = d[d["trait"] == tr]
        for _, r in sub.iterrows():
            col = T.ANC.get(r["ancestry"], T.OI["grey"])
            ax.scatter(r["dist"], r["R2"], s=90, color=col, edgecolor="white",
                       linewidth=1.2, zorder=3)
            ax.annotate(r["short"], (r["dist"], r["R2"]), xytext=(6, 4),
                        textcoords="offset points", fontsize=8.5)
        ax.set_title(tr)
        ax.set_xlabel("Genetic distance from European training sample\n(PC distance)")
        ax.grid(color="#EEEEEE")
    axes[0].set_ylabel("Incremental $R^2$ (%)")
    # annotate the two same-ancestry Bangladeshi cohorts
    handles = [Line2D([0], [0], marker="o", lw=0, markersize=9, markerfacecolor=T.ANC[a],
               markeredgecolor="white", label={"SAS": "South Asian", "AFR": "African"}[a])
               for a in ["SAS", "AFR"]]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.99, 0.99), title="Cohort ancestry")
    fig.tight_layout()
    T.save(fig, "SF2_distance_vs_r2", OUT)
    print("  SF2_distance_vs_r2")

# =====================================================================
# S6 -- robustness across the eight blood-pressure definitions (SAS)
# =====================================================================
def fig_S6():
    d = tsv("S6_definition_sensitivity")
    d = d[d["stratum"] == "SAS"].copy()
    present = set(d["definition"])
    defs = [x for x in ["resid", "mean", "tri1", "tri2", "tri3", "last", "lt20", "ge20"] if x in present]
    deflab = {"resid": "GA-resid", "mean": "mean", "tri1": "T1", "tri2": "T2",
              "tri3": "T3", "last": "last", "lt20": "<20wk", "ge20": "≥20wk"}
    outs = ["BWT", "PTB", "LBW", "SGA"]
    fig, axes = plt.subplots(1, 4, figsize=(12, 3.8))
    for ax, oc in zip(axes, outs):
        sub = d[d["outcome"] == oc]
        binary = oc != "BWT"
        null = 1.0 if binary else 0.0
        for i, dd in enumerate(defs):
            rr = sub[sub["definition"] == dd]
            if not len(rr):
                continue
            rr = rr.iloc[0]
            if binary:
                est, lo, hi = num(pd.Series([rr["effect"]]))[0], np.nan, np.nan
                th, se = num(pd.Series([rr["theta_per10"]]))[0], num(pd.Series([rr["se"]]))[0]
                lo, hi = np.exp(th - 1.96 * se), np.exp(th + 1.96 * se)
            else:
                est = num(pd.Series([rr["theta_per10"]]))[0]
                lo, hi = num(pd.Series([rr["lo"]]))[0], num(pd.Series([rr["hi"]]))[0]
            col = T.OI["orange"] if dd == "resid" else T.OI["blue"]
            ax.plot([lo, hi], [i, i], color=col, lw=1.5, zorder=2)
            ax.scatter([est], [i], s=42, color=col, edgecolor="white", zorder=3)
        ax.axvline(null, color="#555555", ls="--", lw=1)
        ax.set_yticks(range(len(defs)))
        ax.set_yticklabels([deflab[x] for x in defs] if ax is axes[0] else [])
        ax.set_ylim(-0.6, len(defs) - 0.4)
        ax.invert_yaxis()
        ax.set_title(oc)
        ax.set_xlabel("OR / 10 mmHg" if binary else "g / 10 mmHg")
        ax.grid(axis="x", color="#EEEEEE"); ax.grid(axis="y", visible=False)
    axes[0].set_ylabel("Blood-pressure definition")
    fig.tight_layout()
    T.save(fig, "S6_definition_spread", OUT)
    print("  S6_definition_spread")

# =====================================================================
# SF4 -- per-cohort forest of the MR estimates (binary perinatal + BWT)
# =====================================================================
def _forest(outcomes, to_or, xlabel, logx, null, fname, title, panel="perinatal",
            col_titles=None):
    d = tsv("SF4_forest_data")
    d = d[d["panel"] == panel].copy()
    d["short"] = d["label"].map(lambda x: T.cohort_short(x))
    d["theta"], d["lo"], d["hi"] = num(d["theta"]), num(d["lo"]), num(d["hi"])
    d["pooled"] = d["is_pooled"].str.upper().eq("TRUE")
    traits = ["SBP", "DBP"]
    fig, axes = plt.subplots(len(outcomes), 2, figsize=(10.5, 2.1 * len(outcomes) + 1),
                             sharex="col", squeeze=False)
    for r, oc in enumerate(outcomes):
        for c, tr in enumerate(traits):
            ax = axes[r][c]
            sub = d[(d["outcome"] == oc) & (d["trait"] == tr)].copy()
            # cohorts first (by canonical order), pooled last
            sub["ordk"] = sub.apply(lambda x: (99 if x["pooled"] else 0) +
                                    (T.COHORT_ORDER.index(x["short"]) if x["short"] in T.COHORT_ORDER else 50), axis=1)
            sub = sub.sort_values("ordk")
            labs = []
            for y, (_, rr) in enumerate(sub.iterrows()):
                yy = len(sub) - 1 - y
                est = np.exp(rr["theta"]) if to_or else rr["theta"]
                lo = np.exp(rr["lo"]) if to_or else rr["lo"]
                hi = np.exp(rr["hi"]) if to_or else rr["hi"]
                col = T.OI["black"] if rr["pooled"] else T.OI["blue"]
                mk = "D" if rr["pooled"] else "o"
                ax.plot([lo, hi], [yy, yy], color=col, lw=1.5, zorder=2)
                ax.scatter([est], [yy], s=54 if rr["pooled"] else 40, marker=mk,
                           color=col, edgecolor="white", zorder=3)
                labs.append(rr["label"].strip() if rr["pooled"] else rr["short"])
            ax.axvline(null, color="#555555", ls="--", lw=1)
            if logx:
                ax.set_xscale("log")
            ax.set_yticks(range(len(labs)))
            ax.set_yticklabels(list(reversed(labs)) if c == 0 else [])
            ax.set_ylim(-0.6, len(labs) - 0.4)
            if r == 0:
                ax.set_title(tr)
            if c == 1:
                ax.text(1.03, 0.5, oc, transform=ax.transAxes, rotation=-90,
                        va="center", fontweight="bold")
            ax.grid(axis="x", color="#EEEEEE"); ax.grid(axis="y", visible=False)
    for c in range(2):
        axes[-1][c].set_xlabel(xlabel)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    T.save(fig, fname, OUT)
    print(f"  {fname}")

def fig_SF4():
    _forest(["PTB", "LBW", "SGA"], True, "Odds ratio per 10 mmHg (log scale)", True, 1.0,
            "SF4_forest_mr", "Per-cohort Mendelian-randomization estimates (binary outcomes)")

def fig_SF4c():
    _forest(["BWT"], False, "Change in birth weight (g) per 10 mmHg", False, 0.0,
            "SF4c_forest_bwt", "Per-cohort Mendelian-randomization estimates (birth weight)")

# =====================================================================
# Positive-control forest: score -> preeclampsia / chronic hypertension
# =====================================================================
def fig_poscontrol():
    d = tsv("S12_controls")
    d = d[d["role"] == "positive"].copy()
    if not len(d):
        print("  poscontrol_forest SKIPPED (no positive-control rows)"); return
    d = d[d["trait"] == "SBP"]                       # the ancestry-matched SBP instrument
    d["b"], d["se"] = num(d["beta_pcadj"]), num(d["se_pcadj"])
    d["OR"] = np.exp(d["b"]); d["lo"] = np.exp(d["b"] - 1.96 * d["se"]); d["hi"] = np.exp(d["b"] + 1.96 * d["se"])
    d["pooled"] = d["cohort"].eq("POOLED")
    d["short"] = d.apply(lambda r: "Pooled" if r["pooled"] else T.cohort_short(CANON.get(r["cohort"], r["cohort"])), axis=1)
    conds = [("PE", "Preeclampsia / HDP"), ("CHRON_HTN", "Chronic hypertension")]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), squeeze=False)
    for c, (key, title) in enumerate(conds):
        ax = axes[0][c]
        sub = d[d["variable"] == key].copy()
        sub["ordk"] = sub.apply(lambda x: 99 if x["pooled"] else
                     (T.COHORT_ORDER.index(x["short"]) if x["short"] in T.COHORT_ORDER else 50), axis=1)
        sub = sub.sort_values("ordk"); labs = []
        for y, (_, r) in enumerate(sub.iterrows()):
            yy = len(sub) - 1 - y
            if not np.isfinite(r["OR"]): continue
            col = T.OI["black"] if r["pooled"] else T.OI["green"]
            mk = "D" if r["pooled"] else "o"
            ax.plot([r["lo"], r["hi"]], [yy, yy], color=col, lw=1.5, zorder=2)
            ax.scatter([r["OR"]], [yy], s=54 if r["pooled"] else 40, marker=mk, color=col, edgecolor="white", zorder=3)
            labs.append(r["short"])
        ax.axvline(1.0, color="#555555", ls="--", lw=1); ax.set_xscale("log")
        ax.set_yticks(range(len(labs))); ax.set_yticklabels(list(reversed(labs)) if c == 0 else [])
        ax.set_ylim(-0.6, len(labs) - 0.4); ax.set_title(title); ax.set_xlabel("OR per SD of score")
        ax.grid(axis="x", color="#EEEEEE"); ax.grid(axis="y", visible=False)
    fig.tight_layout()
    T.save(fig, "poscontrol_forest", OUT)
    print("  poscontrol_forest")

def fig_SF4b():
    _forest(["spontaneous", "indicated"], True, "Odds ratio per 10 mmHg (log scale)", True, 1.0,
            "SF4b_forest_subtype",
            "Per-cohort estimates by preterm-birth subtype", panel="PTB subtype")

# =====================================================================
# SF4a -- observational (adjusted) associations per cohort, from S9
# =====================================================================
def fig_SF4a():
    d = tsv("S9_confounding")
    d = d[(d["model"] == "M3_plus_EDU") & (d["sample"] == "own")].copy()
    d["OR"], d["lo"], d["hi"] = num(d["OR"]), num(d["lo"]), num(d["hi"])
    # one row per (trait, outcome, scope): take the first definition present
    d = d.sort_values("definition").drop_duplicates(["trait", "outcome", "scope"])
    d["pooled"] = d["scope"].eq("pooled")
    d["short"] = d["scope"].replace(CANON).map(lambda x: "Pooled" if x == "pooled" else T.cohort_short(x))
    outcomes = ["PTB", "LBW", "SGA"]
    traits = ["SBP", "DBP"]
    fig, axes = plt.subplots(len(outcomes), 2, figsize=(10.5, 2.1 * len(outcomes) + 1),
                             sharex="col", squeeze=False)
    for r, oc in enumerate(outcomes):
        for c, tr in enumerate(traits):
            ax = axes[r][c]
            sub = d[(d["outcome"] == oc) & (d["trait"] == tr)].copy()
            sub["ordk"] = sub.apply(lambda x: (99 if x["pooled"] else
                          (T.COHORT_ORDER.index(x["short"]) if x["short"] in T.COHORT_ORDER else 50)), axis=1)
            sub = sub.sort_values("ordk")
            labs = []
            for y, (_, rr) in enumerate(sub.iterrows()):
                yy = len(sub) - 1 - y
                col = T.OI["black"] if rr["pooled"] else T.OI["blue"]
                mk = "D" if rr["pooled"] else "o"
                ax.plot([rr["lo"], rr["hi"]], [yy, yy], color=col, lw=1.5, zorder=2)
                ax.scatter([rr["OR"]], [yy], s=54 if rr["pooled"] else 40, marker=mk,
                           color=col, edgecolor="white", zorder=3)
                labs.append(rr["short"])
            ax.axvline(1.0, color="#555555", ls="--", lw=1)
            ax.set_xscale("log")
            ax.set_yticks(range(len(labs)))
            ax.set_yticklabels(list(reversed(labs)) if c == 0 else [])
            ax.set_ylim(-0.6, len(labs) - 0.4)
            if r == 0:
                ax.set_title(tr)
            if c == 1:
                ax.text(1.03, 0.5, oc, transform=ax.transAxes, rotation=-90,
                        va="center", fontweight="bold")
            ax.grid(axis="x", color="#EEEEEE"); ax.grid(axis="y", visible=False)
    for c in range(2):
        axes[-1][c].set_xlabel("Odds ratio per 10 mmHg (log scale)")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    T.save(fig, "SF4a_forest_observational", OUT)
    print("  SF4a_forest_observational")

# =====================================================================
# SF1 -- cohorts in PC space (centroid + 1-SD ellipse per cohort)
# =====================================================================
def fig_SF1():
    # warm hues = South Asian cohorts, cool hues = African cohorts
    PAL = {"Sylhet": T.OI["orange"], "Karachi": T.OI["vermillion"], "Matlab": T.OI["purple"],
           "Pemba": T.OI["blue"], "Lusaka": T.OI["skyblue"]}
    pts = os.path.join(TBL, "SF1_pca_points.tsv")
    fig, ax = plt.subplots(figsize=(6.8, 6))
    if os.path.exists(pts):
        # standard PCA: one point per individual, PC1 vs PC2, coloured by cohort
        p = tsv("SF1_pca_points")
        p["short"] = p["cohort"].replace(CANON).map(T.cohort_short)
        p["PC1"], p["PC2"] = num(p["PC1"]), num(p["PC2"])
        for sh in T.COHORT_ORDER:
            s = p[p["short"] == sh]
            if len(s):
                ax.scatter(s["PC1"], s["PC2"], s=6, alpha=0.35, linewidths=0,
                           color=PAL.get(sh, T.OI["grey"]), label=sh, rasterized=True)
        v1 = p["PC1_pct"].iloc[0] if "PC1_pct" in p.columns else None
        v2 = p["PC2_pct"].iloc[0] if "PC2_pct" in p.columns else None
        ax.set_xlabel(f"PC1 ({num(pd.Series([v1]))[0]:.1f}%)" if v1 else "PC1")
        ax.set_ylabel(f"PC2 ({num(pd.Series([v2]))[0]:.1f}%)" if v2 else "PC2")
        leg = ax.legend(title="Cohort", markerscale=3, loc="best", handletextpad=0.3)
    else:
        # fallback until deliv_SF1_pca.R exports per-individual points: clean centroids + 1-SD bars
        d = tsv("SF1_pca_centroids")
        d["short"] = d["cohort"].replace(CANON).map(T.cohort_short)
        for k in ("PC1", "PC2", "PC1_sd", "PC2_sd"):
            d[k] = num(d[k])
        for _, r in d.iterrows():
            col = PAL.get(r["short"], T.OI["grey"])
            ax.errorbar(r["PC1"], r["PC2"], xerr=r["PC1_sd"], yerr=r["PC2_sd"],
                        fmt="o", ms=8, color=col, ecolor=col, elinewidth=1, capsize=3, label=r["short"])
        ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
        ax.legend(title="Cohort", loc="best")
    ax.grid(color="#EEEEEE")
    T.save(fig, "SF1_pca", OUT)
    print("  SF1_pca" + ("" if os.path.exists(pts) else "  [centroid fallback: rerun deliv_SF1_pca.R for per-individual points]"))

# =====================================================================
# SF3 -- BP distributions by cohort (from optional density export)
# =====================================================================
def fig_SF3():
    if not os.path.exists(os.path.join(TBL, "SF3_bp_density.tsv")):
        print("  SF3_bpdist  [awaits SF3 density export from deliv_S16_bpdist.R]")
        return
    de = tsv("SF3_bp_density")
    de["short"] = de["cohort"].replace(CANON).map(T.cohort_short)
    de["x"], de["dens"] = num(de["sbp"]), num(de["density"])
    whichs = list(dict.fromkeys(de["which"])); traits = ["SBP", "DBP"]
    palette = {"Sylhet": T.OI["vermillion"], "Karachi": T.OI["orange"], "Matlab": T.OI["blue"],
               "Pemba": T.OI["green"], "Lusaka": T.OI["purple"]}
    fig, axes = plt.subplots(len(traits), len(whichs),
                             figsize=(4.6 * len(whichs), 3.1 * len(traits)), squeeze=False)
    for r, tr in enumerate(traits):
        for c, wh in enumerate(whichs):
            ax = axes[r][c]
            sub = de[(de["trait"] == tr) & (de["which"] == wh)]
            for sh in T.COHORT_ORDER:
                s = sub[sub["short"] == sh].sort_values("x")
                if len(s):
                    ax.plot(s["x"], s["dens"], color=palette.get(sh, T.OI["grey"]), label=sh, lw=1.3)
            if r == 0:
                ax.set_title(wh)
            if c == 0:
                ax.set_ylabel(f"{tr}\nDensity")
            ax.set_xlabel(f"{tr} (mmHg)")
            ax.grid(color="#EEEEEE")
    axes[0][-1].legend(fontsize=8, loc="upper right")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    T.save(fig, "SF3_bpdist", OUT)
    print("  SF3_bpdist")

# =====================================================================
# F3 -- measurement diagnostic (Sylhet vs Matlab). Panel B always; A/C
# from optional summary exports written by deliv_F3_diagnostic.R.
# =====================================================================
def fig_F3():
    have_ga = os.path.exists(os.path.join(TBL, "F3_bp_by_ga.tsv"))
    have_de = os.path.exists(os.path.join(TBL, "F3_bp_density.tsv"))
    if not (have_ga or have_de):
        # exports not yet generated -> leave the existing (R) 3-panel F3 in place
        print("  F3_diagnostic  [skipped: rerun deliv_F3_diagnostic.R for themed A/B/C panels]")
        return
    npanel = 1 + int(have_ga) + int(have_de)
    fig = plt.figure(figsize=(8.6, 3.0 * npanel + 0.4))
    gs = fig.add_gridspec(npanel, 1, hspace=0.6, top=0.90)
    row = 0
    C = {"Sylhet": T.OI["vermillion"], "Matlab": T.OI["blue"]}

    if have_ga:
        ax = fig.add_subplot(gs[row]); row += 1
        g = tsv("F3_bp_by_ga"); g["short"] = g["cohort"].replace(CANON).map(T.cohort_short)
        g["ga"], g["m"] = num(g["ga_week"]), num(g["mean_sbp"])
        for sh in ["Sylhet", "Matlab"]:
            s = g[g["short"] == sh].sort_values("ga")
            ax.plot(s["ga"], s["m"], "-o", ms=3, color=C[sh], label=sh)
        ax.axvline(20, ls="--", color="#999999"); ax.set_title("A", loc="left", fontweight="bold")
        ax.set_xlabel("Gestational age at reading (weeks)"); ax.set_ylabel("Mean SBP (mmHg)")
        ax.legend()

    # Panel B: R2 by definition
    ax = fig.add_subplot(gs[row]); row += 1
    b = tsv("F3_r2_by_definition"); b["short"] = b["cohort"].replace(CANON).map(T.cohort_short)
    b["R2"], b["lo"], b["hi"] = num(b["R2pct"]), num(b["lo"]), num(b["hi"])
    defs = list(dict.fromkeys(b["definition"]))
    x = np.arange(len(defs)); w = 0.38
    for k, sh in enumerate(["Sylhet", "Matlab"]):
        s = b[b["short"] == sh].set_index("definition").reindex(defs)
        ax.bar(x + (k - 0.5) * w, s["R2"], w, color=C[sh], label=sh,
               yerr=[s["R2"] - s["lo"], s["hi"] - s["R2"]], capsize=2,
               error_kw=dict(lw=0.8, ecolor="#555555"))
    ax.set_xticks(x); ax.set_xticklabels(defs, fontsize=8)
    ax.set_ylabel("Incremental $R^2$ (%)"); ax.legend()
    ax.set_title("B", loc="left", fontweight="bold")
    ax.grid(axis="x", visible=False)

    if have_de:
        ax = fig.add_subplot(gs[row]); row += 1
        de = tsv("F3_bp_density")
        de["short"] = de["cohort"].replace(CANON).map(T.cohort_short)
        de["x"], de["dens"] = num(de["sbp"]), num(de["density"])
        for sh in ["Sylhet", "Matlab"]:
            s = de[de["short"] == sh].sort_values("x")
            ax.fill_between(s["x"], s["dens"], color=C[sh], alpha=0.25)
            ax.plot(s["x"], s["dens"], color=C[sh], label=sh)
        ax.set_title("C", loc="left", fontweight="bold"); ax.set_xlabel("SBP (mmHg)")
        ax.set_ylabel("Density"); ax.legend()

    T.save(fig, "F3_diagnostic", OUT)
    print("  F3_diagnostic" + ("" if (have_ga and have_de) else "  [panels A/C await F3 export]"))

# =====================================================================
# F1 -- participant flow (detailed CONSORT-style, drawn with Graphviz)
# =====================================================================
def fig_F1():
    try:
        import graphviz
        graphviz.Digraph().pipe(format="pdf")   # verify `dot` is installed
    except Exception as e:
        print(f"  F1_flow SKIPPED (graphviz/dot not available: {e}); keeping existing PDF")
        return
    f = tsv("F1_flow_counts")
    c = {str(s).strip(): n for s, n in zip(f["step"], f["n"])}
    def N(k):
        try: return int(c[k])
        except Exception: return None
    def cm(n): return f"{n:,}" if n is not None else "?"
    # per-site breakdown of the transferability sample
    sites = [(T.cohort_short(CANON.get(k.replace("partI_", ""), k.replace("partI_", ""))), N(k))
             for k in c if k.startswith("partI_")]
    sites = [s for s in sites if s[1]]
    sites.sort(key=lambda x: -x[1])
    site_lines = r"\l".join(f"{sh}: {cm(n)}" for sh, n in sites) + r"\l"

    bp = N("of_which_postnatal_BP_only"); nb = N("of_which_no_BP_at_all")
    has_enrol = N("0_enrolled_first_preg") is not None

    g = graphviz.Digraph("flow", format="pdf")
    g.attr(rankdir="TB", nodesep="0.55", ranksep="0.45", splines="ortho",
           fontname="Helvetica", bgcolor="white", pad="0.25")
    g.attr("node", shape="box", style="filled", fillcolor="white", color="black",
           fontname="Helvetica", fontsize="11", margin="0.22,0.15", penwidth="1.1")
    g.attr("edge", color="black", penwidth="1.0", arrowsize="0.7")
    EX = dict(shape="box", fillcolor="white", color="black", fontname="Helvetica",
              fontsize="9.5", margin="0.18,0.11", penwidth="1.0")

    def mainbox(nid, title, key):
        g.node(nid, f'<{title}<br/><b>N = {cm(N(key))}</b>>')
    def exclbox(nid, total_key, reasons):
        rows = "".join(f'{r}<br align="left"/>' for r in reasons)
        g.node(nid, f'<<b>Excluded (N = {cm(N(total_key))})</b><br align="left"/>{rows}>', **EX)

    # main vertical spine (id, title, count-key)
    spine = []
    if has_enrol:
        spine.append(("b0", "Pregnant women enrolled (first pregnancy)", "0_enrolled_first_preg"))
        spine.append(("b1", "Singleton pregnancies", "1_first_preg_mothers"))
    else:
        spine.append(("b1", "Pregnant women enrolled (first pregnancy)", "1_first_preg_mothers"))
    spine.append(("b2", "Genotyped", "2_genotyped"))
    spine.append(("b3", "Transferability sample<br/>(genotyped, &#8805;1 valid antenatal BP)", "3_partI_analytic"))
    spine.append(("b4", "Causal-analysis sample<br/>(non-missing perinatal outcome)", "4_partII_MR"))
    for nid, title, key in spine:
        mainbox(nid, title, key)

    # exclusions branching off each transition (aligned to the spine order)
    excls = []
    if has_enrol:
        excls.append(("excluded_non_singleton", ["Non-singleton (twin / multiple) pregnancies"]))
    excls.append(("excluded_not_genotyped", ["No genotype data"]))
    excls.append(("excluded_no_valid_BP", [f"Postnatal readings only: {cm(bp)}",
                                           f"No blood-pressure reading: {cm(nb)}"]))
    excls.append(("excluded_missing_PTB", ["Missing perinatal outcome"]))

    for i in range(len(spine) - 1):
        j = f"j{i}"
        g.node(j, shape="point", width="0.001", style="invis")
        g.edge(spine[i][0], j, arrowhead="none")
        g.edge(j, spine[i + 1][0])
        e = f"e{i}"
        exclbox(e, excls[i][0], excls[i][1])
        g.edge(j, e)
        with g.subgraph() as s:
            s.attr(rank="same"); s.node(j); s.node(e)

    import tempfile, shutil as _sh
    os.makedirs(OUT, exist_ok=True)
    tmp = tempfile.mkdtemp(); base = os.path.join(tmp, "F1_flow")
    g.render(base, cleanup=True)
    g.format = "png"; g.attr(dpi="300"); g.render(base, cleanup=True)
    for ext in ("pdf", "png"):
        _sh.copy2(base + "." + ext, os.path.join(OUT, "F1_flow." + ext))
    print("  F1_flow")

if __name__ == "__main__":
    fig_F1()
    fig_F2()
    fig_F4()
    fig_F4b()
    fig_SF2()
    fig_S6()
    fig_SF4()
    fig_SF4c()
    fig_SF4b()
    fig_SF4a()
    fig_poscontrol()
    fig_SF1()
    fig_SF3()
    fig_F3()
    print("done ->", OUT)
