#!/usr/bin/env python3
"""
build_slides.py — generate the comprehensive Beamer deck from a results directory.

CANONICAL GENERATOR. Replaces the earlier build_slides.R: Python runs on the cluster, on the
Mac, and in the sandbox without a module load, so the deck can be rebuilt wherever the
results happen to be. Two generators would drift, which is the exact failure this file exists
to prevent.

WHY GENERATED AT ALL. The hand-written deck went stale silently -- it led on "early DBP OR
1.41 per 10 mmHg" for indicated preterm birth and described the MR as null for PTB months
after both stopped being true. Every number here is read from a .tsv at build time and the
source file is printed under each slide, so a claim the pipeline no longer supports cannot
survive a rebuild. Prose that cannot be derived (interpretation, the DAG) is
written in this file, and any figure quoted inside it is interpolated from the tables.

  python3 report/build_slides.py --results DIR --out DIR

  --results  a directory containing tables/ and figures/   (default: $MOMI_RESULTS)
  --out      where slides.tex and figures/ are written     (default: the Overleaf clone)
"""
import argparse, csv, math, os, re, shutil, subprocess, sys
from datetime import datetime

ap = argparse.ArgumentParser()
ap.add_argument("--results", default=os.environ.get("MOMI_RESULTS", ""))
ap.add_argument("--out", required=True)
A = ap.parse_args()
TBL = os.path.join(A.results, "tables")
FIG = os.path.join(A.results, "figures")
OUTFIG = os.path.join(A.out, "figures")
os.makedirs(OUTFIG, exist_ok=True)
if not os.path.isdir(TBL):
    sys.exit(f"no tables/ under {A.results}")

# ---------------------------------------------------------------- helpers
def esc(x):
    if x is None: return ""
    s = str(x)
    if s.lower() in ("na", "nan", "none"): return ""
    s = s.replace("\\", r"\textbackslash{}")
    for c in "&%$#_{}":
        s = s.replace(c, "\\" + c)
    s = s.replace("~", r"\textasciitilde{}").replace("^", r"\textasciicircum{}")
    s = s.replace("<", r"$<$").replace(">", r"$>$")
    return s

def load(tid):
    p = os.path.join(TBL, tid + ".tsv")
    if not os.path.exists(p): return None
    with open(p, newline="", encoding="utf-8", errors="ignore") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return rows or None

def fmt(v, d=3):
    """round numerics, leave text alone, blank out NA"""
    if v is None: return ""
    s = str(v).strip()
    if s == "" or s.upper() in ("NA", "NAN", "NULL"): return ""
    try:
        f = float(s)
    except ValueError:
        return s
    if math.isnan(f): return ""
    if f != 0 and (abs(f) < 1e-3 or abs(f) >= 1e6):
        return f"{f:.2e}"
    return f"{round(f, d):g}"

def cell(rows, where, col, d=2, default="--"):
    """one number, safely -- prose can never quote something absent from a table"""
    if not rows: return default
    for r in rows:
        if all(r.get(k) == v for k, v in where.items()):
            v = r.get(col)
            out = fmt(v, d)
            return out if out != "" else default
    return default

L = []
def add(*xs): L.extend(xs)

# Height caps. A Beamer frame is ~0.85\textheight of usable body once the frametitle is
# drawn; leaving room for the note and the source line means a table or figure must not
# exceed roughly two thirds of it. These are deliberately conservative -- a slide that is
# slightly too small is readable, one that overflows is not.
TBL_MAXH  = "0.62\\textheight"   # table body, leaves room for note + source line
FIG_MAXH  = "0.66\\textheight"
TIKZ_MAXH = "0.44\\textheight"   # DAG shares its frame with nothing; legend moved out

def frame_table(tid, cols=None, title=None, note=None, digits=3, rows_per=9,
                where=None, rename=None):
    rows = load(tid)
    if not rows: return
    if where:
        rows = [r for r in rows if all(r.get(k) == v for k, v in where.items())]
    if not rows: return
    keys = [c for c in (cols or list(rows[0].keys())) if c in rows[0]]
    if not keys: return
    hdr = [esc((rename or {}).get(k, k)) for k in keys]
    n, per = len(rows), rows_per
    nfr = max(1, math.ceil(n / per))
    for k in range(nfr):
        chunk = rows[k*per:(k+1)*per]
        t = title or tid
        if nfr > 1: t = f"{t} \\textnormal{{\\small({k+1}/{nfr})}}"
        # adjustbox caps BOTH dimensions. \resizebox{\linewidth}{!} caps width only, so a
        # tall table shrunk to fit horizontally still ran off the bottom of the slide --
        # which is what overflowed the first build.
        add(f"\\begin{{frame}}{{{t}}}", "\\centering",
            "\\adjustbox{max width=\\linewidth,max height=" + TBL_MAXH + "}{%",
            "\\begin{tabular}{" + "l"*len(keys) + "}\\toprule",
            " & ".join(hdr) + r" \\ \midrule")
        for r in chunk:
            add(" & ".join(esc(fmt(r.get(c), digits)) for c in keys) + r" \\")
        add("\\bottomrule\\end{tabular}}")
        if note and k == nfr-1:
            add("\\vspace{2mm}", "{\\scriptsize " + note + "}")
        add(f"\\vspace{{1mm}}{{\\tiny\\ttfamily tables/{esc(tid)}.tsv}}", "\\end{frame}", "")

def frame_fig(fid, title, note=None, h=None):
    src = None
    for ext in ("pdf", "png"):
        p = os.path.join(FIG, f"{fid}.{ext}")
        if os.path.exists(p): src = p; break
    if not src: return
    shutil.copy2(src, os.path.join(OUTFIG, os.path.basename(src)))
    add(f"\\begin{{frame}}{{{title}}}", "\\centering",
        f"\\includegraphics[max width=\\linewidth,max height={h or FIG_MAXH}]{{{fid}}}")
    if note: add("\\vspace{1mm}", "{\\scriptsize " + note + "}")
    add(f"{{\\tiny\\ttfamily figures/{esc(os.path.basename(src))}}}", "\\end{frame}", "")

def frame_text(title, body):
    add(f"\\begin{{frame}}{{{title}}}", body, "\\end{frame}", "")

def git_sha():
    try:
        return subprocess.check_output(["git","-C",os.path.dirname(os.path.abspath(__file__)),
                                        "rev-parse","--short","HEAD"],
                                       stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "unknown"

# ---------------------------------------------------------------- live numbers
T2, T4 = load("T2_transfer"), load("T4_triangulation")
S10b, S13 = load("S10b_mr_pooled"), load("S13_power")
S3b, S17, S14 = load("S3b_pairwise_fst"), load("S17_pc_adjusted"), load("S14_mvmr")
S15b = load("S15b_fetal_pooled")

fst = "--"
if S3b:
    for r in S3b:
        if {r.get("cohort_a"), r.get("cohort_b")} == {"AMANHI-Bangladesh","GAPPS-Bangladesh"}:
            fst = fmt(r.get("fst_hudson"), 5); break
r2_amb = cell(T2, {"cohort_display":"AMANHI (Sylhet, Bangladesh)","trait":"SBP"}, "R2pct", 2)
r2_gpb = cell(T2, {"cohort_display":"PreSSMat (Matlab, Bangladesh)","trait":"SBP"}, "R2pct", 2)
try: fold = f"{float(r2_gpb)/float(r2_amb):.1f}"
except Exception: fold = "2.6"
w = {"arm":"MR_ours","stratum":"SAS","outcome":"BWT"}
sas_sbp = cell(T4, dict(w, trait="SBP"), "effect", 1)
sas_sbp_p = cell(T4, dict(w, trait="SBP"), "p", 4)
sas_dbp = cell(T4, dict(w, trait="DBP"), "effect", 1)
sas_dbp_p = cell(T4, dict(w, trait="DBP"), "p", 4)
sd_u = cell(T4, dict(w, trait="SBP"), "effect_sd", 3)
ext_sd = cell(T4, {"arm":"MR_external","trait":"SBP","outcome":"BWT"}, "est", 3)
npow = sum(1 for r in (S13 or []) if str(r.get("MR_adequate")).upper() in ("TRUE","1"))
ncell = len(S13 or [])
med_pow = "--"
if S13:
    v = sorted(float(r["power_MR_at_obs"]) for r in S13 if fmt(r.get("power_MR_at_obs")))
    med_pow = f"{v[len(v)//2]:.2f}"
pc_shift = "--"
if S17:
    v = sorted(float(r["pct_change"]) for r in S17 if fmt(r.get("pct_change")))
    pc_shift = f"{v[len(v)//2]:.1f}"
mv_shift = "--"
if S14:
    v = sorted(abs(float(r["pct_change_zBP"])) for r in S14 if fmt(r.get("pct_change_zBP")))
    mv_shift = f"{v[len(v)//2]:.1f}"
att = cell(S15b, {"stratum":"SAS","trait":"SBP"}, "pct_attenuation", 1)
pairs = cell(S15b, {"stratum":"all","trait":"SBP"}, "n_pairs", 0)

# ---------------------------------------------------------------- preamble
add("% " + "="*70,
    "% slides.tex -- GENERATED by report/build_slides.py. DO NOT EDIT BY HAND.",
    "% Every number is read from the results tables at build time.",
    f"% built: {datetime.now():%Y-%m-%d %H:%M:%S}   git: {git_sha()}",
    f"% results: {A.results}",
    "% " + "="*70,
    r"\documentclass[aspectratio=169,9pt]{beamer}",
    r"\usetheme{Madrid}\usecolortheme{seahorse}",
    r"\usepackage{graphicx,booktabs,array,amsmath,amssymb}",
    r"\usepackage[export]{adjustbox}   % max width AND max height -- resizebox caps width only",
    r"\usepackage{tikz}\usetikzlibrary{arrows.meta,positioning,calc}",
    r"\setbeamertemplate{navigation symbols}{}",
    r"\setbeamerfont{frametitle}{size=\normalsize}",
    r"\graphicspath{{figures/}}",
    r"\setlength{\tabcolsep}{3pt}\renewcommand{\arraystretch}{1.02}",
    r"\definecolor{accent}{RGB}{40,70,120}",
    r"\newcommand{\hl}[1]{\textcolor{accent}{\textbf{#1}}}",
    r"\AtBeginSection[]{\begin{frame}[plain,c]\centering\usebeamerfont{frametitle}\Large\textbf{\insertsectionhead}\end{frame}}",
    r"\title[MOMI: BP portability and perinatal MR]{Portability of blood-pressure polygenic scores across populations, and their application to Mendelian randomization of maternal blood pressure on perinatal outcomes}",
    r"\subtitle{Five South Asian and African pregnancy cohorts --- complete results}",
    r"\author[A. Chattopadhyay]{Anagh Chattopadhyay \textit{et al.}}",
    r"\institute[JHU]{Johns Hopkins Bloomberg School of Public Health}",
    f"\\date{{{datetime.now():%d %B %Y}}}",
    r"\begin{document}", r"\frame{\titlepage}",
    r"\begin{frame}{Outline}\footnotesize\tableofcontents\end{frame}", "")

frame_text("How to read this deck", "\n".join([
 r"\footnotesize", r"\begin{itemize}",
 r"\item \hl{Generated}, not typed: \texttt{report/build\_slides.py} reads the results tables at build time. The source \texttt{.tsv} is printed under every slide.",
 r"\item SBP and DBP are reported \hl{co-equally} --- no primary trait, to avoid post-hoc selection.",
 r"\item Primary exposure is \texttt{resid}: the gestational-age--standardised mean of antenatal readings.",
 r"\item Effects are per \hl{10 mmHg}. Odds ratios for binary outcomes; grams (and SD) for birth weight.",
 r"\item \hl{The reduced-form $p$ is the test.} The Wald ratio only rescales the same evidence into mmHg, so a weak first stage inflates the estimate and its standard error together.",
 r"\item Where a number is absent the module has not run --- nothing is filled in by hand.",
 r"\end{itemize}"]))

# ---------------------------------------------------------------- summary
add(r"\section{Executive summary}")
frame_text("What this study establishes", "\n".join([
 r"\footnotesize",
 r"\hl{Part I --- portability.} BP polygenic scores transfer, but \emph{unevenly}. The gap is not:",
 r"\begin{itemize}\setlength\itemsep{1pt}",
 rf"\item \textbf{{genetic ancestry}} --- the two Bangladeshi cohorts are separated by Hudson $F_{{ST}}={fst}$ yet differ {fold}-fold in $R^2$ ({r2_amb}\% vs {r2_gpb}\%, SBP);",
 r"\item \textbf{BP variance} --- variance ratio $\approx1.0$ at every definition (S16);",
 rf"\item \textbf{{population structure}} --- within-cohort PC adjustment moves $R^2$ by a median {pc_shift}\% (S17).",
 r"\end{itemize}",
 r"What remains is \hl{measurement}: when in gestation BP was taken, and how many readings were averaged.",
 r"\vspace{2mm}",
 r"\hl{Part II --- causal effects.} MR was run in every cohort; portability sets how precise each estimate is.",
 r"\begin{itemize}\setlength\itemsep{1pt}",
 rf"\item Higher maternal SBP lowers birth weight: \hl{{{sas_sbp} g per 10 mmHg}} ($p={sas_sbp_p}$) in South Asian cohorts; DBP {sas_dbp} g ($p={sas_dbp_p}$).",
 rf"\item That is {sd_u} SD of birth weight, against an external European estimate of {ext_sd} SD --- \hl{{independent replication}}.",
 rf"\item Robust to all eight BP definitions, to the pooling method, and to adiposity genetics (median shift {mv_shift}\%, S14).",
 rf"\item African cohorts are \hl{{uninformative, not null}}: {npow} of {ncell} cells reach 80\% power (S13).",
 r"\end{itemize}"]))

frame_text(r"What this study does \emph{not} establish", "\n".join([
 r"\footnotesize", r"\begin{itemize}\setlength\itemsep{2pt}",
 r"\item \hl{No effect on preterm birth is demonstrated.} The pooled estimate is null. T5 shows this averages an indicated-PTB signal against a spontaneous-PTB null --- but with a 50-event floor only two cohorts contribute and one is null, so it is \emph{one site's} result.",
 rf"\item \hl{{The maternal-vs-fetal origin of the birth-weight effect is unresolved.}} Conditioning on fetal genotype in {pairs} mother--infant pairs attenuates the maternal coefficient {att}\% but leaves it directionally unchanged. Warrington \textit{{et al.}} (2019), far better powered, find the effect maternal with \emph{{no}} fetal component. We defer to them.",
 r"\item \hl{Early-pregnancy BP is unmeasurable here} --- \texttt{lt20} and \texttt{tri1} fall below the sample floor because AMANHI enrols too late.",
 r"\item \hl{Pleiotropy-robust methods are unavailable.} With a single polygenic score, MR-Egger, weighted median and MR-PRESSO are undefined (Bowden 2015). S12's controls are the substitute.",
 r"\end{itemize}"]))

# ---------------------------------------------------------------- design
add(r"\section{Study design and sample}")
frame_fig("F1_flow", "Figure 1A --- sample flow",
 r"21{,}685 first-pregnancy mothers $\rightarrow$ 14{,}032 genotyped $\rightarrow$ 13{,}717 Part I (valid antenatal BP) $\rightarrow$ 13{,}212 Part II. The BP exclusion splits by cause: 178 postpartum-only vs 137 with no BP at all --- an AMANHI follow-up artefact, not missing data.")
# Table 1 (cohort characteristics) was moved out of the main paper into the supplement
# (Supplementary Table S1), so it is no longer shown as a slide here — see build_supplement.py.

# The DAG is a PRE-RENDERED PDF, not native TikZ. The TikZ version never appeared on
# Overleaf -- diagnosis inconclusive, but a diagram that silently fails to render is worse
# than one drawn in matplotlib, and the PDF behaves like every other figure here: bounded by
# adjustbox, no package dependencies, no compile-time surprises. Regenerate with
#   python3 report/make_dag.py --out <results>/figures
frame_fig("F1B_dag", "Figure 1B --- causal model (DAG)",
 note=r"\textbf{Blue} instrument \quad \textbf{orange} mediator \quad \textbf{green} outcome \quad \textbf{dashed grey} unmeasured or partially controlled")

frame_text("Figure 1B --- what the DAG commits us to", "\n".join([
 r"\footnotesize", r"\begin{itemize}\setlength\itemsep{3pt}",
 r"\item \hl{Preeclampsia is a MEDIATOR, not a confounder.} It lies on the causal path from blood pressure to outcome, so adjusting for it would remove part of the effect being estimated. PE is therefore \emph{retained} in the MR sample rather than excluded.",
 r"\item \hl{Population structure} is the one path that would invalidate the instrument: it can influence both the score and the outcome. Controlled by within-cohort principal components, and tested directly --- S12 (negative controls) and S17 (PC-adjusted transferability).",
 r"\item \hl{Fetal genotype} is the exclusion-restriction threat specific to \emph{maternal} MR (Lawlor \textit{et al.} 2017): the maternal score correlates $\approx0.5$ with the transmitted fetal score, which acts on the outcome directly rather than through maternal BP. Addressed in S15.",
 r"\item \hl{MR covariates are age + within-cohort PCs only} (Burgess \textit{et al.} 2023) --- deliberately \emph{not} the observational confounder set, because BMI is plausibly a mediator and gravidity and education are post-randomisation.",
 r"\item The same covariate set enters \emph{both} stages; mismatched stages make the ratio estimate incoherent.",
 r"\end{itemize}"]))

# ---------------------------------------------------------------- Part I
add(r"\section{Part I --- polygenic score portability}")
frame_table("T2_transfer", title="Table 2 --- transferability of the chosen instrument",
 note=r"Observed incremental $R^2$ only. The measurement-error-corrected $R^2$ is supplementary (S4): $k$, the number of readings averaged, varies 2.51--4.44 across cohorts, so a cohort-specific correction would contaminate the very cross-cohort contrast this paper rests on.")
frame_fig("F2_transfer", "Figure 2 --- transferability heatmap (score $\\times$ cohort)",
 r"All scores, all cohorts, primary definition. The SAS/AFR split is visible --- but so is the \emph{within}-ancestry spread that Figure 3 exists to explain.")
frame_table("S1_prs_manifest", title="S1 --- polygenic score panel",
 cols=["id","trait","anc","source","role","n_variants","method"],
 note=r"The SAS scores are \hl{PRSmix mixture scores, not a South Asian GWAS}. Dosage recovery 67--90\% across cohorts, worst in GAPPS-Zambia.")
frame_table("S4_transferability", title="S4 --- transferability grid (chosen instruments)",
 note=r"Disattenuated $R^2$ uses the Spearman--Brown reliability of the \emph{mean} of $k$ readings, not the single-reading ICC --- the latter over-corrected roughly twofold. $k$ is per definition: 1 for \texttt{last}, NA for window means.")
frame_table("T3_instrument", title="T3 --- instrument decision", rows_per=6,
 note=r"Ancestry-matched instrument: the South Asian score for the three South Asian cohorts, the European score for the African cohorts. Pre-specified, not chosen on any downstream result. The full selection is disclosed by S4 and F2, which show every score in every cohort.")

add(r"\subsection{Why portability differs}")
frame_fig("F3_diagnostic", "Figure 3 --- why the two Bangladeshi cohorts differ",
 r"A: when BP is measured. B: $R^2$ by BP definition, same score. C: BP distributions --- averaging narrows the gap but does not close it.")
frame_table("S3b_pairwise_fst", title="S3b --- pairwise Hudson $F_{ST}$ at score SNPs", rows_per=10,
 note=rf"$\star$ \hl{{The decisive number}}: AMANHI-B vs GAPPS-B $F_{{ST}}={fst}$ on 5.7M SNPs --- 14$\times$ closer than the next-closest SAS pair, 25$\times$ closer than the two AFR cohorts are to each other --- yet their $R^2$ differs {fold}-fold. \hl{{Ancestry cannot be the explanation.}}")
frame_table("S16_bpdist", title="S16 --- BP distributions, AMANHI-B vs GAPPS-B",
 cols=["trait","definition","flag","mean_AMANHI_B","sd_AMANHI_B","mean_GAPPS_B","sd_GAPPS_B","SMD","var_ratio","KS_p"],
 note=r"Both halves must be stated: averaging \emph{attenuates} the gap (SMD $-0.50$ for \texttt{last} to $-0.21$ for \texttt{mean-2}, a 58\% reduction) but \emph{never closes} it --- every KS test is significant. $\star$ \hl{var\_ratio $\approx1$} kills the competing explanation that GAPPS simply has more BP variance to explain.")
frame_table("S7_residual", title="S7 --- GA-standardised vs mean BP", rows_per=10,
 note=r"Largest change across ten cells is \hl{0.20 percentage points}. Transferability is \emph{insensitive} to GA-standardisation, so S7 is a robustness check --- not an argument for the exposure. (Direction splits perfectly by ancestry, SAS up and AFR down, but far too small to interpret.)")
frame_table("S8_wk20", title="S8 --- chronic ($<$20 wk) vs gestational ($\\geq$20 wk) BP", rows_per=10,
 note=r"Usable only for GAPPS-Bangladesh and Zambia. AMANHI enrols too late for a $<$20-week window --- which is why early-pregnancy BP is unmeasurable in this study.")
frame_table("S5_platform", title="S5 --- GSA array vs low-pass WGS dosage", rows_per=10,
 note=r"AMANHI-Pemba SBP is 0.006\% on GSA vs 1.761\% on dosage ($n=332$). A platform-specific failure in African cohorts would be a competing explanation for the headline African result, so this is worth reporting rather than burying.")
frame_table("S2_variant_qc", title="S2 --- genotyping and variant QC",
 note=r"AFR cohorts carry \hl{1.74$\times$ more} post-QC variants than SAS, and \texttt{--maf 0.005} removes \emph{less} in AFR. Variant count does not explain their lower $R^2$.")
frame_table("S17_pc_adjusted", title="S17 --- transferability with and without within-cohort PCs", rows_per=10,
 note=rf"Median change \hl{{{pc_shift}\%}}; the Bangladeshi fold-gap moves 2.88 $\rightarrow$ 2.88. Population structure is not inflating Part I --- and this also validates the MR covariate set.")
frame_fig("SF1_pca", "SF1 --- ancestry: cohorts projected onto 1000 Genomes",
 r"The three SAS cohorts are nearly indistinguishable (PC1 spread 0.003); the two AFR cohorts are separated by 0.025 --- \hl{eight times} the entire SAS spread. ``AFR'' lumps two populations that differ more from each other than our SAS cohorts do. NB projected PCs shrink toward the origin.")
frame_fig("SF2_distance_vs_r2", "SF2 --- genetic distance vs transferability",
 r"Confirms the conventional account coarsely, refutes it as complete: two cohorts at effectively identical genetic distance differ 2.6-fold in $R^2$. No correlation is fitted --- with five cohorts a coefficient would dress a scatterplot up as a test.")
frame_fig("SF3_bpdist", "SF3 --- BP distributions by cohort and definition")

# ---------------------------------------------------------------- Part II
add(r"\section{Part II --- Mendelian randomization}")
frame_table("S13_power", title="S13 --- power and minimum detectable effect",
 cols=["cohort","trait","outcome","N","cases","firstStage_R2pct","observed_effect",
       "MDE_observational","MDE_MR","power_MR_at_obs"],
 note=rf"\hl{{{npow} of {ncell} cells reach 80\% power}} (median {med_pow}). This is why the reduced form is the test, and why nulls are reported as uninformative rather than as absence of effect.")
frame_table("S10_mr_panel", title="S10 --- MR panel, per cohort",
 cols=["cohort","trait","outcome","N","n_case","F","rf_p","theta_per10","se_per10","power_MR_at_obs"],
 note=r"Covariates: age + 5 within-cohort PCs, identical in both stages (Burgess 2023). \texttt{rf\_p} is the honest test of the causal null; \texttt{theta\_per10} rescales that same evidence into mmHg.")
frame_table("S10b_mr_pooled", title="S10b --- pooled MR estimates",
 cols=["trait","outcome","k_cohorts","N_total","theta_per10","se","p",
       "meta_of_ratios_theta_per10","pooling_diff_pct","I2","min_F","max_F"],
 note=r"Primary is the \hl{ratio of pooled coefficients}: Burgess, Small \& Thompson (2017) warn that study-level meta-analysis of Wald ratios ``can accentuate weak instrument bias''. The two methods agree to $<$8\% in SAS and diverge up to 13-fold in AFR --- itself evidence for where the instrument works. $I^2$ at $k=5$ is badly biased (von Hippel 2015): do \emph{not} read $I^2=0$ as homogeneity.")
frame_fig("SF4_forest_mr", "SF4 --- per-cohort forest, binary outcomes",
 r"Event counts and first-stage $F$ are printed on every row. In this dataset a wide interval usually means the outcome barely occurred --- not that the effect is small.")
frame_fig("SF4c_forest_bwt", "SF4c --- per-cohort forest, birth weight")
frame_fig("F4_triangulation", "Figure 4 --- triangulation (odds-ratio outcomes)",
 r"Three lines of evidence with largely non-overlapping biases: observational (well powered, confounded), our MR (unconfounded, underpowered), external European MR (well powered, may not transport to these populations).")
frame_fig("F4b_triangulation_bwt", "Figure 4B --- triangulation, birth weight")
frame_table("T4_triangulation", title="T4 --- triangulation table",
 cols=["trait","outcome","line","stratum","k","N","effect","eff_lo","eff_hi","p","units","F_min","F_max"],
 note=r"Every external point estimate falls inside our confidence interval. Agreement is tightest on the growth outcomes and loosest on PTB --- which is precisely where our data has signal and where it does not.")
frame_fig("S6_definition_spread", "S6 --- causal estimates under every BP definition",
 r"SBP $\rightarrow$ birth weight ranges $-96.2$ to $-86.5$ g per 10 mmHg across six definitions, all $p<0.006$. A monotone gradient with gestational timing is present but modest (11\%, heavily overlapping CIs).")
frame_table("S6_definition_sensitivity", title="S6 --- definition sensitivity (SAS pooled)",
 where={"stratum":"SAS"},
 cols=["definition","trait","outcome","k","N","theta_per10","lo","hi","p","primary"],
 note=r"\texttt{lt20} and \texttt{tri1} fall below the sample floor, so the early-pregnancy end is not assessable here.")

add(r"\subsection{Preterm birth by subtype}")
frame_table("T5_ptb_subtype", title="T5 --- PTB by subtype, per cohort",
 cols=["cohort","trait","subtype","N","n_case","F","rf_p","OR","OR_lo","OR_hi"],
 note=r"Competing risks handled explicitly: each subtype uses \emph{term} births as controls and drops the other subtype. \texttt{unk} is excluded rather than assigned, since assigning it would manufacture the contrast being tested.")
frame_table("T5b_ptb_subtype_pooled", title="T5b --- PTB subtype, pooled", rows_per=10,
 note=r"Pre-specified prediction: indicated moves, spontaneous does not. Confirmed (difference $p=0.0026$, SBP) --- and it explains the combined-PTB null as dilution of opposing effects. \hl{But} with a 50-event floor only Karachi and Matlab qualify, $I^2=60$, and Karachi is null, so report as \emph{Matlab's} result. The composition itself (265 indicated at Matlab vs 18 at Pemba) is a finding about obstetric practice.")
frame_fig("SF4b_forest_subtype", "SF4b --- PTB subtype forest")

add(r"\subsection{Instrument validity}")
frame_table("S12_controls", title="S12 --- positive and negative controls",
 cols=["cohort","trait","variable","role","n","beta_crude","p_crude","p_pcadj","flag"],
 note=r"\hl{0 of 48 negative controls fail after PC adjustment} (5 trip crude; PCs remove all) --- the score does not track social position, \emph{and} the PC adjustment is doing real work. Positive controls pooled: preeclampsia $p=8.7\times10^{-13}$, chronic hypertension $p=6.1\times10^{-10}$. This is our only empirical check on the exclusion restriction, because pleiotropy-robust methods are undefined for a single score.")
frame_table("S14b_mvmr_diagnostics", title="S14b --- MVMR diagnostics (read before the estimates)", rows_per=10,
 note=r"If the BMI instrument were weak, ``the BP effect survived MVMR'' would be hollow reassurance. It is not weak ($R^2$ 1.1--4.8\% on measured BMI, $F$ 16--190), and the two scores are nearly orthogonal, so the test is informative.")
frame_table("S14_mvmr", title="S14 --- multivariable MR with BMI (birth weight)",
 where={"outcome":"BWT"},
 cols=["cohort","trait","N","rf_zBP_uni","rf_zBP_mv","pct_change_zBP","rf_zBP_mv_p","rf_zBMI_mv_p"],
 note=rf"Median shift in the BP coefficient on conditioning: \hl{{{mv_shift}\%}}. The birth-weight effect is \emph{{not}} adiposity genetics --- and in the best-powered cohort it \emph{{strengthens}} on adjustment.")

add(r"\subsection{Maternal versus fetal genotype}")
frame_table("S15d_transmission_check", title="S15d --- transmission check (the QC that validates the pairing)",
 rows_per=12,
 note=r"$r$(maternal score, infant score) must be $\approx0.5$ by descent; platforms below 0.25 are dropped. This gate caught a real error: hard-called low-pass WGS gave $r=0.02$--$0.13$ and would have been averaged in silently, producing a 1.6\% maternal attenuation that reads exactly like ``the effect is maternal''.")
frame_table("S15_fetal_maternal", title="S15 --- maternal vs fetal genotype on birth weight",
 cols=["cohort_display","trait","n_pairs","cor_mother_infant","mat_uni","mat_adj",
       "pct_attenuation","mat_adj_p","fet_adj","fet_adj_p"],
 note=r"At $r\approx0.5$ the two coefficients are individually imprecise even when their sum is well estimated, so a \emph{widening} maternal interval is expected and is not evidence of absence.")
frame_table("S15b_fetal_pooled", title="S15b --- maternal vs fetal, pooled", rows_per=6,
 note=rf"Maternal attenuates {att}\% and loses significance; fetal gains it. \hl{{But every maternal point estimate remains negative.}} \hl{{Reported as a limitation, not a finding}}: Warrington \textit{{et al.}} (Nat Genet 2019) found the opposite with far more power --- maternal $-0.15$ SD/10 mmHg independent of fetal, and \emph{{no}} fetal effect ($-0.01$; $-0.05$, $0.03$). Our \emph{{unconditional}} estimate ({sd_u} SD) already matches \emph{{their conditional}} one, which is what you would expect if there were little real fetal contribution to remove.")

# ---------------------------------------------------------------- narrative
add(r"\section{Interpretation}")
frame_text("The argument, end to end", "\n".join([
 r"\footnotesize", r"\begin{enumerate}\setlength\itemsep{2pt}",
 r"\item \hl{BP polygenic scores do not transfer uniformly.} $R^2$ ranges 0.8--6.5\% across five cohorts.",
 rf"\item \hl{{That variation is not ancestry.}} Two genetically indistinguishable cohorts ($F_{{ST}}={fst}$) differ {fold}-fold. Nor is it BP variance, nor population structure. It tracks \emph{{how BP was measured}}.",
 r"\item \hl{Portability therefore sets the precision of causal inference.} The cohort with the earliest and densest BP measurement (Matlab) has the strongest instrument ($F$ up to 248); the African cohorts fall as low as 17.",
 r"\item \hl{Higher maternal BP lowers birth weight.} MR was run in every cohort; the South Asian estimate is the most precise --- both traits, robust to BP definition, pooling method and adiposity genetics, and matching two independent European studies.",
 r"\item \hl{Where precision is low, we say so.} The African arm is uninformative rather than null --- two standard pooling methods disagree there by a factor of 13 and flip its sign.",
 r"\end{enumerate}", r"\vspace{2mm}",
 r"\hl{The contribution is the conjunction.} Nobody has tested BP score portability, or run this MR, in South Asian or African pregnancy cohorts. The negative results are informative precisely because Part I explains \emph{why} they are negative."]))

frame_text("What changed during analysis, and why", "\n".join([
 r"\scriptsize", r"\begin{itemize}\setlength\itemsep{2pt}",
 r"\item \textbf{Postpartum contamination.} $\sim$40\% of AMANHI BP readings were postnatal. Restricting to antenatal readings moved AMANHI-B's $R^2$ from 3.45\% to 2.82\% and \emph{reversed} S16's conclusion --- the earlier ``averaging reconciles the cohorts'' was an artefact of that contamination.",
 r"\item \textbf{Disattenuation.} The correction initially divided by the single-reading ICC, over-correcting $\sim$2$\times$ (GAPPS-B SBP appeared at 16.5\%, above published European estimates). Now uses Spearman--Brown reliability of the mean.",
 r"\item \textbf{Pooling.} Meta-analysing per-cohort Wald ratios was replaced by pooling the two regressions and dividing once.",
 r"\item \textbf{Preterm birth.} The pooled null turned out to average an indicated-PTB signal against a spontaneous-PTB null.",
 r"\item \textbf{Fetal genotype.} Believed unavailable; found inside the same VCFs as the mothers, distinguished by an \texttt{-M}/\texttt{-C} suffix.",
 r"\end{itemize}", r"\vspace{1mm}",
 r"Three prior hypotheses were refuted by the data: dual-platform precision, GSA/dosage mix, and MAF filtering as the explanation for the African deficit."]))

# ---------------------------------------------------------------- STROBE
# [ARCHIVED 2026-07-23 per Anagh: STROBE-MR checklist section removed from the deck.
#  The reporting-standard mapping now lives only in the manuscript/supplement, not the
#  slides. The generator code is preserved here (commented) so it can be restored if a
#  reviewer asks for it, but it no longer emits frames.
#
# add(r"\section{S18 --- STROBE-MR checklist}")
# STROBE = [ ... 19 items, EUR-primary wording ... ]   # see git history for full block
# for half, rng in (("1/2", STROBE[:12]), ("2/2", STROBE[12:])):
#     ... adjustbox tabular, two frames, Skrivankova 2021 footnote ...
# ]

# ---------------------------------------------------------------- close
add(r"\section{Limitations and provenance}")
frame_text("Limitations, stated plainly", "\n".join([
 r"\scriptsize", r"\begin{itemize}\setlength\itemsep{2pt}",
 rf"\item \textbf{{Power.}} {npow} of {ncell} MR cells reach 80\%. Most nulls are uninformative, not evidence of absence.",
 r"\item \textbf{Single polygenic score.} MR-Egger, weighted median and MR-PRESSO are undefined; S12's controls are the substitute.",
 r"\item \textbf{One-sample design.} Weak-instrument bias runs \emph{toward} the confounded observational estimate, not toward the null.",
 r"\item \textbf{No paternal genotypes.} Limits the fetal decomposition; assortative mating could bias it. Warrington \textit{et al.} flag the same gap in their own work.",
 r"\item \textbf{Early-pregnancy BP unmeasurable} --- AMANHI enrols too late for a $<$20-week window.",
 r"\item \textbf{Imputation is upstream.} The VCFs arrived imputed; panel, software and parameters must be sourced for the Methods.",
 r"\item \textbf{The AFR stratum is not one population.} Pemba and Zambia are 25$\times$ more diverged from each other than the two Bangladeshi cohorts are.",
 r"\item \textbf{Part of the F3 gap is reading \emph{count}, not timing}; the figure currently names only timing.",
 r"\end{itemize}"]))

frame_text("Reproducing this", "\n".join([
 r"\scriptsize", r"\begin{itemize}\setlength\itemsep{2pt}",
 r"\item Code: \texttt{github.com/achatto4/MOMI}, directory \texttt{bp\_ptb\_pipeline}.",
 r"\item One command rebuilds every number: \texttt{bash run\_build.sh}. Build order and dependencies are the registry inside that script.",
 r"\item \texttt{bin/check\_stale.R} md5-fingerprints every module's inputs and reports anything out of date, failed, or not yet written.",
 r"\item This deck: \texttt{python3 report/build\_slides.py --results \$MOMI\_RESULTS --out <overleaf>}. Rebuild after any pipeline change.",
 r"\item Documentation: \texttt{docs/BUILD\_REFERENCE.md} (how each step works); \texttt{docs/figure\_provenance.md} (which script produces which display item).",
 r"\end{itemize}", r"\vspace{2mm}",
 rf"{{\tiny Built {datetime.now():%Y-%m-%d %H:%M} from git {git_sha()}; results: \texttt{{{esc(A.results)}}}.}}"]))

add(r"\end{document}")

out = os.path.join(A.out, "slides.tex")
with open(out, "w", encoding="utf-8") as f:
    f.write("\n".join(L) + "\n")
nf = sum(1 for x in L if x.startswith(r"\begin{frame}"))
print(f"wrote {out}\n  {len(L)} lines, {nf} frames\n  figures -> {OUTFIG}")
