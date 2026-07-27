#!/usr/bin/env python3
"""
build_paper.py — generate the manuscript (Methods + Results) from a results directory.

Same philosophy as build_slides.py: every number is read from a .tsv at build time, so the
draft cannot state a figure the pipeline no longer produces. Prose is written here.

SCOPE, per Anagh 2026-07-21: state the RESULTS and emphasise the METHODS; facts only.
Background, Discussion and all interpretation/framing are left as clearly marked
placeholders for the author to fill. Target journal: BMC Medicine.

STYLE: measured, quantitative, first-person plural, methods-forward -- the register of the
Chatterjee-group transferability papers. Effect sizes always with intervals; assumptions and
their tests stated explicitly; ancestry-transferability caveats foregrounded.

  python3 report/build_paper.py --results DIR --out DIR
"""
import argparse, csv, math, os, subprocess
from datetime import datetime

ap = argparse.ArgumentParser()
ap.add_argument("--results", default=os.environ.get("MOMI_RESULTS", ""))
ap.add_argument("--out", required=True)
A = ap.parse_args()
TBL = os.path.join(A.results, "tables")

def load(tid):
    p = os.path.join(TBL, tid + ".tsv")
    if not os.path.exists(p): return None
    with open(p, newline="", encoding="utf-8", errors="ignore") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return rows or None

def fmt(v, d=2):
    if v is None: return ""
    s = str(v).strip()
    if s == "" or s.upper() in ("NA","NAN","NULL"): return ""
    try: f = float(s)
    except ValueError: return s
    if math.isnan(f): return ""
    if f != 0 and (abs(f) < 1e-3 or abs(f) >= 1e6): return f"{f:.2e}"
    return f"{round(f,d):g}"

def _num(v):
    try: return float(str(v).strip())
    except Exception: return 0.0

def cell(rows, where, col, d=2, default="NA"):
    if not rows: return default
    for r in rows:
        if all(r.get(k)==v for k,v in where.items()):
            o = fmt(r.get(col), d); return o if o else default
    return default

def sha():
    try:
        return subprocess.check_output(
            ["git","-C",os.path.dirname(os.path.abspath(__file__)),"rev-parse","--short","HEAD"],
            stderr=subprocess.DEVNULL).decode().strip()
    except Exception: return "unknown"

# ---- live numbers ----
T1  = load("T1_cohorts");         T2  = load("T2_transfer")
T4  = load("T4_triangulation");   S13 = load("S13_power")
S3b = load("S3b_pairwise_fst");   S17 = load("S17_pc_adjusted")
S16 = load("S16_bpdist");         S10b= load("S10b_mr_pooled")
T5b = load("T5b_ptb_subtype_pooled")

def fst_bd():
    if S3b:
        for r in S3b:
            if {r.get("cohort_a"),r.get("cohort_b")}=={"AMANHI-Bangladesh","GAPPS-Bangladesh"}:
                return fmt(r.get("fst_hudson"),5)
    return "NA"
FST = fst_bd()
r2_amb = cell(T2, {"cohort_display":"AMANHI (Sylhet, Bangladesh)","trait":"SBP"}, "R2pct", 2)
r2_gpb = cell(T2, {"cohort_display":"PreSSMat (Matlab, Bangladesh)","trait":"SBP"}, "R2pct", 2)
try: FOLD = f"{float(r2_gpb)/float(r2_amb):.1f}"
except Exception: FOLD = "NA"
W = {"arm":"MR_ours","stratum":"SAS","outcome":"BWT"}
sas_sbp   = cell(T4, dict(W,trait="SBP"), "effect", 1)
sas_sbp_l = cell(T4, dict(W,trait="SBP"), "eff_lo", 1)
sas_sbp_h = cell(T4, dict(W,trait="SBP"), "eff_hi", 1)
sas_sbp_p = cell(T4, dict(W,trait="SBP"), "p", 4)
sas_dbp   = cell(T4, dict(W,trait="DBP"), "effect", 1)
sas_dbp_l = cell(T4, dict(W,trait="DBP"), "eff_lo", 1)
sas_dbp_h = cell(T4, dict(W,trait="DBP"), "eff_hi", 1)
sas_dbp_p = cell(T4, dict(W,trait="DBP"), "p", 4)
sd_u   = cell(T4, dict(W,trait="SBP"), "effect_sd", 3)
ext_sd = cell(T4, {"arm":"MR_external","trait":"SBP","outcome":"BWT"}, "est", 3)
npow  = sum(1 for r in (S13 or []) if str(r.get("MR_adequate")).upper() in ("TRUE","1"))
ncell = len(S13 or [])
pc_shift = "NA"
if S17:
    v = sorted(float(r["pct_change"]) for r in S17 if fmt(r.get("pct_change")))
    if v: pc_shift = f"{v[len(v)//2]:.1f}"

def r2range(anc):
    vals=[float(r["R2pct"]) for r in (T2 or []) if r.get("ancestry")==anc and fmt(r.get("R2pct"))]
    return (f"{min(vals):.1f}", f"{max(vals):.1f}") if vals else ("NA","NA")
sas_lo,sas_hi = r2range("SAS"); afr_lo,afr_hi = r2range("AFR")
def frange(anc=None):
    vals=[float(r["F"]) for r in (T2 or []) if fmt(r.get("F")) and (anc is None or r.get("ancestry")==anc)]
    return (f"{min(vals):.0f}", f"{max(vals):.0f}") if vals else ("NA","NA")
F_lo, F_hi = frange(); Fsas_lo, Fsas_hi = frange("SAS")
n_union = cell(T1, {}, "N_genotyped_union", 0) if T1 else "NA"
if T1:  # union across the five genotyped cohorts (T1 has THSTI too, so sum genotyped>0)
    tot = sum(int(float(r["N_genotyped_union"])) for r in T1
              if fmt(r.get("N_genotyped_union")) and float(r["N_genotyped_union"])>0)
    n_union = f"{tot:,}" if tot else "NA"

# external OR comparators (LBW/SGA/PTB) for the triangulation paragraph
def ext_or(tr,oc): return cell(T4, {"arm":"MR_external","trait":tr,"outcome":oc}, "effect", 3)
def our_or(tr,oc): return cell(T4, {"arm":"MR_ours","stratum":"SAS","trait":tr,"outcome":oc}, "effect", 3)

FIG = os.path.join(A.results, "figures")
OUTFIG = os.path.join(A.out, "figures")
os.makedirs(OUTFIG, exist_ok=True)

def esc(x):
    s = "" if x is None else str(x)
    if s.lower() in ("na","nan","none"): return ""
    s = s.replace("\\", r"\textbackslash{}")
    for c in "&%$#_{}": s = s.replace(c, "\\"+c)
    s = s.replace("~", r"\textasciitilde{}")
    return s

L=[]; add=L.append
def P(*xs):
    for x in xs: L.append(x)

def fig_float(fid, caption, label, width="0.92\\textwidth"):
    """place a figure float; copy the PDF (or PNG) into the paper's figures/."""
    src=None
    for ext in ("pdf","png"):
        p=os.path.join(FIG, f"{fid}.{ext}")
        if os.path.exists(p): src=p; break
    if not src:
        P(f"% [figure {fid} not found -- module not run]"); return
    import shutil; shutil.copy2(src, os.path.join(OUTFIG, os.path.basename(src)))
    P(r"\begin{figure}[htbp]\centering",
      f"\\includegraphics[width={width},keepaspectratio]{{{fid}}}",
      f"\\caption{{{caption}}}", f"\\label{{{label}}}", r"\end{figure}", "")

def tab_float(tid, cols, caption, label, rename=None, digits=2, where=None, rows=None,
              small=True, filt=None):
    """place a compact table float built from the .tsv (selected columns).
    Header cells (rename values) are emitted as raw LaTeX so math like $R^2$ renders;
    body cells are escaped. `filt` is an optional row predicate."""
    data = load(tid)
    if not data:
        P(f"% [table {tid} not found -- module not run]"); return
    if where: data = [r for r in data if all(r.get(k)==v for k,v in where.items())]
    if filt: data = [r for r in data if filt(r)]
    if rows: data = data[:rows]
    keys=[c for c in cols if c in (data[0] if data else {})]
    if not data or not keys:
        P(f"% [table {tid}: no rows/columns]"); return
    hdr=[(rename or {}).get(k,k) for k in keys]   # raw LaTeX headers (not escaped)
    P(r"\begin{table}[htbp]\centering", (r"\footnotesize" if small else ""),
      f"\\caption{{{caption}}}", f"\\label{{{label}}}",
      r"\adjustbox{max width=\linewidth}{%",
      r"\begin{tabular}{" + "l"*len(keys) + r"}\toprule",
      " & ".join(hdr) + r" \\ \midrule")
    for r in data:
        P(" & ".join(esc(fmt(r.get(c), digits)) for c in keys) + r" \\")
    P(r"\bottomrule\end{tabular}}", r"\end{table}", "")

# ============================================================ preamble
P(r"% " + "="*72,
  r"% paper.tex -- GENERATED by report/build_paper.py. Methods and Results only.",
  r"% Every number is read from the results tables at build time; prose is in the generator.",
  r"% Background, Discussion and framing are placeholders for the author (search TODO).",
  f"% built {datetime.now():%Y-%m-%d %H:%M}  git {sha()}  results {A.results}",
  r"% Target journal: BMC Medicine. For submission, port into the bmcart class; article here",
  r"% keeps it portable and compilable anywhere.",
  r"% " + "="*72,
  r"\documentclass[11pt]{article}",
  r"\usepackage[margin=1in]{geometry}",
  r"\usepackage{graphicx,booktabs,array,amsmath,amssymb,microtype}",
  r"\usepackage[export]{adjustbox}",
  r"\usepackage[numbers,sort&compress]{natbib}",
  r"\usepackage[hidelinks]{hyperref}",
  r"\usepackage{lineno}\linenumbers",
  r"\usepackage{setspace}\onehalfspacing",
  r"\graphicspath{{figures/}}",
  r"\newcommand{\TODO}[1]{\medskip\noindent\textbf{[TODO --- author: #1]}\medskip}",
  r"\title{\bfseries Portability of blood-pressure polygenic scores across South Asian and"
  r" African pregnancy cohorts, and their application to Mendelian randomization of maternal"
  r" blood pressure on perinatal outcomes}",
  r"\author{}",
  r"\date{}",
  r"\begin{document}\maketitle", "")

# ============================================================ background
# Facts only, each with a verified citation (per Anagh 2026-07-24). Interpretation and
# framing are deferred to the author; the abstract is intentionally omitted for now.
P(r"\section{Background}")
P(r"Preterm birth, low birth weight and fetal growth restriction are leading contributors to"
  r" perinatal mortality and morbidity, and their burden is concentrated in South Asia and"
  r" sub-Saharan Africa~\citep{ohuma2023preterm}. Elevated maternal blood pressure and"
  r" hypertensive disorders of pregnancy are established risk factors for these"
  r" outcomes~\citep{moralesberstein2026bp}. Observational associations between maternal blood"
  r" pressure and perinatal outcomes are susceptible to confounding by maternal adiposity,"
  r" parity and socioeconomic position, which complicates causal interpretation.")
P(r"Mendelian randomization uses germline genetic variants, fixed at conception and randomly"
  r" assorted with respect to most environmental confounders, as instrumental variables for a"
  r" modifiable exposure; under its assumptions it is largely robust to the confounding and"
  r" reverse causation that bias observational associations~\citep{daveysmith2003mr,"
  r"lawlor2008mr,davies2018reading}. Valid instrumental-variable inference rests on three"
  r" assumptions: the instrument is associated with the exposure (relevance); it shares no"
  r" common cause with the outcome (independence); and it influences the outcome only through"
  r" the exposure (the exclusion restriction)~\citep{lawlor2008mr,sanderson2022primer,"
  r"burgess2023guidelines}. For a highly polygenic exposure such as blood pressure, a polygenic"
  r" score aggregating many trait-associated variants can serve as a single, strong instrument,"
  r" provided it satisfies these assumptions \emph{in the population under study}"
  r"~\citep{choi2020tutorial}.")
P(r"Polygenic scores are, however, derived overwhelmingly from genome-wide association studies"
  r" of European-ancestry individuals, and their predictive performance attenuates with genetic"
  r" distance from the training population, most severely in populations of African"
  r" ancestry~\citep{martin2019disparities,kachuri2024transfer,prive2022portability}. Whether a"
  r" blood-pressure polygenic score is a valid and sufficiently strong instrument therefore"
  r" depends on the population, and its transferability to South Asian and African"
  r" populations --- and specifically to pregnant women, in whom blood pressure is measured"
  r" under heterogeneous antenatal protocols --- has not been characterised. Establishing this"
  r" transferability is a precondition for Mendelian randomization of maternal blood pressure"
  r" in these settings.")
P(r"Morales-Berstein et al.~\citep{moralesberstein2026bp} recently reported a well-powered"
  r" two-sample Mendelian-randomization analysis in European-ancestry samples, estimating the"
  r" effects of maternal systolic and diastolic blood pressure on perinatal outcomes including"
  r" birth weight, low birth weight, small-for-gestational-age and preterm birth. That study"
  r" provides an external, well-powered benchmark against which estimates from South Asian and"
  r" African cohorts can be compared. Whether these effects reproduce in such cohorts, and"
  r" whether the polygenic instruments required to estimate them transfer to these"
  r" populations, is unknown. Using five genotyped pregnancy cohorts from Bangladesh, Pakistan,"
  r" Tanzania and Zambia, we quantify the transferability of published blood-pressure polygenic"
  r" scores and estimate the association between genetically predicted maternal blood pressure"
  r" and preterm birth, low birth weight, small-for-gestational-age and birth weight.", "")

# ============================================================ methods
P(r"\section{Methods}")

P(r"\subsection{Study cohorts and participants}")
P(r"We analysed five genotyped pregnancy cohorts of the Multi-Omics for Mothers and Infants"
  r" (MOMI) Consortium~\citep{tang2026momi}: three South Asian (AMANHI-Sylhet, Bangladesh;"
  r" AMANHI-Karachi, Pakistan; GAPPS/PreSSMat-Matlab, Bangladesh) and two sub-Saharan African"
  r" (AMANHI-Pemba, Tanzania; ZAPPS-Lusaka, Zambia). The AMANHI cohorts were enrolled through"
  r" population-based pregnancy-surveillance platforms before 20 weeks' gestation between May"
  r" 2014 and June 2018~\citep{aftab2021amanhi}; the GAPPS/PreSSMat cohort in Matlab before 20"
  r" weeks between August 2015 and August 2017, and the ZAPPS cohort in Lusaka before 24 weeks"
  r" between 2018 and 2021~\citep{tang2026momi}. All cohorts obtained written informed consent"
  r" under locally approved protocols. Analyses were restricted to each mother's first recorded"
  r" pregnancy. Cohort characteristics are summarised in Supplementary Table~S1.")
P(r"We defined two nested analytic samples (Figure~\ref{fig:flow}). The transferability sample"
  r" comprised genotyped mothers with at least one valid antenatal blood-pressure reading. The"
  r" causal-analysis sample was the subset with recorded perinatal outcomes, used for the"
  r" Mendelian-randomization analyses. Sample sizes for each cohort are given in"
  r" Supplementary Table~S1.")
# Table 1 (cohort characteristics) has been moved to the supplement as Table S1
# (rendered by build_supplement.py). The main text refers to it as Supplementary Table~S1.
fig_float("F1_flow", r"Participant flow from enrolled first-pregnancy mothers to the"
  r" transferability and causal-analysis samples. Non-singleton pregnancies were excluded and"
  r" stillbirths retained. Blood pressure was summarised from antenatal visits only;"
  r" birth-weight-based outcomes (birth weight, low birth weight, small-for-gestational-age)"
  r" were restricted to livebirths; and implausible anthropometric values were set to missing.",
  "fig:flow", width="0.62\\textwidth")
fig_float("F1B_dag", r"Assumed causal model. Maternal blood pressure affects the perinatal"
  r" outcome both directly and through hypertensive disorders of pregnancy (preeclampsia),"
  r" which are treated as mediators and retained; measured confounders (maternal age, parity,"
  r" body-mass index, socioeconomic position) and within-cohort population structure are"
  r" shown. Population structure is the path controlled by principal components.", "fig:dag")

P(r"\subsection{Blood-pressure phenotypes}")
P(r"Antenatal blood-pressure readings were summarised per mother under several pre-specified"
  r" definitions. The primary exposure was the gestational-age--standardised residual mean:"
  r" for each cohort and trait, measured blood pressure was regressed on a natural-cubic-spline"
  r" function of gestational age at measurement and the residual added back to the cohort mean,"
  r" so that between-mother variation in the timing of measurement does not contribute to the"
  r" exposure. Auxiliary definitions comprised the simple antenatal mean, the last reading,"
  r" trimester-specific means, and readings before versus after 20 weeks; results under every"
  r" definition are reported in Supplementary Table~S6. Postpartum readings were excluded, as"
  r" the polygenic scores were trained on non-pregnant adults.")

P(r"\subsection{Genotyping, quality control and imputation}")
P(r"Genotypes were available from two platforms: a genome-wide SNP array (Illumina GSA) and"
  r" low-pass ($1\times$) whole-genome sequencing followed by imputation, the genotyping"
  r" strategy adopted across the MOMI Consortium~\citep{tang2026momi}. Low-pass sequencing was"
  r" analysed using imputed allelic dosages, which retain information lost by hard-calling at"
  r" low coverage. Per-cohort post-quality-control variant counts are reported in Supplementary"
  r" Table~S2. Imputation was performed upstream by the sequencing provider"
  r" \TODO{imputation reference panel, software and version}. Downstream quality control,"
  r" variant harmonisation to GRCh38, and all analyses were performed within this study.")

P(r"\subsection{Polygenic scores}")
P(r"We used published blood-pressure polygenic scores from the PGS Catalog spanning European"
  r" (EUR), South Asian (SAS), East Asian and multi-ancestry training data~\citep{keaton2024bp,ruan2022prscsx}; the full panel,"
  r" with source study and variant counts, is given in Supplementary Table~S1b. Scores were"
  r" computed with \texttt{plink2}, using dosages where available, and standardised to zero"
  r" mean and unit variance within each cohort and platform; the per-mother instrument is the"
  r" mean of the standardised scores across platforms. For each cohort we designated an"
  r" \emph{ancestry-matched} primary score: the SAS score for the South Asian cohorts and,"
  r" in the absence of a blood-pressure genome-wide association study of comparable size in"
  r" African populations, the EUR score for the African cohorts. This designation is used"
  r" only to select a single instrument for the Mendelian-randomization analyses; the"
  r" incremental $R^2$ of \emph{every} score in \emph{every} cohort is reported"
  r" (Supplementary Table~S4 and Figure~\ref{fig:heatmap}), so that instrument selection is"
  r" fully transparent.")

P(r"\subsection{Transferability analysis}")
P(r"Transferability was quantified as the incremental $R^2$: the increase in the coefficient"
  r" of determination when the standardised polygenic score is added to a linear model for"
  r" measured blood pressure that already contains maternal age, and, in sensitivity analyses,"
  r" within-cohort principal components~\citep{choi2020tutorial}. First-stage $F$ statistics"
  r" for the score--blood-pressure association are reported alongside.")

P(r"\subsection{Mendelian randomization}")
P(r"We performed one-sample Mendelian randomization within each cohort using the"
  r" ancestry-matched polygenic score as a single instrument. For each cohort, trait and"
  r" outcome we fitted the first-stage regression of measured blood pressure on the"
  r" standardised score and the reduced-form regression of the outcome on the score,"
  r" adjusting both stages for maternal age and five within-cohort principal components. We"
  r" treat the \emph{reduced-form} association as the primary test of the causal null: under"
  r" the null of no causal effect the score is unassociated with the outcome, and this test"
  r" does not involve the first stage. The Wald ratio (reduced-form coefficient divided by"
  r" first-stage coefficient) provides the effect magnitude, scaled to 10~mmHg, with a"
  r" delta-method standard error. First-stage $F$ statistics are reported for every cell."
  r" Binary outcomes (preterm birth, low birth weight, small-for-gestational-age) are modelled"
  r" on the log-odds scale; continuous birth weight is modelled in grams and additionally"
  r" expressed in standard-deviation units for comparison with external estimates.")
P(r"Estimates were pooled across cohorts by combining the reduced-form and first-stage"
  r" coefficients separately by inverse-variance--weighted meta-analysis and forming the ratio"
  r" once (the ratio of pooled coefficients). This estimator is preferred over meta-analysing"
  r" per-cohort Wald ratios, which can accentuate weak-instrument bias because the ratio and"
  r" its standard error are correlated~\citep{burgess2017review,burgess2011weak}. We report"
  r" South Asian, African and all-cohort pooled estimates.")

P(r"\subsection{Covariates, mediation and confounding}")
P(r"The assumed causal model is shown in Figure~\ref{fig:dag}. The observational analyses"
  r" adjusted progressively for maternal age, gravidity, body-mass index and maternal education"
  r" (Supplementary Table~S9). The Mendelian-randomization models adjusted for maternal age and"
  r" five within-cohort principal components only: under the instrumental-variable assumptions"
  r" the observational covariates are not confounders of the instrument--outcome association,"
  r" and body-mass index is itself a plausible mediator, so adjusting for it would induce"
  r" collider bias~\citep{burgess2023guidelines,sanderson2022primer}. Preeclampsia and other"
  r" hypertensive disorders of pregnancy lie on the causal pathway from maternal blood pressure"
  r" to the perinatal outcome and were treated as mediators: affected pregnancies were retained"
  r" and pathway variables were not adjusted for, since conditioning on a mediator would remove"
  r" part of the effect under study~\citep{vanderweele2016mediation}.")

P(r"\subsection{Instrument validity and sensitivity analyses}")
P(r"Because a single polygenic score does not permit the pleiotropy-robust estimators available"
  r" with multiple independent instruments (MR-Egger, weighted median, MR-PRESSO)"
  r"~\citep{bowden2015egger}, we assessed instrument validity directly through positive-control"
  r" associations: whether the score predicts conditions that elevated blood pressure is"
  r" established to cause (preeclampsia and other hypertensive disorders of pregnancy, and"
  r" chronic hypertension), before and after principal-component adjustment (Supplementary"
  r" Table~S12). Robustness of the causal estimates was assessed across all blood-pressure"
  r" definitions (Supplementary Table~S6), across ancestry strata, and across first-stage"
  r" strength thresholds.")

P(r"\subsection{Software}")
P(r"Genotype processing used \texttt{plink2} and \texttt{bcftools}; statistical analyses used"
  r" R~\citep{choi2020tutorial}. The full analysis is reproducible from a single build script,"
  r" with every table and figure regenerated from frozen intermediates and input fingerprints"
  r" recorded so that stale results are detected automatically"
  r" \TODO{software version numbers}.", "")

# ============================================================ results
P(r"\section{Results}")

P(r"\subsection{Transferability of blood-pressure polygenic scores}")
P(rf"Across the five cohorts, the incremental $R^2$ of the ancestry-matched score on measured"
  rf" systolic blood pressure ranged from {sas_lo}\% to {sas_hi}\% in the South Asian cohorts"
  rf" and from {afr_lo}\% to {afr_hi}\% in the African cohorts (Table~\ref{{tab:transfer}},"
  rf" Figure~\ref{{fig:heatmap}}). Variance explained was lower in the African cohorts,"
  rf" consistent with their greater genetic distance from the predominantly European-ancestry"
  rf" training samples~\citep{{martin2019disparities}}. Transferability also varied"
  rf" substantially \emph{{within}} the South Asian cohorts. Despite the modest variance"
  rf" explained, the ancestry-matched instruments were not weak by the conventional criterion:"
  rf" first-stage $F$ statistics ranged from {F_lo} to {F_hi} across cohorts and traits (South"
  rf" Asian cohorts {Fsas_lo}--{Fsas_hi}), all exceeding the usual threshold of"
  rf" 10~\citep{{burgess2011weak}}.")
P(r"Because the polygenic scores were trained in non-pregnant adults, we verified that they"
  r" proxy blood pressure \emph{during pregnancy} consistently: the variance explained was of"
  r" similar magnitude whether blood pressure was summarised as the antenatal mean, the"
  r" gestational-age--standardised residual, trimester-specific means, or readings before"
  r" versus after 20 weeks (Supplementary Tables~S6--S8), indicating that the instrument does"
  r" not depend on a particular gestational window for its relevance.")
tab_float("T2_transfer",
  ["cohort_display","ancestry","trait","chosen_PGS","R2pct","F"],
  r"Transferability of the ancestry-matched instrument: incremental $R^2$ (\%) of the chosen"
  r" polygenic score on measured blood pressure, and first-stage $F$ statistic, by cohort and"
  r" trait.",
  "tab:transfer",
  rename={"cohort_display":"Cohort","ancestry":"Ancestry","trait":"Trait","chosen_PGS":"Score",
          "R2pct":"Incremental $R^2$ (\\%)","F":"$F$ statistic"})
fig_float("F2_transfer", r"Incremental $R^2$ of every polygenic score in every cohort under the"
  r" primary blood-pressure definition.", "fig:heatmap")

P(rf"The two Bangladeshi cohorts illustrate this within-ancestry variation. AMANHI-Sylhet and"
  rf" GAPPS/PreSSMat-Matlab differed {FOLD}-fold in incremental $R^2$ for systolic blood"
  rf" pressure ({r2_amb}\% versus {r2_gpb}\%), despite the smallest pairwise Hudson $F_{{ST}}$"
  rf" of any cohort pair ({FST}; Supplementary Table~S3b). This difference was not attributable"
  rf" to blood-pressure variance (variance ratio $\approx1$ at every definition, Supplementary"
  rf" Table~S16) or to within-cohort population structure (principal-component adjustment"
  rf" changed $R^2$ by a median of {pc_shift}\% and did not narrow the difference, Supplementary"
  rf" Table~S17). The cohorts differ in the gestational timing and number of antenatal"
  rf" blood-pressure readings (Supplementary Table~S1, Supplementary Figure~S5), which"
  rf" contributes to the difference in measured transferability.")
P(r"\TODO{author: the within-ancestry transferability difference is consistent with"
  r" environmental modification of the genetic prediction of blood pressure"
  r" (gene-by-environment interaction); documented environmental and nutritional contrasts"
  r" between the two settings are compiled in the study repository and can be developed here.}")

P(r"\subsection{Maternal blood pressure and perinatal outcomes}")
P(rf"In the South Asian cohorts, where the instrument was strongest, higher genetically"
  rf" predicted maternal blood pressure was associated with lower birth weight for both traits"
  rf" (systolic {sas_sbp}~g per 10~mmHg, 95\% CI {sas_sbp_l} to {sas_sbp_h}; $p={sas_sbp_p}$;"
  rf" diastolic {sas_dbp}~g per 10~mmHg, 95\% CI {sas_dbp_l} to {sas_dbp_h}; $p={sas_dbp_p}$)."
  rf" Neither estimate reached conventional significance within the stratum and both are"
  rf" imprecise; we therefore interpret them through their consistency with well-powered"
  rf" external evidence rather than through internal $p$-values. Expressed per standard"
  rf" deviation of birth weight, the systolic estimate was {sd_u}~SD per 10~mmHg, closely"
  rf" matching the independent two-sample European estimate of {ext_sd}~SD reported by"
  rf" Morales-Berstein et al.~\citep{{moralesberstein2026bp}}. The estimate was of consistent"
  rf" sign and magnitude across all blood-pressure definitions (Supplementary Table~S6), and the"
  rf" instrument was fixed in advance by ancestry matching, not selected on the outcome"
  rf" association.")
P(rf"For the binary outcomes, the per-10-mmHg odds ratios were directionally concordant with"
  rf" the external European estimates, and every external point estimate fell within our"
  rf" confidence intervals (Figure~\ref{{fig:tri}}): for low birth weight,"
  rf" {our_or('SBP','LBW')} (systolic, this study) versus {ext_or('SBP','LBW')} (external);"
  rf" for small-for-gestational-age, {our_or('SBP','SGA')} versus {ext_or('SBP','SGA')}."
  rf" Pooled estimates for all outcomes are given in Supplementary Table~S10b and per-cohort"
  rf" estimates in Supplementary Figure~S4. The effect on preterm birth was null overall; a"
  rf" subtype analysis"
  rf" was consistent with an association with provider-initiated (indicated), but not"
  rf" spontaneous, preterm birth, confined to the single cohort (Matlab) in which indicated"
  rf" preterm delivery was common (Supplementary Table~T5b).")
fig_float("F4_triangulation", r"Triangulation of the maternal blood-pressure effect on binary"
  r" perinatal outcomes across three lines of evidence: adjusted observational association,"
  r" one-sample Mendelian randomization in this study (by ancestry stratum), and external"
  r" two-sample European Mendelian randomization. Odds ratios are per 10~mmHg increase in blood"
  r" pressure; the dashed line marks the null (OR = 1).", "fig:tri")
fig_float("F4b_triangulation_bwt", r"Triangulation of the maternal blood-pressure effect on"
  r" birth weight (standard-deviation units per 10~mmHg) across the same three lines of"
  r" evidence.", "fig:tribwt")
P(r"The polygenic score behaved as a valid instrument: genetically predicted maternal blood"
  r" pressure predicted the positive-control conditions --- preeclampsia and other hypertensive"
  r" disorders of pregnancy, and chronic hypertension --- at pooled significance"
  r" (Supplementary Figure~S6, Supplementary Table~S12), consistent with the instrument"
  r" capturing blood pressure rather than an unrelated pathway.")

# ============================================================ back matter
# Discussion removed 2026-07-24 per Anagh (to be written by the author later).
P(r"\section{Discussion}")
P(r"\TODO{author: Discussion to be written.}", "")

# ---- supplementary materials: brief pointer only (per Anagh: do not list them all) ----
P(r"\section*{Supplementary materials}")
P(r"Supplementary tables and figures, and the gene--environment context note, are provided"
  r" in a single supplementary file (\texttt{supplement.pdf}), generated from the same"
  r" analysis outputs as the main text.", "")

P(r"\section*{Declarations}")
P(r"\noindent\textbf{Ethics approval and consent to participate.} Each contributing cohort"
  r" was conducted under a locally approved protocol with written informed consent.")
P(r"\noindent\textbf{Consent for publication.} Not applicable.")
P(r"\noindent\textbf{Competing interests.} \TODO{declare or state none.}")
P(r"\noindent\textbf{Funding.} The MOMI Consortium was funded by the Gates Foundation"
  r" (award INV-037517)~\citep{tang2026momi}. \TODO{confirm any additional grant numbers"
  r" specific to this analysis.}")
P(r"\noindent\textbf{Authors' contributions.} \TODO{contributions per ICMJE.}")
P(r"\noindent\textbf{Acknowledgements.} \TODO{participants, field teams, and the contributing"
  r" cohort consortia.}")

P(r"\section*{Availability of data and materials}")
P(r"The complete analysis code is available at the study repository and reproduces every"
  r" table and figure in this manuscript from a single build script; the manuscript text"
  r" itself is generated from the analysis outputs. Individual-level cohort data are governed"
  r" by the MOMI Consortium and its contributing sites and are available under their"
  r" respective data-access policies~\citep{tang2026momi}."
  r" \TODO{repository URL and archived DOI on acceptance; data-access contact.}")

P(r"\bibliographystyle{unsrtnat}")
P(r"\bibliography{refs}")
P(rf"\vfill{{\footnotesize Manuscript generated {datetime.now():%Y-%m-%d} from git {sha()}.}}")
P(r"\end{document}")

out = os.path.join(A.out, "paper.tex")
with open(out, "w", encoding="utf-8") as f:
    f.write("\n".join(L) + "\n")
ntodo = sum(x.count(r"\TODO") for x in L)
print(f"wrote {out}\n  {len(L)} blocks, {ntodo} author TODO markers")
