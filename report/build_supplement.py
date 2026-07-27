#!/usr/bin/env python3
"""
build_supplement.py -- generate the supplementary file from the results directory.

Assembles the full supplementary tables and figures plus the gene-environment context note
(S18), as a single compilable document. No STROBE-MR checklist (removed 2026-07-23 per Anagh).
Full tables use longtable so many-row tables break cleanly across pages; only sensible column
subsets are shown to keep width on the page.

  python3 report/build_supplement.py --results DIR --out DIR
"""
import argparse, csv, math, os, shutil
from datetime import datetime

ap = argparse.ArgumentParser()
ap.add_argument("--results", default=os.environ.get("MOMI_RESULTS",""))
ap.add_argument("--out", required=True)
A = ap.parse_args()
TBL = os.path.join(A.results,"tables"); FIG = os.path.join(A.results,"figures")
OUTFIG = os.path.join(A.out,"figures"); os.makedirs(OUTFIG, exist_ok=True)

def esc(x):
    s = "" if x is None else str(x)
    if s.lower() in ("na","nan","none"): return ""
    s = s.replace("\\", r"\textbackslash{}")
    for c in "&%$#_{}": s = s.replace(c,"\\"+c)
    s = s.replace("<", r"\textless{}").replace(">", r"\textgreater{}")
    return s.replace("~", r"\textasciitilde{}")

def fmt(v,d=3):
    if v is None: return ""
    s=str(v).strip()
    if s=="" or s.upper() in ("NA","NAN","NULL"): return ""
    try: f=float(s)
    except ValueError: return s
    if math.isnan(f): return ""
    if f!=0 and (abs(f)<1e-3 or abs(f)>=1e6): return f"{f:.2e}"
    return f"{round(f,d):g}"

def load(tid):
    p=os.path.join(TBL,tid+".tsv")
    if not os.path.exists(p): return None
    return list(csv.DictReader(open(p,newline='',encoding='utf-8',errors='ignore'),delimiter='\t')) or None

L=[]; add=L.extend
def supp_table(sid, tid, cols, caption, rename=None, digits=3, where=None, landscape=False):
    """One supplementary table, its own page. Manual 'Table SX.' label (matches the main
    text), booktabs, small font, optional landscape for wide tables."""
    rows=load(tid)
    if not rows: L.append(f"% [{sid} {tid} absent]"); return
    if where: rows=[r for r in rows if all(r.get(k)==v for k,v in where.items())]
    keys=[c for c in cols if c in (rows[0] if rows else {})]
    if not rows or not keys: L.append(f"% [{sid} {tid} empty]"); return
    hdr=[(rename or {}).get(k,k) for k in keys]   # raw LaTeX headers (math renders)
    if landscape: add([r"\begin{landscape}"])
    add([r"{\small",
         r"\begin{longtable}{" + "l"*len(keys) + "}",
         rf"\caption*{{\textbf{{Table {sid}.}} {caption}}}\label{{tab:{sid.lower()}}}\\",
         r"\toprule", " & ".join(hdr)+r" \\ \midrule\endfirsthead",
         r"\toprule", " & ".join(hdr)+r" \\ \midrule\endhead", r"\bottomrule\endfoot"])
    for r in rows:
        add([" & ".join(esc(fmt(r.get(c),digits)) for c in keys)+r" \\"])
    add([r"\end{longtable}", r"}"])
    if landscape: add([r"\end{landscape}"])
    add([r"\clearpage", ""])            # one table per page

def supp_fig(sid, fid, caption, width=r"0.95\textwidth"):
    """One supplementary figure, its own page. Non-floating so it never drifts past its
    section heading; manual 'Figure SX.' label; no title baked into the plot."""
    src=None
    for e in ("pdf","png"):
        p=os.path.join(FIG,f"{fid}.{e}")
        if os.path.exists(p): src=p; break
    if not src: L.append(f"% [{sid} {fid} absent]"); return
    shutil.copy2(src, os.path.join(OUTFIG, os.path.basename(src)))
    add([r"\begin{center}",
         rf"\includegraphics[width={width},height=0.8\textheight,keepaspectratio]{{{fid}}}",
         r"\end{center}",
         r"\vspace{0.5em}",
         rf"\noindent\textbf{{Figure {sid}.}} {caption}\label{{fig:{sid.lower()}}}",
         r"\clearpage", ""])            # one figure per page

def table_cohorts(sid="S1"):
    """Transposed cohort-characteristics table (the former main-text Table 1): characteristics
    down the rows, cohorts (+ Overall) across the columns. Non-floating, supplement style."""
    data=load("T1_cohorts")
    if not data: L.append(f"% [{sid} T1_cohorts absent]"); return
    disp2short={
        "AMANHI (Sylhet, Bangladesh)":"Sylhet","AMANHI (Karachi, Pakistan)":"Karachi",
        "PreSSMat (Matlab, Bangladesh)":"Matlab","AMANHI (Pemba, Tanzania)":"Pemba",
        "ZAPPS (Lusaka, Zambia)":"Lusaka","Overall":"Overall"}
    sub={disp2short[r["cohort_display"]]:r for r in data if r.get("cohort_display") in disp2short}
    cols=[c for c in ["Sylhet","Karachi","Matlab","Pemba","Lusaka","Overall"] if c in sub]
    if not cols: L.append(f"% [{sid} T1_cohorts empty]"); return
    def cellval(short,key):
        v=sub[short].get(key,""); v="" if v is None else str(v).strip()
        if v in ("","NA","nan"): return "--"
        return v.replace("±",r"$\pm$").replace("%",r"\%").replace("_",r"\_")
    grp=[("\\textit{Sample sizes}",None),
         ("Genotyped, $n$","N_genotyped_union"),
         ("Transferability sample, $n$","N_partI"),
         ("Causal-analysis sample, $n$","N_partII"),
         ("\\textit{Maternal characteristics}",None),
         ("Age, years","age"),
         ("Gravidity","gravidity"),
         ("Body-mass index, kg/m$^2$","BMI"),
         ("Systolic blood pressure, mmHg","SBP"),
         ("Diastolic blood pressure, mmHg","DBP"),
         ("Gestational age at first BP reading, wk","GA_first_med_wk"),
         ("\\textit{Pregnancy and perinatal outcomes}",None),
         ("Gestational age at delivery, wk","GAdel_wk"),
         ("Birth weight, g","BWT_g"),
         ("Preterm birth, $n$ (\\%)","PTB_pct"),
         ("Low birth weight, $n$ (\\%)","LBW_pct"),
         ("Small-for-gestational-age, $n$ (\\%)","SGA_pct"),
         ("Stillbirth, $n$ (\\%)","still_pct"),
         ("Preeclampsia / HDP, $n$ (\\%)","PE_n_pct"),
         ("Chronic hypertension, $n$ (\\%)","cHTN_pct")]
    ncol=len(cols)
    add([r"{\small",
         rf"\noindent\textbf{{Table {sid}.}} Characteristics of the five genotyped cohorts"
         r" (continuous variables as mean~$\pm$~standard deviation; categorical as $n$ (\%))."
         r" Descriptives are computed on the transferability sample (genotyped mothers with a"
         r" valid antenatal blood-pressure reading); birth weight and gestational age at delivery"
         r" are among livebirths and stillbirth is over all pregnancies."
         rf"\label{{tab:{sid.lower()}}}",
         r"\vspace{0.6em}",
         r"\begin{center}",
         r"\adjustbox{max width=\linewidth}{%",
         r"\begin{tabular}{l"+"r"*ncol+r"}\toprule",
         "Characteristic & "+" & ".join(cols)+r" \\ \midrule"])
    for label,key in grp:
        if key is None:
            add([rf"\multicolumn{{{ncol+1}}}{{l}}{{{label}}} \\"])
        else:
            add([label+" & "+" & ".join(cellval(c,key) for c in cols)+r" \\"])
    add([r"\bottomrule\end{tabular}}", r"\end{center}", r"}", r"\clearpage", ""])

# ---------------------------------------------------------------- preamble
add([r"% supplement.tex -- GENERATED by report/build_supplement.py. Do not edit by hand.",
     f"% built {datetime.now():%Y-%m-%d %H:%M}  results {A.results}",
     r"\documentclass[11pt]{article}",
     r"\usepackage[margin=1in]{geometry}",
     r"\usepackage{graphicx,booktabs,longtable,array,amsmath,amssymb}",
     r"\usepackage[export]{adjustbox}",
     r"\usepackage{pdflscape}",
     r"\usepackage[numbers,sort&compress]{natbib}",
     r"\usepackage[hidelinks]{hyperref}",
     r"\graphicspath{{figures/}}",
     r"\title{Supplementary material\\ \large Portability of blood-pressure polygenic scores"
     r" and Mendelian randomization of maternal blood pressure on perinatal outcomes}",
     r"\author{}\date{}",
     r"\begin{document}\maketitle\tableofcontents\clearpage", ""])

# ---------------------------------------------------------------- supp tables
add([r"\section{Supplementary tables}"])

# Table S1: cohort characteristics (the detailed descriptive table, moved out of the main text).
table_cohorts("S1")

# Static outcome/exposure definitions table (Morales ST1 analogue). Hand-written, not from a
# results .tsv, so authored directly here.
add([r"{\small",
     r"\begin{longtable}{p{3.2cm} p{11.5cm}}",
     r"\caption*{\textbf{Table S1a.} Definitions of the exposure and perinatal outcomes.}\\",
     r"\toprule Variable & Definition \\ \midrule \endfirsthead",
     r"\toprule Variable & Definition \\ \midrule \endhead \bottomrule \endfoot",
     r"Maternal blood pressure (exposure) & Gestational-age--standardised residual mean of antenatal systolic (SBP) and diastolic (DBP) blood pressure, in mmHg; auxiliary definitions (simple mean, last reading, trimester means, $<$20 vs $\geq$20 weeks) are reported in Supplementary Table~S6. Postpartum readings excluded. \\",
     r"Preterm birth (PTB) & Live birth before 37 completed weeks of gestation. \\",
     r"Low birth weight (LBW) & Birth weight below 2500\,g (livebirths). \\",
     r"Small-for-gestational-age (SGA) & Birth weight below the 10th percentile for gestational age and sex (INTERGROWTH-21st standard). \\",
     r"Birth weight (BWT) & Continuous birth weight in grams (livebirths); also expressed in standard-deviation units. \\",
     r"Preterm-birth subtype & Spontaneous versus provider-initiated onset, where the MOMI labour-onset variable is available (GAPPS cohorts). \TODO{confirm subtype coding per site} \\",
     r"Hypertensive disorders of pregnancy & Preeclampsia and gestational hypertension; used as positive-control outcomes for the instrument. \\",
     r"\end{longtable}", r"}", r"\clearpage", ""])

supp_table("S1b","S1_prs_manifest",
  ["id","trait","anc","source","role","n_variants","method"],
  "Polygenic score panel: source study, training ancestry, variant count and construction method.",
  rename={"id":"PGS","trait":"Trait","anc":"Anc.","source":"Source","role":"Role","n_variants":"Variants","method":"Method"})
supp_table("S2","S2_variant_qc",
  ["cohort","platform","n_samples_final","variants_preQC_derived","removed_maf","removed_missing","removed_hwe","variants_postQC_file","pct_removed_total"],
  "Per-cohort, per-platform variant quality control: samples, variants before and after QC, and variants removed by filter.",
  landscape=True,
  rename={"cohort":"Cohort","platform":"Platform","n_samples_final":"Samples","variants_preQC_derived":"Pre-QC variants","removed_maf":"Rm.\\ MAF","removed_missing":"Rm.\\ miss.","removed_hwe":"Rm.\\ HWE","variants_postQC_file":"Post-QC variants","pct_removed_total":"\\% removed"})
supp_table("S3b","S3b_pairwise_fst",
  ["cohort_a","cohort_b","anc_a","anc_b","n_snps","fst_hudson"],
  "Pairwise Hudson $F_{ST}$ between cohorts at score SNPs.",
  rename={"cohort_a":"Cohort A","cohort_b":"Cohort B","anc_a":"Anc A","anc_b":"Anc B","n_snps":"SNPs","fst_hudson":"$F_{ST}$"}, digits=5)
supp_table("S4","S4_transferability",
  ["score_id","score_anc","cohort","trait","definition","R2pct","F","N"],
  "Full transferability grid: incremental $R^2$ of every score in every cohort, primary definition.",
  rename={"score_id":"Score","score_anc":"Anc.","cohort":"Cohort","trait":"Trait","definition":"Def.","R2pct":"$R^2$\\%","F":"$F$","N":"$N$"})
supp_table("S5","S5_platform",
  ["cohort_display","trait","chosen_PGS","R2pct_gsa","N_gsa","R2pct_dosage","N_dosage"],
  "Transferability by genotyping platform: SNP array versus low-pass sequencing dosages.",
  rename={"cohort_display":"Cohort","trait":"Trait","chosen_PGS":"Score","R2pct_gsa":"$R^2$ array","N_gsa":"$n$ array","R2pct_dosage":"$R^2$ dosage","N_dosage":"$n$ dosage"})
supp_table("S6","S6_definition_sensitivity",
  ["definition","trait","outcome","k","N","theta_per10","lo","hi","p"],
  "Causal estimates under every blood-pressure definition (South Asian stratum).",
  where={"stratum":"SAS"},
  rename={"definition":"Def.","trait":"Trait","outcome":"Outcome","k":"$k$","N":"$N$","theta_per10":"$\\theta$/10mmHg","lo":"Lo","hi":"Hi","p":"$p$"})
supp_table("S7","S7_residual",
  ["cohort_display","trait","chosen_PGS","R2pct_mean","R2pct_resid","N"],
  "Gestational-age--standardised (residual) versus mean blood pressure.",
  rename={"cohort_display":"Cohort","trait":"Trait","chosen_PGS":"Score","R2pct_mean":"$R^2$ mean","R2pct_resid":"$R^2$ resid","N":"$N$"})
supp_table("S8","S8_wk20",
  ["cohort_display","trait","R2pct_lt20","N_lt20","R2pct_ge20","N_ge20","note"],
  "Chronic ($<$20 wk) versus gestational ($\\geq$20 wk) blood pressure.",
  rename={"cohort_display":"Cohort","trait":"Trait","R2pct_lt20":"$R^2$ $<$20","N_lt20":"$n$ $<$20","R2pct_ge20":"$R^2$ $\\geq$20","N_ge20":"$n$ $\\geq$20","note":"Note"})
supp_table("S9","S9_confounding",
  ["trait","definition","outcome","scope","model","sample","OR","lo","hi","p"],
  "Observational associations under progressive confounder adjustment (pooled, own-sample).",
  where={"scope":"pooled","sample":"own"}, landscape=True,
  rename={"trait":"Trait","definition":"Def.","outcome":"Outcome","scope":"Scope","model":"Model","sample":"Sample","OR":"OR","lo":"Lo","hi":"Hi","p":"$p$"})
supp_table("S12","S12_controls",
  ["cohort","trait","variable","n","p_crude","p_pcadj"],
  "Positive-control associations: the blood-pressure polygenic score against conditions that elevated blood pressure is established to cause, crude and principal-component-adjusted.",
  where={"role":"positive"},
  rename={"cohort":"Cohort","trait":"Trait","variable":"Condition","n":"$n$","p_crude":"$p$ crude","p_pcadj":"$p$ PC-adj"})
supp_table("S10b","S10b_mr_pooled",
  ["trait","outcome","k_cohorts","theta_per10","se","p","I2"],
  "Pooled Mendelian-randomization estimates (all cohorts), per 10~mmHg, by the ratio of pooled coefficients.",
  rename={"trait":"Trait","outcome":"Outcome","k_cohorts":"Cohorts","theta_per10":"Effect per 10 mmHg","se":"SE","p":"$p$","I2":"$I^2$"})
supp_table("T4","T4_triangulation",
  ["trait","outcome","line","stratum","effect","eff_lo","eff_hi","p","units"],
  "Triangulation estimates by outcome, evidence line and ancestry stratum: adjusted observational, one-sample Mendelian randomization (this study), and external European two-sample Mendelian randomization, per 10~mmHg.",
  landscape=True,
  rename={"trait":"Trait","outcome":"Outcome","line":"Evidence","stratum":"Stratum","effect":"Effect","eff_lo":"Lo","eff_hi":"Hi","p":"$p$","units":"Units"})
supp_table("S16","S16_bpdist",
  ["trait","definition","mean_AMANHI_B","sd_AMANHI_B","mean_GAPPS_B","sd_GAPPS_B","SMD","var_ratio","KS_p"],
  "Blood-pressure distributions, AMANHI-Sylhet versus GAPPS-Matlab, by definition.",
  rename={"trait":"Trait","definition":"Def.","mean_AMANHI_B":"Mean Syl","sd_AMANHI_B":"SD Syl","mean_GAPPS_B":"Mean Mat","sd_GAPPS_B":"SD Mat","SMD":"SMD","var_ratio":"Var ratio","KS_p":"KS $p$"})
supp_table("S17","S17_pc_adjusted",
  ["cohort_display","trait","R2pct_age_only","R2pct_age_plus_PCs","pct_change"],
  "Transferability with and without within-cohort principal components.",
  rename={"cohort_display":"Cohort","trait":"Trait","R2pct_age_only":"$R^2$ age","R2pct_age_plus_PCs":"$R^2$ +PC","pct_change":"\\% change"})
supp_table("T5","T5_ptb_subtype",
  ["cohort","trait","subtype","N","n_case","F","rf_p","OR","OR_lo","OR_hi"],
  "Preterm birth by subtype (spontaneous versus provider-initiated), per cohort.",
  rename={"cohort":"Cohort","trait":"Trait","subtype":"Subtype","N":"$N$","n_case":"Cases","F":"$F$","rf_p":"RF $p$","OR":"OR","OR_lo":"Lo","OR_hi":"Hi"})

# ---------------------------------------------------------------- supp figures
add([r"\clearpage", r"\section{Supplementary figures}"])
supp_fig("S1","SF1_pca","Cohorts in genetic principal-component space, coloured by cohort. Axes are the first two within-sample principal components; each point is one participant (where per-individual data are available), otherwise the cohort mean with $\\pm$1 standard deviation.")
supp_fig("S2","SF2_distance_vs_r2","Genetic distance from the training population versus transferability.")
supp_fig("S3","SF3_bpdist","Blood-pressure distributions by cohort, earliest reading versus mean.")
supp_fig("S4","SF4_forest_mr","Per-cohort forest plots of the Mendelian-randomization estimates (binary outcomes).")
supp_fig("S4a","SF4a_forest_observational","Per-cohort observational associations, fully adjusted (age, gravidity, BMI, education).")
supp_fig("S4b","SF4b_forest_subtype","Per-cohort estimates by preterm-birth subtype (spontaneous versus provider-initiated).")
supp_fig("S4c","SF4c_forest_bwt","Per-cohort forest plot for birth weight.")
supp_fig("S5","F3_diagnostic","Blood-pressure measurement across the two Bangladeshi cohorts (AMANHI-Sylhet and GAPPS/PreSSMat-Matlab). (A) timing of blood-pressure measurement across gestation; (B) incremental $R^2$ by blood-pressure definition using a single polygenic score; (C) blood-pressure distributions.")
supp_fig("S6","poscontrol_forest","Positive-control associations: the blood-pressure polygenic score (systolic) against preeclampsia/hypertensive disorders of pregnancy and chronic hypertension, per cohort and pooled, principal-component-adjusted. Odds ratios are per standard deviation of the standardised score.")

# ---------------------------------------------------------------- S18 GxE note
add([r"\clearpage", r"\section{Supplementary Note S18: gene--environment context}"])
add([r"\small",
 r"The transferability of the blood-pressure polygenic score differs markedly between the two"
 r" Bangladeshi cohorts (AMANHI-Sylhet, $R^2\approx2.5\%$; GAPPS-Matlab, $R^2\approx6.5\%$)"
 r" despite negligible genetic differentiation between them (Hudson $F_{ST}\approx0.0002$)."
 r" The following documented environmental, nutritional and social contrasts between the two"
 r" regions support a gene-by-environment interpretation. \emph{The specific modifying"
 r" exposures are not identified by these data; the evidence below is contextual.}",
 "",
 r"\paragraph{Micronutrient status.} In early pregnancy in Matlab (MINIMat trial), 46\% of"
 r" women were vitamin~B\textsubscript{12} deficient and 18\% folate deficient"
 r" \citep{lindstrom2011minimat}. In the AMANHI-Sylhet cohort, folate deficiency in early"
 r" pregnancy was measured in a nested case-control study and associated with the risk of"
 r" preterm birth \citep{lazar2024folate}, establishing that folate status is both variable"
 r" and outcome-relevant in the Sylhet population. Assay and gestational timing differ"
 r" between the two studies, limiting direct quantitative comparison.",
 "",
 r"\paragraph{Drinking-water salinity and sodium.} Matlab lies in the tidal southeast, subject"
 r" to seasonal salinity intrusion; drinking-water salinity and urinary sodium in coastal"
 r" Bangladesh are associated with raised blood pressure and with (pre)eclampsia and"
 r" gestational hypertension \citep{khan2014salinity,scheelbeek2017salinity}. Sylhet is inland"
 r" freshwater (the northeastern haor basin).",
 "",
 r"\paragraph{Hypertension prevalence.} Adults in coastal Bangladesh have a higher prevalence"
 r" of hypertension than non-coastal residents (13.4\% versus 9.5\%; adjusted prevalence ratio"
 r" 1.29) \citep{hossain2025coastal}, consistent with a higher-blood-pressure environment in"
 r" the region containing Matlab.",
 "",
 r"\paragraph{Groundwater arsenic.} Chandpur district, which contains Matlab, lies in the"
 r" southeastern arsenic belt and is among the districts requiring arsenic mitigation; Sylhet's"
 r" dominant groundwater problem is iron and manganese rather than arsenic. District-level"
 r" figures should be taken from the BGS/DPHE (2001) national survey.",
 "",
 r"\paragraph{Surveillance and intervention.} Matlab is the site of a demographic-surveillance"
 r" system established in 1966, covering $\sim$230{,}000 residents, within which half the"
 r" population has been enrolled since 1977 in an intensive maternal--child-health and"
 r" family-planning programme with documented fertility and mortality effects"
 r" \citep{alam2017matlab}. Sylhet has no comparable multi-decade intervention-exposed"
 r" platform. Matlab is therefore an unusually well-characterised and health-serviced rural"
 r" population.",
 "",
 r"\emph{Verification and gaps.} Full source details and verification status for each item are"
 r" recorded in the study repository (\texttt{docs/supp\_gene\_environment.md}). Priorities for"
 r" confirmation before publication: a Sylhet-specific vitamin~D and B\textsubscript{12}"
 r" prevalence from AMANHI; exact district arsenic percentages; and any Matlab-specific"
 r" water-salinity or adult blood-pressure measurement.",
 "", r"\clearpage"])

# ---------------------------------------------------------------- Supplementary Methods
add([r"\section{Supplementary methods}", r"\small",
 r"\paragraph{Cohorts and data sources.} The five genotyped cohorts belong to the Multi-Omics"
 r" for Mothers and Infants (MOMI) Consortium~\citep{tang2026momi}: AMANHI-Sylhet,"
 r" AMANHI-Karachi and AMANHI-Pemba~\citep{aftab2021amanhi}, and the GAPPS cohorts"
 r" PreSSMat-Matlab and ZAPPS-Lusaka. Enrolment windows, genotyping platforms and sample sizes"
 r" are summarised in Table~\ref{tab:cohorts}. \TODO{per-cohort recruitment, consent and IRB"
 r" approval numbers.}",
 "",
 r"\paragraph{Genotyping and imputation.} Genotypes derived from a genome-wide SNP array"
 r" (Illumina GSA) and low-pass ($1\times$) whole-genome sequencing with imputation, the"
 r" MOMI genotyping strategy~\citep{tang2026momi}. \TODO{imputation reference panel, software"
 r" and version; post-imputation info-score filter.}",
 "",
 r"\paragraph{Polygenic scores.} Published blood-pressure polygenic scores were obtained from"
 r" the PGS Catalog spanning European, South Asian, East Asian and multi-ancestry training"
 r" data~\citep{keaton2024bp,ruan2022prscsx}; the panel is listed in Table~S1. Scores were"
 r" computed with \texttt{plink2} using allelic dosages where available and standardised within"
 r" each cohort and platform.",
 "",
 r"\paragraph{Principal components.} Within-cohort principal components used for adjustment were"
 r" computed per cohort from LD-pruned low-pass dosage genotypes (single platform, avoiding"
 r" cross-platform batch effects); five components were used, following convention"
 r"~\citep{burgess2023guidelines}.",
 "", r"\clearpage"])

# ---------------------------------------------------------------- STROBE-MR checklist
add([r"\section{STROBE-MR checklist}", r"\small",
 r"Reporting follows the STROBE-MR guideline for Mendelian-randomization studies"
 r"~\citep{strobemr2021}. The table maps each item to where it is addressed.", "",
 r"{\footnotesize\begin{longtable}{p{0.8cm} p{6.5cm} p{6.5cm}}",
 r"\toprule \# & Item & Addressed in \\ \midrule \endfirsthead",
 r"\toprule \# & Item & Addressed in \\ \midrule \endhead \bottomrule \endfoot"])
_strobe = [
 ("1","Title/abstract: state MR design","Title"),
 ("2","Background: rationale and assumptions","Background (relevance, independence, exclusion restriction)"),
 ("3","Objectives incl.\\ hypotheses","Background (final paragraph)"),
 ("4","Study design and framework","Methods: one-sample MR across five cohorts"),
 ("5","Setting, data sources","Methods: Study cohorts; Table~1"),
 ("6","Participants, sample flow","Figure~1 (CONSORT); Methods"),
 ("7","Exposure, outcomes, covariates, definitions","Methods; Table~S1a"),
 ("8","Genetic instrument and its selection","Methods: Polygenic scores; Table~S1"),
 ("9","Sources of bias","Methods: assumptions; positive controls; population structure"),
 ("10","Instrument strength (F, $R^2$)","Methods: Transferability; Table~2, Table~S4"),
 ("11","Statistical methods, pooling","Methods: Mendelian randomization"),
 ("12","Assessment of assumptions","Methods: Instrument validity and sensitivity analyses"),
 ("13","Sensitivity analyses","BP definitions (S6), ancestry strata, $F$-thresholds"),
 ("14","Descriptive data","Table~1"),
 ("15","Instrument--exposure associations","Table~2; Figure~2"),
 ("16","Main MR estimates","Figures~4--4b; Table~S10b"),
 ("17","Other analyses","Triangulation (Table~T4); subtype (T5)"),
 ("18","Key results","Results"),
 ("19","Limitations","Discussion \\TODO{on completion}"),
 ("20","Interpretation, generalisability","Discussion \\TODO{on completion}"),
]
for a,b,c in _strobe:
    add([f"{a} & {b} & {c} \\\\"])
add([r"\end{longtable}}", "",
 r"\clearpage", r"\bibliographystyle{unsrtnat}", r"\bibliography{refs}",
 r"\end{document}"])

out=os.path.join(A.out,"supplement.tex")
open(out,"w",encoding="utf-8").write("\n".join(L)+"\n")
ntab=sum(1 for x in L if x.startswith(r"\begin{longtable}"))
nfig=sum(1 for x in L if x.startswith(r"\begin{figure}"))
print(f"wrote {out}\n  {ntab} tables, {nfig} figures")
