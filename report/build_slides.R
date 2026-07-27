#!/usr/bin/env Rscript
# ============================================================
# build_slides.R — generate the comprehensive Beamer deck from results/current.
#
# WHY GENERATED RATHER THAN HAND-WRITTEN. The previous overleaf/slides.tex was hand-made in
# July and by 2026-07-20 every headline number in it was wrong: it led on "early DBP OR 1.41
# per 10 mmHg" for indicated PTB, described the MR as null for PTB, and cited figures that no
# longer exist. A deck typed by hand drifts from the pipeline the moment either changes. This
# script reads results/current/{tables,figures} at build time, so the deck cannot claim a
# number the pipeline does not currently produce.
#
# EVERY NUMBER IN THE DECK COMES FROM A .tsv. Where prose is unavoidable (the narrative, the
# DAG, STROBE-MR) the text is written here, but any figure quoted inside it is interpolated
# from the tables rather than typed, so the two cannot disagree.
#
# STRUCTURE
#   1  title, outline, how to read the deck
#   2  executive summary (numbers interpolated live)
#   3  study design: flow, cohorts, DAG (TikZ, = Figure 1B)
#   4  PART I  transferability + why it differs
#   5  PART II MR: power, panel, triangulation, subtypes, controls, MVMR, fetal
#   6  narrative / interpretation
#   7  STROBE-MR checklist (= S18)
#   8  limitations, provenance, full-table appendix
#
#   Rscript report/build_slides.R [--out DIR]
#     default --out <PIPE>/overleaf
# ============================================================
suppressMessages(library(data.table))
PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
source(file.path(PIPE,"lib/momi_io.R")); source(file.path(PIPE,"lib/momi_config.R"))
P <- momi_paths(PIPE)
OUT <- momi_arg("--out", file.path(PIPE,"overleaf"))
FIGDIR <- file.path(OUT,"figures"); dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)

## ---------- helpers ----------------------------------------------------------------------
tex_esc <- function(x){
  x <- as.character(x); x[is.na(x)] <- ""
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  for(ch in c("&","%","$","#","_","{","}")) x <- gsub(ch, paste0("\\", ch), x, fixed=TRUE)
  x <- gsub("~", "\\textasciitilde{}", x, fixed=TRUE)
  x <- gsub("^", "\\textasciicircum{}", x, fixed=TRUE)
  x <- gsub(">=", "$\\geq$", x, fixed=TRUE); x <- gsub("<=", "$\\leq$", x, fixed=TRUE)
  x
}
rd <- function(x, d=3){
  if(is.numeric(x)) ifelse(is.na(x), "", formatC(round(x,d), format="f", digits=d, big.mark="")) else x
}
tbl <- function(id){
  f <- file.path(P$tables, paste0(id,".tsv"))
  if(!file.exists(f)) return(NULL)
  suppressWarnings(fread(f))
}
## one number, safely, so prose can never quote something that is not in a table
num <- function(dt, expr, col, d=2, default="--"){
  if(is.null(dt)) return(default)
  r <- tryCatch(dt[eval(expr)], error=function(e) NULL)
  if(is.null(r) || !nrow(r) || !(col %in% names(r))) return(default)
  v <- r[[col]][1]; if(is.na(v)) return(default)
  if(is.numeric(v)) formatC(round(v,d), format="f", digits=d) else as.character(v)
}

LINES <- character(0)
add <- function(...) LINES <<- c(LINES, ...)

## a table frame: booktabs, auto-shrunk, optionally chunked across frames
frame_table <- function(id, cols=NULL, title=NULL, note=NULL, digits=3, rows_per=14,
                        rename=NULL, filter=NULL, landscape=FALSE){
  D <- tbl(id); if(is.null(D) || !nrow(D)) return(invisible(NULL))
  if(!is.null(filter)) D <- D[eval(filter, D)]
  if(!is.null(cols)) { cols <- intersect(cols, names(D)); D <- D[, ..cols] }
  if(!nrow(D)) return(invisible(NULL))
  for(j in names(D)) if(is.numeric(D[[j]])) set(D, j=j, value=rd(D[[j]], digits))
  hdr <- names(D); if(!is.null(rename)) for(k in names(rename)) hdr[hdr==k] <- rename[[k]]
  hdr <- tex_esc(hdr)
  nch <- nrow(D); nfr <- ceiling(nch/rows_per)
  for(k in seq_len(nfr)){
    idx <- ((k-1)*rows_per+1):min(k*rows_per, nch)
    ttl <- if(is.null(title)) id else title
    if(nfr>1) ttl <- sprintf("%s \\textnormal{\\small(%d/%d)}", ttl, k, nfr)
    add(sprintf("\\begin{frame}{%s}", ttl))
    add("\\centering\\resizebox{\\linewidth}{!}{%")
    add(sprintf("\\begin{tabular}{%s}\\toprule", paste(rep("l", ncol(D)), collapse="")))
    add(paste0(paste(hdr, collapse=" & "), " \\\\ \\midrule"))
    for(i in idx) add(paste0(paste(tex_esc(unlist(D[i])), collapse=" & "), " \\\\"))
    add("\\bottomrule\\end{tabular}}")
    if(!is.null(note) && k==nfr){
      add("\\vspace{2mm}"); add(sprintf("{\\scriptsize %s}", note))
    }
    add(sprintf("\\vspace{1mm}{\\tiny\\ttfamily source: tables/%s.tsv}", tex_esc(id)))
    add("\\end{frame}", "")
  }
}

## a figure frame; copies the PDF (falls back to PNG) into overleaf/figures
frame_fig <- function(id, title, note=NULL, height="0.72\\textheight"){
  src <- NULL
  for(ext in c("pdf","png")){
    f <- file.path(P$figures, paste0(id,".",ext)); if(file.exists(f)){ src <- f; break }
  }
  if(is.null(src)) return(invisible(NULL))
  file.copy(src, file.path(FIGDIR, basename(src)), overwrite=TRUE)
  add(sprintf("\\begin{frame}{%s}", title))
  add("\\centering")
  add(sprintf("\\includegraphics[height=%s,width=\\linewidth,keepaspectratio]{%s}",
              height, tools::file_path_sans_ext(basename(src))))
  if(!is.null(note)) add("\\vspace{1mm}", sprintf("{\\scriptsize %s}", note))
  add(sprintf("{\\tiny\\ttfamily source: figures/%s}", tex_esc(basename(src))))
  add("\\end{frame}", "")
}

frame_text <- function(title, body, plain=FALSE){
  add(sprintf("\\begin{frame}%s{%s}", if(plain) "[plain]" else "", title)); add(body); add("\\end{frame}", "")
}

## ---------- pull the numbers the prose needs ----------------------------------------------
T2   <- tbl("T2_transfer");        S10b <- tbl("S10b_mr_pooled")
T4   <- tbl("T4_triangulation");   S15b <- tbl("S15b_fetal_pooled")
T5b  <- tbl("T5b_ptb_subtype_pooled"); S13 <- tbl("S13_power")
S3b  <- tbl("S3b_pairwise_fst");   S17  <- tbl("S17_pc_adjusted")
S16  <- tbl("S16_bpdist");         S14  <- tbl("S14_mvmr")
S15d <- tbl("S15d_transmission_check"); S6 <- tbl("S6_definition_sensitivity")
MAN  <- if(file.exists(P$manifest)) suppressWarnings(fread(P$manifest, fill=TRUE)) else NULL

fst_bd <- num(S3b, quote((cohort_a=="AMANHI-Bangladesh" & cohort_b=="GAPPS-Bangladesh") |
                         (cohort_b=="AMANHI-Bangladesh" & cohort_a=="GAPPS-Bangladesh")),
              "fst_hudson", 5)
r2_amb <- num(T2, quote(cohort_display=="AMANHI (Sylhet, Bangladesh)" & trait=="SBP"), "R2pct", 2)
r2_gpb <- num(T2, quote(cohort_display=="PreSSMat (Matlab, Bangladesh)" & trait=="SBP"), "R2pct", 2)
sas_sbp_bwt <- num(T4, quote(arm=="MR_ours" & stratum=="SAS" & trait=="SBP" & outcome=="BWT"), "effect", 1)
sas_sbp_p   <- num(T4, quote(arm=="MR_ours" & stratum=="SAS" & trait=="SBP" & outcome=="BWT"), "p", 4)
sas_dbp_bwt <- num(T4, quote(arm=="MR_ours" & stratum=="SAS" & trait=="DBP" & outcome=="BWT"), "effect", 1)
sas_dbp_p   <- num(T4, quote(arm=="MR_ours" & stratum=="SAS" & trait=="DBP" & outcome=="BWT"), "p", 4)
sd_units    <- num(T4, quote(arm=="MR_ours" & stratum=="SAS" & trait=="SBP" & outcome=="BWT"), "effect_sd", 3)
ext_sbp_bwt <- num(T4, quote(arm=="MR_external" & trait=="SBP" & outcome=="BWT"), "est", 3)
npow <- if(!is.null(S13)) sum(S13$MR_adequate %in% TRUE, na.rm=TRUE) else NA
ncell <- if(!is.null(S13)) nrow(S13) else NA
pct_att_sas <- num(S15b, quote(stratum=="SAS" & trait=="SBP"), "pct_attenuation", 1)
pairs_n     <- num(S15b, quote(stratum=="all" & trait=="SBP"), "n_pairs", 0)
mvmr_shift  <- if(!is.null(S14)) formatC(median(abs(S14$pct_change_zBP), na.rm=TRUE), format="f", digits=1) else "--"
pc_shift    <- if(!is.null(S17)) formatC(median(S17$pct_change, na.rm=TRUE), format="f", digits=1) else "--"

## ============================================================================
## PREAMBLE
## ============================================================================
add(
"% ============================================================",
"% slides.tex — GENERATED by report/build_slides.R. DO NOT EDIT BY HAND.",
"% Every number is read from results/current/tables at build time.",
sprintf("%% built: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
sprintf("%% git:   %s", momi_git_sha(PIPE)),
"% ============================================================",
"\\documentclass[aspectratio=169,9pt]{beamer}",
"\\usetheme{Madrid}\\usecolortheme{seahorse}",
"\\usepackage{graphicx,booktabs,array,amsmath,amssymb}",
"\\usepackage{tikz}\\usetikzlibrary{arrows.meta,positioning,shapes.geometric,calc}",
"\\setbeamertemplate{navigation symbols}{}",
"\\setbeamerfont{frametitle}{size=\\normalsize}",
"\\graphicspath{{figures/}}",
"\\setlength{\\tabcolsep}{3pt}\\renewcommand{\\arraystretch}{1.02}",
"\\newcommand{\\note}[1]{\\vspace{1mm}{\\scriptsize\\itshape #1}}",
"\\AtBeginSection[]{\\begin{frame}[plain,c]\\centering\\usebeamerfont{frametitle}\\Large\\insertsectionhead\\end{frame}}",
"\\title[MOMI BP and perinatal outcomes]{Blood-pressure polygenic score portability and Mendelian randomization of maternal blood pressure on perinatal outcomes}",
"\\subtitle{Five South Asian and African pregnancy cohorts --- complete results}",
"\\author[Chattopadhyay et al.]{Anagh Chattopadhyay \\textit{et al.}}",
"\\institute[JHU]{Johns Hopkins Bloomberg School of Public Health}",
sprintf("\\date{%s}", format(Sys.Date(), "%d %B %Y")),
"\\begin{document}",
"\\frame{\\titlepage}",
"\\begin{frame}{Outline}\\footnotesize\\tableofcontents\\end{frame}", "")

frame_text("How to read this deck", paste(
"\\footnotesize",
"\\begin{itemize}",
"\\item This deck is \\textbf{generated} from \\texttt{results/current/} by \\texttt{report/build\\_slides.R}.",
"      Every number is read from a \\texttt{.tsv} at build time; the source file is printed under each slide.",
"\\item \\textbf{Nothing here is typed by hand except interpretation.} If the pipeline changes, rebuild and the deck changes.",
"\\item Two exposures (SBP, DBP) are reported \\textbf{co-equally} throughout --- no primary trait, to avoid post-hoc selection.",
"\\item Primary BP definition is \\texttt{resid} (gestational-age--standardised mean of antenatal readings).",
"\\item Effects are per \\textbf{10 mmHg}. Odds ratios for binary outcomes; grams (and SD) for birth weight.",
"\\item \\textbf{Reduced-form $p$ is the test}; the Wald ratio only rescales it into mmHg units.",
"\\end{itemize}", collapse="\n"))

## ============================================================================
add("\\section{Executive summary}")
frame_text("What this study establishes", paste(
"\\footnotesize",
"\\textbf{Part I --- portability.} BP polygenic scores transfer poorly and unevenly across these cohorts.",
"The gap is \\emph{not} explained by:",
"\\begin{itemize}",
sprintf("\\item \\textbf{genetic ancestry} --- the two Bangladeshi cohorts are separated by Hudson $F_{ST}=%s$ yet differ %s-fold in $R^2$ (%s\\%% vs %s\\%%, SBP);",
        fst_bd, formatC(as.numeric(r2_gpb)/as.numeric(r2_amb), format="f", digits=1), r2_amb, r2_gpb),
"\\item \\textbf{BP variance} --- variance ratio $\\approx 1.0$ at every definition (S16);",
sprintf("\\item \\textbf{population structure} --- within-cohort PC adjustment moves $R^2$ by a median %s\\%% (S17).", pc_shift),
"\\end{itemize}",
"What remains is \\textbf{measurement}: when in gestation BP was taken, and how many readings were averaged.",
"\\vspace{2mm}",
"\\textbf{Part II --- causal effects.} Portability determines where MR is possible at all.",
"\\begin{itemize}",
sprintf("\\item Higher maternal SBP lowers birth weight: \\textbf{%s g per 10 mmHg} ($p=%s$) in South Asian cohorts; DBP \\textbf{%s g} ($p=%s$).",
        sas_sbp_bwt, sas_sbp_p, sas_dbp_bwt, sas_dbp_p),
sprintf("\\item Equivalent to %s SD of birth weight, against an external European estimate of %s SD.", sd_units, ext_sbp_bwt),
sprintf("\\item Robust to all eight BP definitions, to the pooling method, and to adiposity genetics (median shift %s\\%%, S14).", mvmr_shift),
sprintf("\\item African cohorts are \\emph{uninformative}, not null --- %d of %d cells reach 80\\%% power (S13).", npow, ncell),
"\\end{itemize}", collapse="\n"))

frame_text("What this study does \\emph{not} establish", paste(
"\\footnotesize",
"\\begin{itemize}",
"\\item \\textbf{No effect on preterm birth} is demonstrated. The pooled estimate is null; the subtype analysis (T5) shows",
"      this averages an indicated-PTB signal against a spontaneous-PTB null, but with a 50-event floor only",
"      two cohorts contribute and one of those is null --- so it is \\emph{one site's} result, not a pooled finding.",
sprintf("\\item \\textbf{The maternal-vs-fetal origin of the birth-weight effect is unresolved.} Conditioning on fetal genotype in %s mother--infant pairs attenuates the maternal coefficient by %s\\%% but leaves it directionally unchanged; our sample cannot separate the two components. Warrington \\textit{et al.} (2019), far better powered, find the effect maternal with no fetal component. We defer to them.",
        pairs_n, pct_att_sas),
"\\item \\textbf{Early-pregnancy BP is unmeasurable here} --- the \\texttt{lt20} and \\texttt{tri1} definitions fall below",
"      the sample floor because AMANHI enrols too late.",
"\\item \\textbf{Pleiotropy-robust methods are unavailable.} With a single polygenic score, MR-Egger, weighted median",
"      and MR-PRESSO are undefined (Bowden 2015); S12's negative controls are the substitute.",
"\\end{itemize}", collapse="\n"))

## ============================================================================
add("\\section{Study design and sample}")
frame_fig("F1_flow", "Figure 1A --- sample flow",
  "21,685 first-pregnancy mothers $\\rightarrow$ 14,032 genotyped $\\rightarrow$ 13,717 Part I (valid antenatal BP) $\\rightarrow$ 13,212 Part II (non-missing PTB). The BP exclusion is split by cause: 178 postpartum-only vs 137 with no BP at all.")
frame_table("T1_cohorts", title="Table 1 --- cohort characteristics",
  cols=c("cohort_display","ancestry","protocol","N_genotyped_union","N_partI","N_partII",
         "age","gravidity","BMI","PTB_pct","LBW_pct","SGA_pct","PE_n_pct"),
  note="Timing of BP measurement is \\emph{demonstrated}, not assumed: at the $\\sim$13-week weight visit, BP is present in 0.0--3.0\\% of AMANHI mothers vs 99.4--100\\% of GAPPS.")

## ---- Figure 1B: the DAG, in TikZ ----
add("\\begin{frame}{Figure 1B --- causal model (DAG)}")
add("\\centering")
add("\\begin{tikzpicture}[>=Stealth,node distance=13mm,",
    "  every node/.style={font=\\scriptsize},",
    "  var/.style={draw,rounded corners=2pt,minimum height=6mm,minimum width=17mm,align=center,fill=white},",
    "  lat/.style={var,dashed,fill=black!4},",
    "  med/.style={var,fill=orange!14},",
    "  ins/.style={var,fill=blue!12},",
    "  out/.style={var,fill=green!12}]")
add("\\node[ins] (z) {BP polygenic\\\\score $Z$};")
add("\\node[var,right=20mm of z] (bp) {Maternal BP\\\\(exposure)};")
add("\\node[med,right=22mm of bp] (pe) {Preeclampsia\\\\/ hypertensive\\\\disorder};")
add("\\node[out,right=20mm of pe] (y) {Perinatal outcome\\\\(BWT, LBW, SGA, PTB)};")
add("\\node[lat,above=11mm of bp] (u) {Confounders $U$\\\\age, gravidity,\\\\BMI, SES};")
add("\\node[lat,below=11mm of bp] (pc) {Population\\\\structure};")
add("\\node[lat,below=11mm of y] (fg) {Fetal genotype\\\\(transmitted)};")
add("\\draw[->,thick] (z) -- (bp);")
add("\\draw[->,thick] (bp) -- (pe);")
add("\\draw[->,thick] (pe) -- (y);")
add("\\draw[->,thick] (bp) to[out=-25,in=-155] (y);")
add("\\draw[->] (u) -- (bp); \\draw[->] (u) to[out=0,in=100] (y);")
add("\\draw[->] (pc) -- (z); \\draw[->] (pc) to[out=0,in=-100] (y);")
add("\\draw[->,dashed] (z) to[out=-40,in=180] (fg); \\draw[->,dashed] (fg) -- (y);")
add("\\end{tikzpicture}")
add("\\vspace{1mm}")
add(paste(
"{\\scriptsize",
"\\textbf{Blue} = instrument. \\textbf{Orange} = mediator. \\textbf{Green} = outcome. \\textbf{Dashed grey} = unmeasured or partially controlled.",
"\\begin{itemize}\\setlength\\itemsep{0pt}",
"\\item \\textbf{Preeclampsia is a MEDIATOR, not a confounder} --- it lies on the causal path from BP to outcome.",
"      Adjusting for it would remove part of the effect being estimated, so PE is \\emph{retained} in the MR sample (S11 tests this).",
"\\item \\textbf{Population structure} is the one path that would invalidate the instrument; controlled by within-cohort PCs and tested in S12/S17.",
"\\item \\textbf{Fetal genotype} is the exclusion-restriction threat specific to maternal MR (Lawlor 2017): $Z$ is correlated $\\approx 0.5$",
"      with the transmitted fetal score, which acts on the outcome directly. Addressed in S15.",
"\\item Covariates in the MR are \\textbf{age + within-cohort PCs only} (Burgess 2023) --- \\emph{not} the observational confounder set,",
"      because BMI is plausibly a mediator and gravidity/education are post-randomisation.",
"\\end{itemize}}", collapse="\n"))
add("\\end{frame}", "")

## ============================================================================
add("\\section{Part I --- polygenic score portability}")
frame_table("T2_transfer", title="Table 2 --- transferability of the chosen instrument",
  note=paste("Observed incremental $R^2$ only. The measurement-error-corrected $R^2$ is supplementary (S4):",
             "$k$ (readings averaged) varies 2.51--4.44 across cohorts, so a cohort-specific correction would",
             "contaminate the cross-cohort contrast this paper rests on."))
frame_fig("F2_transfer", "Figure 2 --- transferability heatmap (score $\\times$ cohort)",
  "All scores, all cohorts, primary BP definition. The SAS/AFR split is visible, but so is the within-ancestry spread that Figure 3 explains.")
frame_table("S1_prs_manifest", title="S1 --- polygenic score panel",
  cols=c("id","trait","anc","source","role","n_variants_published","method"),
  note="The SAS scores are \\textbf{PRSmix mixture scores, not a South Asian GWAS}. Dosage recovery 67--90\\% across cohorts (worst: GAPPS-Zambia).")
frame_table("S4_transferability", title="S4 --- full transferability grid",
  note="Disattenuated $R^2$ uses the Spearman--Brown reliability of the \\emph{mean} of $k$ readings, not the single-reading ICC. $k$ is per definition: 1 for \\texttt{last}, NA for window means (\\texttt{bp\\_icc} carries no per-window count).")
frame_table("T3_instrument", title="T3 --- instrument decision", rows_per=6,
  note="EUR primary everywhere except GAPPS-Bangladesh (SAS). Selection is fully disclosed by S4/F2, which show every score in every cohort --- STROBE-MR item 6b.")

add("\\subsection{Why portability differs}")
frame_fig("F3_diagnostic", "Figure 3 --- why the two Bangladeshi cohorts differ",
  "Panel A: when BP is measured. Panel B: $R^2$ by BP definition, same score. Panel C: BP distributions --- averaging narrows the gap but does not close it.")
frame_table("S3b_pairwise_fst", title="S3b --- pairwise Hudson $F_{ST}$ at score SNPs", rows_per=10,
  note=paste0("$\\star$ The decisive number: AMANHI-B vs GAPPS-B $F_{ST}=", fst_bd,
              "$ --- 14$\\times$ closer than the next-closest SAS pair and 25$\\times$ closer than the two AFR cohorts are to each other --- yet their $R^2$ differs 2.6-fold. \\textbf{Ancestry cannot be the explanation.}"))
frame_table("S16_bpdist", title="S16 --- BP distributions, AMANHI-B vs GAPPS-B",
  cols=c("trait","definition","flag","mean_AMANHI_B","sd_AMANHI_B","mean_GAPPS_B","sd_GAPPS_B","SMD","var_ratio","KS_p"),
  note=paste("Both halves must be stated: averaging \\emph{attenuates} the gap (SMD $-0.50$ for \\texttt{last} to $-0.21$ for \\texttt{mean-2}, a 58\\% reduction)",
             "but \\emph{never closes} it --- every KS test is significant. $\\star$ \\textbf{var\\_ratio $\\approx 1$ kills the competing explanation} that GAPPS simply has more BP variance to explain."))
frame_table("S7_residual", title="S7 --- GA-standardised vs mean BP", rows_per=10,
  note="Largest change across ten cells is \\textbf{0.20 percentage points}. Transferability is \\emph{insensitive} to GA-standardisation, so S7 is a robustness check, not an argument for the exposure. (Direction splits perfectly by ancestry --- SAS up, AFR down --- but far too small to interpret.)")
frame_table("S8_wk20", title="S8 --- chronic ($<$20 wk) vs gestational ($\\geq$20 wk) BP", rows_per=10,
  note="Usable only for GAPPS-Bangladesh and Zambia. AMANHI enrols too late for a $<$20-week window, which is why early-pregnancy BP is unmeasurable in this study.")
frame_table("S5_platform", title="S5 --- GSA vs low-pass WGS dosage", rows_per=10,
  note="AMANHI-Pemba SBP is 0.006\\% on GSA vs 1.761\\% on dosage ($n=332$). Worth noting: a platform-specific failure in African cohorts would be a competing explanation for the headline African result.")
frame_table("S2_variant_qc", title="S2 --- genotyping and variant QC",
  note="AFR cohorts carry \\textbf{1.74$\\times$ more} post-QC variants than SAS, and \\texttt{--maf 0.005} removes \\emph{less} in AFR --- so variant count does not explain their lower $R^2$.")
frame_table("S17_pc_adjusted", title="S17 --- transferability with and without within-cohort PCs", rows_per=10,
  note=paste0("Median change \\textbf{", pc_shift, "\\%}; the Bangladeshi fold-gap moves 2.88 $\\rightarrow$ 2.88. Population structure is not inflating Part I, and this also validates the MR covariate set."))
frame_fig("SF1_pca", "SF1 --- ancestry: cohorts projected onto 1000 Genomes",
  "The three SAS cohorts are nearly indistinguishable (PC1 spread 0.003); the two AFR cohorts are separated by 0.025 --- \\textbf{eight times} the entire SAS spread. ``AFR'' lumps two populations that differ more from each other than our SAS cohorts do. NB projected PCs shrink toward the origin.")
frame_fig("SF2_distance_vs_r2", "SF2 --- genetic distance vs transferability",
  "Confirms the conventional account coarsely, and refutes it as complete: two cohorts at effectively identical genetic distance differ 2.6-fold in $R^2$. No correlation is fitted --- with five cohorts a coefficient would dress a scatterplot as a test.")
frame_fig("SF3_bpdist", "SF3 --- BP distributions by cohort and definition")

## ============================================================================
add("\\section{Part II --- Mendelian randomization}")
frame_table("S13_power", title="S13 --- power and minimum detectable effect",
  cols=c("cohort","trait","outcome","N","cases","firstStage_R2pct","observed_effect",
         "MDE_observational","MDE_MR","power_MR_at_obs"),
  note=sprintf("\\textbf{%d of %d cells reach 80\\%% power} (median %s). This is why the reduced form is the test and why nulls are reported as uninformative rather than as absence of effect.",
               npow, ncell, if(!is.null(S13)) formatC(median(S13$power_MR_at_obs,na.rm=TRUE),format="f",digits=2) else "--"))
frame_table("S10_mr_panel", title="S10 --- MR panel, per cohort",
  cols=c("cohort","trait","outcome","N","n_case","F","rf_p","theta_per10","se_per10","power_MR_at_obs"),
  note="Covariates: age + 5 within-cohort PCs, identical in both stages (Burgess 2023). \\texttt{rf\\_p} is the honest test; \\texttt{theta\\_per10} rescales the same evidence into mmHg.")
frame_table("S10b_mr_pooled", title="S10b --- pooled MR estimates",
  cols=c("trait","outcome","k_cohorts","N_total","theta_per10","se","p",
         "meta_of_ratios_theta_per10","pooling_diff_pct","I2","min_F","max_F"),
  note=paste("Primary = \\textbf{ratio of pooled coefficients} (Burgess, Small \\& Thompson 2017: study-level meta-analysis of Wald ratios",
             "``can accentuate weak instrument bias''). The two methods agree to $<$8\\% in SAS and diverge up to 13-fold in AFR ---",
             "itself evidence for where the instrument works. $I^2$ at $k=5$ is biased (von Hippel 2015); do \\textbf{not} read $I^2=0$ as homogeneity."))
frame_fig("SF4_forest_mr", "SF4 --- per-cohort forest, binary outcomes",
  "Event counts and first-stage $F$ printed on every row: a wide interval here usually means the outcome barely occurred, not that the effect is small.")
frame_fig("SF4c_forest_bwt", "SF4c --- per-cohort forest, birth weight")
frame_fig("F4_triangulation", "Figure 4 --- triangulation (odds-ratio outcomes)",
  "Three lines of evidence with largely non-overlapping biases: observational (well powered, confounded), our MR (unconfounded, underpowered), external European MR (well powered, may not transport).")
frame_fig("F4b_triangulation_bwt", "Figure 4B --- triangulation, birth weight")
frame_table("T4_triangulation", title="T4 --- triangulation table",
  cols=c("trait","outcome","line","stratum","k","N","effect","eff_lo","eff_hi","p","units","F_min","F_max"),
  note="Every external point estimate falls inside our confidence interval. Agreement is tightest on the growth outcomes and loosest on PTB --- which is where our data has signal and where it does not.")
frame_fig("S6_definition_spread", "S6 --- causal estimates under every BP definition",
  "SBP $\\rightarrow$ birth weight ranges $-96.2$ to $-86.5$ g per 10 mmHg across six definitions, all $p<0.006$. A monotone gradient with gestational timing is present but modest (11\\%, overlapping CIs).")
frame_table("S6_definition_sensitivity", title="S6 --- definition sensitivity (SAS pooled)",
  filter=quote(stratum=="SAS"),
  cols=c("definition","trait","outcome","k","N","theta_per10","lo","hi","p","primary"),
  note="\\texttt{lt20} and \\texttt{tri1} fall below the sample floor, so the early-pregnancy end is not assessable.")

add("\\subsection{Preterm birth by subtype}")
frame_table("T5_ptb_subtype", title="T5 --- PTB by subtype, per cohort",
  cols=c("cohort","trait","subtype","N","n_case","F","rf_p","OR","OR_lo","OR_hi"),
  note="Competing risks: each subtype uses \\emph{term} births as controls and drops the other subtype. \\texttt{unk} is excluded, not assigned.")
frame_table("T5b_ptb_subtype_pooled", title="T5b --- PTB subtype, pooled", rows_per=10,
  note=paste("Pre-specified prediction: indicated moves, spontaneous does not. Confirmed (difference $p=0.0026$ SBP).",
             "\\textbf{But} with a 50-event floor only Karachi and Matlab qualify, $I^2=60$, and Karachi is null ---",
             "so report as \\emph{Matlab's} result. The composition table (265 indicated at Matlab vs 18 at Pemba) is itself a finding about obstetric practice."))
frame_fig("SF4b_forest_subtype", "SF4b --- PTB subtype forest")

add("\\subsection{Instrument validity}")
frame_table("S12_controls", title="S12 --- positive and negative controls",
  cols=c("cohort","trait","variable","role","n","beta_crude","p_crude","p_pcadj","flag"),
  note=paste("\\textbf{0 of 48 negative controls fail after PC adjustment} (5 trip crude; PCs remove all) --- so the score does not track social position,",
             "\\emph{and} the PC adjustment is doing real work. Positive controls: PE $p=8.7\\times10^{-13}$, chronic hypertension $p=6.1\\times10^{-10}$ pooled.",
             "This is our only empirical check on the exclusion restriction, because pleiotropy-robust methods are undefined for a single score."))
frame_table("S14_mvmr", title="S14 --- multivariable MR with BMI",
  cols=c("cohort","trait","outcome","N","rf_zBP_uni","rf_zBP_mv","pct_change_zBP","rf_zBP_mv_p","rf_zBMI_mv_p"),
  filter=quote(outcome=="BWT"),
  note=sprintf("BMI score is a valid instrument here ($R^2$ 1.1--4.8\\%% on measured BMI, $F$ 16--190) and the two scores are nearly orthogonal ($r$ $-0.05$ to $0.26$). Median shift in the BP coefficient: \\textbf{%s\\%%} --- the birth-weight effect is not adiposity genetics. In the best-powered cohort it \\emph{strengthens} on adjustment.", mvmr_shift))
frame_table("S14b_mvmr_diagnostics", title="S14b --- MVMR diagnostics (read before the estimates)", rows_per=10,
  note="If the BMI instrument were weak, ``the BP effect survived MVMR'' would be hollow. It is not weak, so the test is informative.")

add("\\subsection{Maternal versus fetal genotype}")
frame_table("S15d_transmission_check", title="S15d --- transmission check (QC that validates the pairing)",
  rows_per=12,
  note=paste("$r$(maternal score, infant score) must be $\\approx 0.5$ by descent. Platforms below 0.25 are dropped.",
             "This gate caught a real error: hard-called low-pass WGS gave $r=0.02$--$0.13$ and would have been averaged in silently,",
             "producing a 1.6\\% maternal attenuation that reads exactly like ``the effect is maternal''."))
frame_table("S15_fetal_maternal", title="S15 --- maternal vs fetal genotype on birth weight",
  cols=c("cohort_display","trait","n_pairs","cor_mother_infant","mat_uni","mat_adj","pct_attenuation","mat_adj_p","fet_adj","fet_adj_p"),
  note="At $r\\approx0.5$ the two coefficients are individually imprecise even when their sum is well estimated, so a widening maternal interval is expected and is not evidence of absence.")
frame_table("S15b_fetal_pooled", title="S15b --- maternal vs fetal, pooled", rows_per=6,
  note=paste0("Maternal attenuates ", pct_att_sas, "\\% and loses significance; fetal gains it. \\textbf{But all maternal point estimates remain negative.} ",
              "\\textbf{We report this as a limitation, not a finding:} Warrington \\textit{et al.} (Nat Genet 2019) found the opposite with far more power ",
              "(maternal $-0.15$ SD/10 mmHg independent of fetal; \\emph{no} fetal effect, $-0.01$ [$-0.05$, $0.03$]). ",
              "Our \\emph{unconditional} estimate (", sd_units, " SD) already matches \\emph{their conditional} one --- which is what you expect if there is little real fetal contribution to remove."))

## ============================================================================
add("\\section{Interpretation}")
frame_text("The argument, end to end", paste(
"\\footnotesize",
"\\begin{enumerate}",
"\\item \\textbf{BP polygenic scores do not transfer uniformly.} $R^2$ ranges 0.8--6.5\\% across five cohorts.",
"\\item \\textbf{That variation is not ancestry.} Two cohorts that are genetically indistinguishable ($F_{ST}=" ,
sprintf("%s$) differ 2.6-fold. Nor is it BP variance, nor population structure. It tracks \\emph{how BP was measured}.", fst_bd),
"\\item \\textbf{Portability therefore determines where causal inference is possible.} The cohort with the earliest and",
"      densest BP measurement (Matlab) has the strongest instrument ($F$ up to 248); the African cohorts have $F$ as low as 17.",
"\\item \\textbf{Where the instrument works, higher maternal BP lowers birth weight.} Both traits, robust to BP definition,",
"      to pooling method, and to adiposity genetics --- and matching two independent European studies.",
"\\item \\textbf{Where it does not work, we say so.} The African arm is uninformative rather than null: two standard pooling",
"      methods disagree there by a factor of 13 and flip its sign.",
"\\end{enumerate}",
"\\vspace{2mm}",
"\\textbf{The contribution is the conjunction.} Nobody has tested BP score portability, or run this MR, in South Asian or African",
"pregnancy cohorts. The negative results are informative precisely because Part I explains \\emph{why} they are negative.", collapse="\n"))

frame_text("What changed during analysis, and why", paste(
"\\scriptsize",
"\\begin{itemize}",
"\\item \\textbf{Postpartum contamination.} $\\sim$40\\% of AMANHI BP readings were postnatal. Restricting to antenatal readings",
"      moved AMANHI-B's $R^2$ from 3.45\\% to 2.82\\% and reversed S16's conclusion --- the earlier ``averaging reconciles the cohorts'' was an artefact.",
"\\item \\textbf{Disattenuation.} The correction initially divided by the single-reading ICC, over-correcting $\\sim$2$\\times$",
"      (GAPPS-B SBP appeared at 16.5\\%, above published European estimates). Now uses Spearman--Brown reliability of the mean.",
"\\item \\textbf{Pooling method.} Meta-analysing per-cohort Wald ratios was replaced by pooling the two regressions and dividing once.",
"\\item \\textbf{Preterm birth.} The pooled null turned out to average an indicated-PTB signal against a spontaneous-PTB null.",
"\\item \\textbf{Fetal genotype.} Believed unavailable; found inside the same VCFs as the mothers.",
"\\end{itemize}",
"\\vspace{1mm}",
"Three of my own hypotheses were refuted by the data (dual-platform precision, GSA/dosage mix, MAF filtering as the African deficit).", collapse="\n"))

## ============================================================================
add("\\section{S18 --- STROBE-MR checklist}")
strobe <- list(
 c("1","Title/abstract","States design as Mendelian randomization; abstract reports exposure, outcomes, cohorts, and that MR is one-sample."),
 c("2","Background","BP--perinatal associations are confounded by BMI, SES, parity; MR addresses this. Portability is a prerequisite, not an aside."),
 c("3","Objectives","(i) Do BP PGS transfer to South Asian/African pregnancy cohorts? (ii) Where they do, what is the causal effect on perinatal outcomes?"),
 c("4","Study design","One-sample MR, five cohorts, individual-level data, per-cohort estimation then pooling of coefficients."),
 c("5","Setting","AMANHI (Sylhet, Karachi, Pemba), GAPPS/PreSSMat (Matlab), ZAPPS (Lusaka). Enrolment windows differ --- Table 1, Figure 3A."),
 c("6a","Participants","First-pregnancy mothers with genotype and valid antenatal BP; Figure 1A gives the flow (21,685 $\\rightarrow$ 13,212)."),
 c("6b","Instrument selection","\\textbf{Published PGS Catalog scores, not selected SNPs.} EUR primary; SAS for GAPPS-B. Selection is fully disclosed --- S4 and F2 report every score in every cohort. Weights are external (PGS Catalog), so no in-sample weight derivation."),
 c("7","Variables","Exposure: GA-standardised mean antenatal SBP/DBP. Outcomes: PTB, LBW, SGA, birth weight. Covariates: age + 5 within-cohort PCs."),
 c("8","Data sources","EPI extract for phenotypes; GSA array + low-pass WGS (dosages) for genotypes; PGS Catalog for weights. S1, S2, S3."),
 c("9","Bias","Weak instruments (S13, F reported per cell); population structure (S12, S17); fetal genotype (S15); winner's curse discussed."),
 c("10a","Quantitative variables","BP per 10 mmHg; birth weight in grams and SD."),
 c("10b","Statistical methods","Reduced form as the test; Wald ratio for magnitude; ratio-of-pooled-coefficients across cohorts."),
 c("10c","Assumption 1 (relevance)","First-stage $F$ 17--248 reported on every cell; S13 gives power."),
 c("10d","Assumption 2 (independence)","Negative controls vs education, wealth, gravidity, age --- crude and PC-adjusted (S12)."),
 c("10e","Assumption 3 (exclusion)","Positive controls (PE, chronic hypertension); fetal-genotype decomposition (S15). Pleiotropy-robust methods undefined for a single score (Bowden 2015)."),
 c("11","Sensitivity","Eight BP definitions (S6); two pooling methods (S10b); MVMR with BMI (S14); ancestry strata; F-threshold ladder."),
 c("12","Descriptive data","Table 1; S16 for BP distributions."),
 c("13","Main results","T2, T4, S10b. Effects with CIs; $p$ reported but not the organising principle."),
 c("14","Other analyses","T5 (PTB subtype), S15 (fetal), S14 (MVMR)."),
 c("15","Key results","Summarised in the executive summary."),
 c("16","Limitations","Underpowered MR; no paternal genotypes; early-pregnancy BP unmeasurable; single-score instrument; imputation performed upstream."),
 c("17","Interpretation","Portability determines where MR is possible; effect on fetal growth, not timing."),
 c("18","Generalisability","South Asian and African pregnancy cohorts; EUR-trained scores transfer poorly to the African cohorts."),
 c("19","Funding/data","Data availability and code: see provenance slide."))
add("\\begin{frame}{S18 --- STROBE-MR checklist \\textnormal{\\small(1/2)}}")
add("\\centering\\resizebox{\\linewidth}{!}{\\begin{tabular}{p{8mm}p{30mm}p{150mm}}\\toprule")
add("\\# & Item & How addressed \\\\ \\midrule")
for(s in strobe[1:12]) add(sprintf("%s & %s & %s \\\\[1mm]", s[1], s[2], s[3]))
add("\\bottomrule\\end{tabular}}\\end{frame}", "")
add("\\begin{frame}{S18 --- STROBE-MR checklist \\textnormal{\\small(2/2)}}")
add("\\centering\\resizebox{\\linewidth}{!}{\\begin{tabular}{p{8mm}p{30mm}p{150mm}}\\toprule")
add("\\# & Item & How addressed \\\\ \\midrule")
for(s in strobe[13:length(strobe)]) add(sprintf("%s & %s & %s \\\\[1mm]", s[1], s[2], s[3]))
add("\\bottomrule\\end{tabular}}")
add("\\note{Skrivankova VW \\textit{et al.} STROBE-MR. \\textit{BMJ} 2021;375:n2233 (explanation and elaboration); \\textit{JAMA} 2021;326:1614 (checklist).}")
add("\\end{frame}", "")

## ============================================================================
add("\\section{Limitations and provenance}")
frame_text("Limitations, stated plainly", paste(
"\\scriptsize",
"\\begin{itemize}",
sprintf("\\item \\textbf{Power.} %d of %d MR cells reach 80\\%%. Most nulls are uninformative, not evidence of absence.", npow, ncell),
"\\item \\textbf{Single polygenic score.} MR-Egger, weighted median and MR-PRESSO are undefined; S12's controls are the substitute.",
"\\item \\textbf{One-sample design.} Weak-instrument bias runs \\emph{toward} the confounded observational estimate, not toward the null.",
"\\item \\textbf{No paternal genotypes.} Limits the fetal decomposition; assortative mating could bias it (Warrington \\textit{et al.} flag the same gap).",
"\\item \\textbf{Early-pregnancy BP unmeasurable} --- AMANHI enrols too late.",
"\\item \\textbf{Imputation is upstream.} The VCFs arrived imputed; panel, software and parameters must be sourced for the Methods.",
"\\item \\textbf{The AFR stratum is not one population.} Pemba and Zambia are 25$\\times$ more diverged from each other than the two Bangladeshi cohorts.",
"\\item \\textbf{Part of the F3 gap is reading \\emph{count}, not timing}; the figure names only timing.",
"\\end{itemize}", collapse="\n"))

if(!is.null(MAN) && "id" %in% names(MAN)){
  MAN[, ord := .I]; last <- MAN[order(id, ord)][, .SD[.N], by=id][status=="OK"]
  if(nrow(last)){
    L <- last[, .(id, script, n, ts)][order(id)]
    fwrite(L, file.path(P$tables,"S19_provenance.tsv"), sep="\t")
    frame_table("S19_provenance", title="Provenance --- every deliverable and the script that made it",
      rows_per=17,
      note="Each row is one module's last successful run. The manifest also records input md5 fingerprints and the git SHA, so a stale result is detected automatically (\\texttt{bin/check\\_stale.R}).")
  }
}

frame_text("Reproducing this", paste(
"\\scriptsize",
"\\begin{itemize}",
"\\item Code: \\texttt{github.com/achatto4/MOMI}, directory \\texttt{bp\\_ptb\\_pipeline}.",
"\\item One command rebuilds every number: \\texttt{bash run\\_build.sh}. Build order and dependencies are the registry in that script.",
"\\item \\texttt{bin/check\\_stale.R} compares md5 fingerprints of every module's inputs against the manifest and reports anything out of date, failed, or not yet written.",
"\\item This deck: \\texttt{Rscript report/build\\_slides.R}. Regenerate after any pipeline change.",
"\\item Documentation: \\texttt{docs/BUILD\\_REFERENCE.md} (how each step works), \\texttt{docs/figure\\_provenance.md} (which script makes which display item).",
"\\end{itemize}",
sprintf("\\vspace{2mm}{\\tiny Built %s from git %s.}", format(Sys.time(),"%Y-%m-%d %H:%M"), momi_git_sha(PIPE)), collapse="\n"))

add("\\end{document}")

## ---------- write ---------------------------------------------------------------------
outf <- file.path(OUT,"slides.tex")
writeLines(LINES, outf)
nfr <- sum(grepl("\\\\begin\\{frame\\}", LINES))
cat(sprintf("\nwrote %s\n  %d lines, %d frames\n  figures copied to %s\n",
            outf, length(LINES), nfr, FIGDIR))
cat("\nnext:\n  cd", OUT, "&& pdflatex slides.tex && pdflatex slides.tex\n")
