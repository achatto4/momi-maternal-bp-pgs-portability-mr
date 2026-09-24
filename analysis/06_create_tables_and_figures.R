## -------------------------------------------------------------------------------
## 06_create_tables_and_figures.R
##
## Builds the display items from the result tables written by 03, 04 and 05:
##   * the main forest plot of the Mendelian-randomization results, four outcome sections
##     with the systolic and diastolic panels sharing each section's axis, cohort
##     estimates as circles and pooled estimates as diamonds, binary outcomes on a
##     logarithmic odds-ratio axis and birth weight on a linear gram axis;
##   * the supplementary table of all cohort estimates with the first-stage regression
##     each one rests on, and of the pooled estimates;
##   * the supplementary figure comparing the four exposure definitions.
##
## Every figure also writes the plotted values as a table, at full precision, so the
## figure can be checked against the numbers it draws.
## -------------------------------------------------------------------------------
if (!exists("MOMI_ROOT")) MOMI_ROOT <- getwd()
if (!exists("config")) source(file.path(MOMI_ROOT,
  if (file.exists(file.path(MOMI_ROOT, "config.R"))) "config.R" else "config.example.R"))
source(file.path(MOMI_ROOT, "analysis", "00_functions.R"))

## ---- main forest plot --------------------------------------------------------
local({
OUTDIR <- FIGURES
C <- fread(file.path(RESULTS, "pooling_input.tsv"))
P <- fread(file.path(RESULTS, "ancestry_pooled_mr_results.tsv"))

ex <- function(x) vapply(x, function(v){
  if(is.na(v)) return(NA_character_); if(!is.finite(v)) return(as.character(v))
  for(d in 1:17){ z <- sprintf(paste0("%.", d, "g"), v); if(identical(as.numeric(z), v)) return(z) }
  sprintf("%.17g", v) }, character(1), USE.NAMES = FALSE)
exact <- function(dt){ dt <- copy(dt); for(cn in names(dt)) if(is.double(dt[[cn]])) set(dt, j=cn, value=ex(dt[[cn]])); dt }

ROWSPEC <- data.table(
  slot  = 1:8,
  kind  = c("cohort","cohort","cohort","pool_sub","cohort","cohort","pool_sub","pool_all"),
  rkey  = c("AMANHI-Bangladesh","AMANHI-Pakistan","GAPPS-Bangladesh","South Asian",
            "AMANHI-Pemba","GAPPS-Zambia","African","Overall"),
  label = c("AMANHI-Sylhet","AMANHI-Karachi","PreSSMat-Matlab","South Asian pooled",
            "AMANHI-Pemba","ZAPPS-Lusaka","African pooled","Overall pooled"))

SECTIONS <- data.table(
  outcome = c("PTB","LBW","SGA","BWT"),
  title   = c("Preterm birth", "Low birth weight (<2,500 g)",
              "Small for gestational age (<10th centile)", "Birth weight"),
  logscale = c(TRUE, TRUE, TRUE, FALSE),
  ref      = c(1, 1, 1, 0),

  lo       = c(0.20, 0.20, 0.08, -550),
  hi       = c(8.00, 25.0, 8.00,  500),
  axlab    = c("Odds ratio per 10 mmHg", "Odds ratio per 10 mmHg",
               "Odds ratio per 10 mmHg", "Grams per 10 mmHg"))
TICKS <- list(PTB = c(0.25, 0.5, 1, 2, 4, 8),
              LBW = c(0.25, 0.5, 1, 2, 5, 10, 20),
              SGA = c(0.1, 0.25, 0.5, 1, 2, 4, 8),
              BWT = c(-400, -200, 0, 200, 400))

MINUS <- intToUtf8(8722)
mn <- function(s){ s <- gsub("-", MINUS, s, fixed = TRUE); Encoding(s) <- "UTF-8"; s }
TLABS <- list(PTB = c("0.25","0.5","1","2","4","8"),
              LBW = c("0.25","0.5","1","2","5","10","20"),
              SGA = c("0.1","0.25","0.5","1","2","4","8"),
              BWT = mn(c("-400","-200","0","200","400")))
get_row <- function(tr, oc, spec){
  bin <- oc != "BWT"
  if(spec$kind == "cohort"){
    r <- C[bp_trait == tr & outcome == oc & cohort == spec$rkey]
    stopifnot(nrow(r) == 1L)
  } else {
    r <- P[bp_trait == tr & outcome == oc & population_group == spec$rkey]
    stopifnot(nrow(r) == 1L)
  }
  est <- if(bin) r$or_per_10mmHg else r$grams_per_10mmHg
  lo  <- if(bin) r$or_lo95      else r$grams_lo95
  hi  <- if(bin) r$or_hi95      else r$grams_hi95
  txt <- if(bin) sprintf("%.2f (%.2f to %.2f)", est, lo, hi)
         else    mn(sprintf("%.0f (%.0f to %.0f)", est, lo, hi))
  data.table(bp_trait = tr, outcome = oc, slot = spec$slot, kind = spec$kind,
             row_key = spec$rkey, row_label = spec$label,
             estimate = est, ci_lo = lo, ci_hi = hi, value_text = txt,
             k_cohorts = if(spec$kind == "cohort") 1L else r$k_cohorts)
}
V <- rbindlist(lapply(seq_len(nrow(SECTIONS)), function(si)
       rbindlist(lapply(c("SBP","DBP"), function(tr)
         rbindlist(lapply(seq_len(nrow(ROWSPEC)), function(j)
           get_row(tr, SECTIONS$outcome[si], as.list(ROWSPEC[j]))))))))
V <- merge(V, SECTIONS[, .(outcome, lo, hi, logscale, ref)], by = "outcome", sort = FALSE)
V[, `:=`(trunc_lo = ci_lo < lo, trunc_hi = ci_hi > hi)]
if(any(V$estimate < V$lo | V$estimate > V$hi))
  stop("a point estimate falls outside its axis; widen the section limits rather than clipping a point")
V[, `:=`(drawn_lo = pmax(ci_lo, lo), drawn_hi = pmin(ci_hi, hi))]
setorderv(V, c("outcome","bp_trait","slot"))
fwrite(exact(V[, .(outcome, section = SECTIONS$title[match(outcome, SECTIONS$outcome)], bp_trait, slot, kind,
             row_key, row_label, k_cohorts, estimate, ci_lo, ci_hi,
             axis_lo = lo, axis_hi = hi, axis_scale = fifelse(logscale, "logarithmic", "linear"),
             reference = ref, truncated_lo = trunc_lo, truncated_hi = trunc_hi,
             drawn_lo, drawn_hi, value_text)]),
       file.path(OUTDIR, "figure3_values.tsv"), sep = "\t")

X_LAB   <- 0.000; X_IND <- 0.013
BLOCK   <- c(SBP = 0.1450, DBP = 0.5875)
PLOTW   <- 0.2450
VALOFF  <- 0.2570

INK  <- "#1B3A5C"; BLACK <- "#000000"; REF <- "grey55"; AXC <- "grey25"; COHF <- "grey35"

R_HEAD <- 1.45
S_TITLE<- 1.15
S_AXIS <- 2.35
S_GAP  <- 1.05
sec_y0 <- numeric(nrow(SECTIONS)); y <- R_HEAD + 0.55
for(i in seq_len(nrow(SECTIONS))){ sec_y0[i] <- y; y <- y + S_TITLE + 8 + S_AXIS + S_GAP }
TOTAL <- y - S_GAP + 0.35

mapx <- function(v, s, xl) {
  z <- if(SECTIONS$logscale[s]) (log(v) - log(SECTIONS$lo[s])) / (log(SECTIONS$hi[s]) - log(SECTIONS$lo[s]))
       else (v - SECTIONS$lo[s]) / (SECTIONS$hi[s] - SECTIONS$lo[s])
  xl + z * PLOTW
}

GUARD <- new.env(); GUARD$bad <- character(0); GUARD$slack <- Inf
guard <- function(txt, x0, xlim, cex, what){
  if(!nzchar(txt)) return(invisible(NULL))
  s <- xlim - (x0 + strwidth(txt, cex = cex)); GUARD$slack <- min(GUARD$slack, s)
  if(s < 0.003) GUARD$bad <- c(GUARD$bad, sprintf("%s: '%s' (%.4f)", what, txt, s))
  invisible(NULL)
}

HEAD_IN <- 0.052; HEAD_ANG <- 22 * pi / 180
head_at <- function(x, y, dir, col, lw){
  ux <- diff(grconvertX(c(0, HEAD_IN), from = "inches", to = "user"))
  uy <- diff(grconvertY(c(0, HEAD_IN), from = "inches", to = "user"))
  bx <- x - dir * ux * cos(HEAD_ANG); by <- uy * sin(HEAD_ANG)
  lines(c(bx, x, bx), c(y - by, y, y + by), col = col, lwd = lw)
}

draw <- function(cex0){
  par(mar = c(0.15, 0.15, 0.15, 0.15), bg = "white", family = "sans", lend = 1, ljoin = 1, xaxs = "i", yaxs = "i")
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(TOTAL, 0))

  for(tr in c("SBP","DBP")){
    xl <- BLOCK[[tr]]
    text(xl + PLOTW/2, R_HEAD - 0.55,
         if(tr == "SBP") "Systolic blood pressure" else "Diastolic blood pressure",
         font = 2, cex = cex0 * 1.05)
    text(xl + VALOFF, R_HEAD - 0.55, "Estimate (95% CI)", adj = c(0, 0.5), cex = cex0 * 0.95)
  }

  for(s in seq_len(nrow(SECTIONS))){
    y0 <- sec_y0[s]; oc <- SECTIONS$outcome[s]
    text(X_LAB, y0, SECTIONS$title[s], adj = c(0, 0.5), font = 2, cex = cex0 * 1.02)
    guard(SECTIONS$title[s], X_LAB, BLOCK[["SBP"]] + PLOTW, cex0 * 1.02, "section title")

    for(j in 1:8){
      sp <- ROWSPEC[j]; yy <- y0 + S_TITLE + j - 0.5
      xx <- if(sp$kind == "cohort") X_IND else X_LAB
      text(xx, yy, sp$label, adj = c(0, 0.5), cex = cex0,
           font = if(sp$kind == "pool_all") 2 else 1,
           col = if(sp$kind == "cohort") "grey20" else BLACK)
      guard(sp$label, xx, BLOCK[["SBP"]] - 0.004, cex0, "row label")
    }

    for(tr in c("SBP","DBP")){
      xl <- BLOCK[[tr]]
      xr <- mapx(SECTIONS$ref[s], s, xl)
      segments(xr, y0 + S_TITLE + 0.10, xr, y0 + S_TITLE + 8.05, col = REF, lty = 2, lwd = 0.7)

      D <- V[outcome == oc & bp_trait == tr][order(slot)]
      for(j in 1:8){
        d <- D[j]; yy <- y0 + S_TITLE + j - 0.5
        xa <- mapx(d$drawn_lo, s, xl); xb <- mapx(d$drawn_hi, s, xl); xe <- mapx(d$estimate, s, xl)
        col <- if(d$kind == "cohort") BLACK else INK
        lw  <- if(d$kind == "pool_all") 1.15 else if(d$kind == "pool_sub") 0.95 else 0.8

        segments(xa, yy, xb, yy, col = col, lwd = lw)
        if(d$trunc_lo) head_at(xa, yy, -1, col, lw) else segments(xa, yy - 0.16, xa, yy + 0.16, col = col, lwd = lw)
        if(d$trunc_hi) head_at(xb, yy, +1, col, lw) else segments(xb, yy - 0.16, xb, yy + 0.16, col = col, lwd = lw)

        if(d$kind == "cohort"){
          points(xe, yy, pch = 21, bg = COHF, col = BLACK, cex = cex0 * 0.62, lwd = 0.6)
        } else {

          hw <- if(d$kind == "pool_all") 0.0062 else 0.0050
          hh <- if(d$kind == "pool_all") 0.36   else 0.30
          polygon(c(xe - hw, xe, xe + hw, xe), c(yy, yy - hh, yy, yy + hh),
                  col = if(d$kind == "pool_all") INK else "white",
                  border = INK, lwd = if(d$kind == "pool_all") 1.0 else 0.85)
        }
        text(xl + VALOFF, yy, d$value_text, adj = c(0, 0.5), cex = cex0 * 0.95,
             font = if(d$kind == "pool_all") 2 else 1)
        guard(d$value_text, xl + VALOFF, xl + VALOFF + 0.16, cex0 * 0.95, "value text")
      }

      ya <- y0 + S_TITLE + 8.35
      segments(xl, ya, xl + PLOTW, ya, col = AXC, lwd = 0.7)
      tk <- TICKS[[oc]]; xt <- mapx(tk, s, xl)
      segments(xt, ya, xt, ya + 0.20, col = AXC, lwd = 0.7)
      text(xt, ya + 0.62, TLABS[[oc]], cex = cex0 * 0.88)
      text(xl + PLOTW/2, ya + 1.45, SECTIONS$axlab[s], cex = cex0 * 0.92)
    }
  }
}

W_IN <- 7.09; H_IN <- 9.06

cairo_pdf(file.path(OUTDIR, "figure3_mr_forest.pdf"), width = W_IN, height = H_IN, pointsize = 8, bg = "white")
draw(1.0); invisible(dev.off())
png(file.path(OUTDIR, "figure3_mr_forest.png"), width = W_IN, height = H_IN, units = "in",
    res = 600, type = "cairo", bg = "white", pointsize = 8)
draw(1.0); invisible(dev.off())

if(length(GUARD$bad))
  stop("figure layout: text would collide\n  ", paste(unique(GUARD$bad), collapse = "\n  "))
cat(sprintf("Figure 3 written to %s (%d plotted estimates; %d intervals truncated; tightest clearance %.4f)\n",
            OUTDIR, nrow(V), sum(V$trunc_lo | V$trunc_hi), GUARD$slack))
print(V[trunc_lo | trunc_hi, .(outcome, bp_trait, row_label, ci_lo = signif(ci_lo,3),
                               ci_hi = signif(ci_hi,3), trunc_lo, trunc_hi)])

})

## ---- supplementary table of all estimates ------------------------------------
local({
OUT <- TABLES
C <- fread(file.path(RESULTS, "pooling_input.tsv"))
P <- fread(file.path(RESULTS, "ancestry_pooled_mr_results.tsv"))

TRLAB <- c(SBP = "Systolic", DBP = "Diastolic")
OCLAB <- c(PTB = "Preterm birth", LBW = "Low birth weight", SGA = "Small for gestational age",
           BWT = "Birth weight")
COH   <- c("AMANHI-Bangladesh","AMANHI-Pakistan","GAPPS-Bangladesh","AMANHI-Pemba","GAPPS-Zambia")
COHLAB<- c("AMANHI-Sylhet","AMANHI-Karachi","PreSSMat-Matlab","AMANHI-Pemba","ZAPPS-Lusaka")
names(COHLAB) <- COH
mn  <- function(s){ s <- gsub("-", intToUtf8(8722), s, fixed = TRUE); Encoding(s) <- "UTF-8"; s }
num <- function(x, d) formatC(x, format = "f", digits = d, big.mark = "")

cnt <- function(x) trimws(formatC(as.integer(x), format = "d", big.mark = ","))
pf  <- function(p) fifelse(p < 0.001, "<0.001", num(p, 3))
est <- function(bin, e, lo, hi)
  fifelse(bin, sprintf("%s (%s to %s)", num(e, 2), num(lo, 2), num(hi, 2)),
               mn(sprintf("%s (%s to %s)", num(e, 0), num(lo, 0), num(hi, 0))))

C[, `:=`(tr_ord = match(bp_trait, c("SBP","DBP")), oc_ord = match(outcome, c("PTB","LBW","SGA","BWT")),
         co_ord = match(cohort, COH))]
setorder(C, tr_ord, oc_ord, co_ord)
bin <- C$outcome != "BWT"
A <- data.table(
  `Blood pressure`      = TRLAB[C$bp_trait],
  `Outcome`             = OCLAB[C$outcome],
  `Cohort`              = COHLAB[C$cohort],
  `Ancestry group`      = C$ancestry_group,
  `Polygenic score`     = C$pgs_id,
  `N`                   = cnt(C$N),
  `Events`              = fifelse(bin, cnt(C$n_events), "—"),
  `First-stage beta`    = num(C$fs_beta_mmHg_per_SD, 2),
  `First-stage SE`      = num(C$fs_se, 3),
  `First-stage F`       = num(C$fs_F, 1),
  `MR estimate (95% CI)`= est(bin, fifelse(bin, C$or_per_10mmHg, C$grams_per_10mmHg),
                              fifelse(bin, C$or_lo95, C$grams_lo95),
                              fifelse(bin, C$or_hi95, C$grams_hi95)),
  `P value`             = pf(C$wald_p))
fwrite(A, file.path(OUT, "tableS4_panelA_cohort.tsv"), sep = "\t")

P[, `:=`(tr_ord = match(bp_trait, c("SBP","DBP")), oc_ord = match(outcome, c("PTB","LBW","SGA","BWT")),
         gp_ord = match(population_group, c("South Asian","African","Overall")))]
setorder(P, tr_ord, oc_ord, gp_ord)
binp <- P$outcome != "BWT"
B <- data.table(
  `Blood pressure`         = TRLAB[P$bp_trait],
  `Outcome`                = OCLAB[P$outcome],
  `Population group`       = fifelse(P$population_group == "Overall", "Overall (all cohorts)",
                                     P$population_group),
  `Cohorts`                = as.character(P$k_cohorts),
  `N`                      = cnt(P$N_total),
  `Events`                 = fifelse(binp, cnt(P$events_total), "—"),
  `Pooled estimate (95% CI)` = est(binp, fifelse(binp, P$or_per_10mmHg, P$grams_per_10mmHg),
                                   fifelse(binp, P$or_lo95, P$grams_lo95),
                                   fifelse(binp, P$or_hi95, P$grams_hi95)),
  `P value`                = pf(P$p),
  `I2`                     = paste0(num(P$I2, 1), "%"),
  `Q P value`              = pf(P$Q_p))
fwrite(B, file.path(OUT, "tableS4_panelB_pooled.tsv"), sep = "\t")

cat(sprintf("Table S4: panel A %d rows, panel B %d rows\n", nrow(A), nrow(B)))
print(head(A, 3)); print(head(B, 3))

})

## ---- exposure-definition figure ----------------------------------------------
local({
OUTDIR <- FIGURES
P <- fread(file.path(RESULTS, "exposure_four_definitions_pooled.tsv"))

ex <- function(x) vapply(x, function(v){
  if(is.na(v)) return(NA_character_); if(!is.finite(v)) return(as.character(v))
  for(d in 1:17){ z <- sprintf(paste0("%.", d, "g"), v); if(identical(as.numeric(z), v)) return(z) }
  sprintf("%.17g", v) }, character(1), USE.NAMES = FALSE)
exact <- function(dt){ dt <- copy(dt); for(cn in names(dt)) if(is.double(dt[[cn]])) set(dt, j=cn, value=ex(dt[[cn]])); dt }

MINUS <- intToUtf8(8722)
mn <- function(s){ s <- gsub("-", MINUS, s, fixed=TRUE); Encoding(s) <- "UTF-8"; s }

NEWDEF <- "ge20last"
DEFS   <- c("mean","ge20","resid", NEWDEF)
DEFLAB <- c(mean  = "Overall mean antenatal BP",
            ge20  = paste0("Mean at ", intToUtf8(8805), "20 weeks"),
            resid = "Gestational-age-standardised",
            ge20last = paste0("Last reading at ", intToUtf8(8805), "20 weeks"))
Encoding(DEFLAB) <- "UTF-8"
DEFPCH <- c(mean = 16, ge20 = 17, resid = 1, ge20last = 5)
DEFCOL <- c(mean = "#1B3A5C", ge20 = "grey45", resid = "#1B3A5C",
            ge20last = "grey45")
DOFF   <- c(mean = -0.30, ge20 = -0.10, resid = 0.10, ge20last = 0.30)

ROWS <- data.table(
  slot     = 1:8,
  panel    = c(rep("A", 6), rep("B", 2)),
  outcome  = c("PTB","PTB","LBW","LBW","SGA","SGA","BWT","BWT"),
  bp_trait = c("SBP","DBP","SBP","DBP","SBP","DBP","SBP","DBP"),
  outtext  = c("Preterm birth","", "Low birth weight","", "Small for gestational age","",
               "Birth weight",""),
  trtext   = c("Systolic","Diastolic","Systolic","Diastolic","Systolic","Diastolic",
               "Systolic","Diastolic"))

V <- rbindlist(lapply(seq_len(nrow(ROWS)), function(i) rbindlist(lapply(DEFS, function(d){
  r <- P[exposure_definition == d & bp_trait == ROWS$bp_trait[i] & outcome == ROWS$outcome[i]]
  if(!nrow(r)) stop("no pooled row for ", d, " ", ROWS$bp_trait[i], " ", ROWS$outcome[i])
  bin <- ROWS$outcome[i] != "BWT"
  data.table(slot=ROWS$slot[i], panel=ROWS$panel[i], outcome=ROWS$outcome[i],
             bp_trait=ROWS$bp_trait[i], exposure_definition=d,
             estimate = if(bin) r$or_per_10mmHg else r$grams_per_10mmHg,
             ci_lo    = if(bin) r$or_lo95       else r$grams_lo95,
             ci_hi    = if(bin) r$or_hi95       else r$grams_hi95,
             k_cohorts = r$k_cohorts, status = r$status)
}))))
MISS <- V[!is.finite(estimate) | !is.finite(ci_lo) | !is.finite(ci_hi)]
if(nrow(MISS)){
  fwrite(exact(MISS), file.path(OUTDIR, "exposure_four_definitions_missing.tsv"), sep="\t")
  cat("NOTE: ", nrow(MISS), " estimate(s) could not be plotted; listed in ",
      "supplementary_exposure_sensitivity_missing.tsv\n", sep="")
  print(MISS[, .(outcome, bp_trait, exposure_definition, k_cohorts, status)])
}

fin <- function(x) x[is.finite(x)]

lims <- function(lo, hi, logscale, ref, minspan){
  lo <- fin(lo); hi <- fin(hi)
  tr  <- function(v) if(logscale) log(v) else v
  inv <- function(v) if(logscale) exp(v) else v
  full <- c(min(tr(lo)), max(tr(hi)))
  bulk <- c(as.numeric(quantile(tr(lo), 0.10)), as.numeric(quantile(tr(hi), 0.90)))
  use  <- if(diff(full) > 2.5 * diff(bulk)) bulk else full
  pad  <- 0.06 * diff(use)
  out  <- c(use[1] - pad, use[2] + pad)
  out  <- c(min(out[1], tr(ref) - minspan/2), max(out[2], tr(ref) + minspan/2))
  inv(out)
}
A <- lims(V[panel == "A"]$ci_lo, V[panel == "A"]$ci_hi, TRUE,  1, log(1.6))
B <- lims(V[panel == "B"]$ci_lo, V[panel == "B"]$ci_hi, FALSE, 0, 120)
A_LO <- A[1]; A_HI <- A[2]; B_LO <- B[1]; B_HI <- B[2]

nice_log <- function(lo, hi){
  k <- seq(floor(log10(lo)), ceiling(log10(hi)))
  t <- sort(unique(as.vector(outer(c(1, 1.5, 2, 3, 5, 7), 10^k))))
  t <- t[t >= lo & t <= hi]
  if(!(1 %in% t) && lo <= 1 && hi >= 1) t <- sort(c(t, 1))
  while(length(t) > 7){
    keep <- t == 1 | seq_along(t) %% 2 == 1
    if(all(keep)) break
    t <- t[keep]
  }
  t
}
fmt_log <- function(t) formatC(t, format="fg")
A_TICK <- nice_log(A_LO, A_HI)
B_TICK <- pretty(c(B_LO, B_HI), n = 5); B_TICK <- B_TICK[B_TICK >= B_LO & B_TICK <= B_HI]
V[, `:=`(trunc_lo = fifelse(panel == "A", ci_lo < A_LO, ci_lo < B_LO),
         trunc_hi = fifelse(panel == "A", ci_hi > A_HI, ci_hi > B_HI))]
V[is.na(trunc_lo), trunc_lo := FALSE]; V[is.na(trunc_hi), trunc_hi := FALSE]
V[, `:=`(drawn_lo = fifelse(panel == "A", pmax(ci_lo, A_LO), pmax(ci_lo, B_LO)),
         drawn_hi = fifelse(panel == "A", pmin(ci_hi, A_HI), pmin(ci_hi, B_HI)))]
if(any(fin(V$estimate) < ifelse(V[is.finite(estimate)]$panel == "A", A_LO, B_LO)) ||
   any(fin(V$estimate) > ifelse(V[is.finite(estimate)]$panel == "A", A_HI, B_HI)))
  stop("a pooled point estimate falls outside its axis; widen the limits rather than clipping a point")
fwrite(exact(V), file.path(OUTDIR, "exposure_four_definitions_values.tsv"), sep="\t")

X_OUT <- 0.000; X_TR <- 0.215; X_PL <- 0.300; X_PR <- 0.985
mapx <- function(v, pan){
  if(pan == "A") X_PL + (log(v) - log(A_LO))/(log(A_HI) - log(A_LO)) * (X_PR - X_PL)
  else           X_PL + (v - B_LO)/(B_HI - B_LO) * (X_PR - X_PL)
}
INK <- "#1B3A5C"; REF <- "grey55"; AXC <- "grey25"

HEAD_IN <- 0.045; HEAD_ANG <- 22 * pi / 180
head_at <- function(x, y, dir, col, lw){
  ux <- diff(grconvertX(c(0, HEAD_IN), from="inches", to="user"))
  uy <- diff(grconvertY(c(0, HEAD_IN), from="inches", to="user"))
  lines(c(x - dir*ux*cos(HEAD_ANG), x, x - dir*ux*cos(HEAD_ANG)),
        c(y - uy*sin(HEAD_ANG), y, y + uy*sin(HEAD_ANG)), col=col, lwd=lw)
}

draw_panel <- function(pan, lab, ticks, tlabs, ref, axtitle, cex0){
  D <- V[panel == pan]; slots <- sort(unique(D$slot)); nr <- length(slots)
  par(mar = c(2.6, 0.2, 1.0, 0.2), xaxs="i", yaxs="i")
  plot.new(); plot.window(xlim=c(0,1), ylim=c(nr + 0.75, 0.15))
  text(X_OUT, 0.42, lab, adj=c(0,0.5), font=2, cex=cex0*1.15, xpd=NA)
  xr <- mapx(ref, pan)
  segments(xr, 0.72, xr, nr + 0.5, col=REF, lty=2, lwd=0.8)
  for(j in seq_along(slots)){
    s <- slots[j]; rw <- ROWS[slot == s]; y <- j
    if(nzchar(rw$outtext)) text(X_OUT, y, rw$outtext, adj=c(0,0.5), cex=cex0, xpd=NA)
    text(X_TR, y, rw$trtext, adj=c(0,0.5), cex=cex0, xpd=NA)
    for(d in DEFS){
      r <- D[slot == s & exposure_definition == d]
      if(!nrow(r) || !is.finite(r$estimate)) next
      yy <- y + DOFF[[d]]
      xa <- mapx(r$drawn_lo, pan); xb <- mapx(r$drawn_hi, pan); xe <- mapx(r$estimate, pan)
      segments(xa, yy, xb, yy, col=DEFCOL[[d]], lwd=0.9)
      if(r$trunc_lo) head_at(xa, yy, -1, DEFCOL[[d]], 0.9)
      else segments(xa, yy-0.10, xa, yy+0.10, col=DEFCOL[[d]], lwd=0.9)
      if(r$trunc_hi) head_at(xb, yy, +1, DEFCOL[[d]], 0.9)
      else segments(xb, yy-0.10, xb, yy+0.10, col=DEFCOL[[d]], lwd=0.9)
      points(xe, yy, pch=DEFPCH[[d]], col=DEFCOL[[d]], bg="white", cex=cex0*0.66, lwd=0.9)
    }
  }
  ya <- nr + 0.55
  segments(X_PL, ya, X_PR, ya, col=AXC, lwd=0.8, xpd=NA)
  xt <- mapx(ticks, pan)
  segments(xt, ya, xt, ya + 0.13, col=AXC, lwd=0.8, xpd=NA)
  text(xt, ya + 0.40, tlabs, cex=cex0*0.88, xpd=NA)
  text((X_PL + X_PR)/2, ya + 0.95, axtitle, cex=cex0*0.95, xpd=NA)
}

render <- function(cex0){
  par(bg="white", family="sans", lend=1, ljoin=1)
  layout(matrix(1:3, ncol=1), heights=c(6.6, 2.9, 1.05))
  draw_panel("A", "A", A_TICK, fmt_log(A_TICK), 1,
             "Odds ratio per 10 mmHg", cex0)
  draw_panel("B", "B", B_TICK, mn(formatC(B_TICK, format="d")), 0,
             "Grams per 10 mmHg", cex0)
  par(mar=c(0.1, 0.2, 0.1, 0.2)); plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1))
  xs <- c(0.02, 0.27, 0.52, 0.74)
  for(i in seq_along(DEFS)){
    d <- DEFS[i]
    points(xs[i], 0.55, pch=DEFPCH[[d]], col=DEFCOL[[d]], bg="white", cex=cex0*0.66, lwd=0.9)
    text(xs[i] + 0.018, 0.55, DEFLAB[[d]], adj=c(0,0.5), cex=cex0*0.92)
  }
}

W_IN <- 7.09; H_IN <- 5.10
cairo_pdf(file.path(OUTDIR, "figureS1_exposure_definitions.pdf"), width=W_IN, height=H_IN,
          pointsize=9, bg="white")
render(1.0); invisible(dev.off())
png(file.path(OUTDIR, "figureS1_exposure_definitions.png"), width=W_IN, height=H_IN, units="in",
    res=600, type="cairo", bg="white", pointsize=9)
render(1.0); invisible(dev.off())
cat(sprintf("four-definition exposure figure written to %s (%d estimates, %d not plottable)\n",
            OUTDIR, nrow(V), nrow(MISS)))

})
cat("06: figures and tables written\n")
