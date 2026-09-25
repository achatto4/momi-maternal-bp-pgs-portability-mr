## -------------------------------------------------------------------------------
## 06_create_tables_and_figures.R
##
## Builds the display items from the aggregate result tables written by 03, 04 and 05 and
## by genetic_pcs/, and from the audited participant-flow counts in data/aggregate/:
##   * Figure 1, the incremental R-squared of the four score families (40 cells);
##   * Figure 2, the forest plot of the Mendelian-randomization results: cohort estimates
##     as circles, fixed-effect pooled estimates as diamonds, and I-squared with the
##     Cochran's Q p value beneath each overall pooled estimate;
##   * the supplementary table of all cohort estimates with the first-stage regression
##     each one rests on, and of the pooled estimates;
##   * Supplementary Figure S2, the participant flow by study site;
##   * Supplementary Figure S3, the descriptive cross-cohort principal-component analysis
##     (drawn when the genetic_pcs/ cross-cohort output is present);
##   * Supplementary Figure S4, the pooled estimates under the four exposure definitions.
##
## Every figure also writes the plotted values as a table, at full precision, so the
## figure can be checked against the numbers it draws.
## -------------------------------------------------------------------------------
if (!exists("MOMI_ROOT")) MOMI_ROOT <- getwd()
if (!exists("config")) source(file.path(MOMI_ROOT,
  if (file.exists(file.path(MOMI_ROOT, "config.R"))) "config.R" else "config.example.R"))
source(file.path(MOMI_ROOT, "analysis", "00_functions.R"))

## numbers are read as text and converted by as.numeric(), which is correctly rounded
rd <- function(f){ d <- fread(f, colClasses = "character", na.strings = NULL)
  for(cc in names(d)){ v <- d[[cc]]; v[v %chin% c("", "NA")] <- NA_character_; nv <- suppressWarnings(as.numeric(v))
    set(d, j = cc, value = if(all(is.na(v) == is.na(nv))) nv else v) }; d }

suppressMessages(library(data.table))
COH5  <- c("AMANHI-Bangladesh", "AMANHI-Pakistan", "GAPPS-Bangladesh", "AMANHI-Pemba", "GAPPS-Zambia")
SHORT <- c("AMANHI-Bangladesh"="Sylhet", "AMANHI-Pakistan"="Karachi", "GAPPS-Bangladesh"="Matlab", "AMANHI-Pemba"="Pemba", "GAPPS-Zambia"="Lusaka")
LONG  <- c("AMANHI-Bangladesh"="AMANHI-Sylhet", "AMANHI-Pakistan"="AMANHI-Karachi", "GAPPS-Bangladesh"="PreSSMat-Matlab",
           "AMANHI-Pemba"="AMANHI-Pemba", "GAPPS-Zambia"="ZAPPS-Lusaka")
FAM4  <- c("European", "South Asian", "East Asian", "Diverse ancestry")
MINUS <- intToUtf8(8722)
mn <- function(s){ s <- gsub("-", MINUS, s, fixed=TRUE); Encoding(s) <- "UTF-8"; s }
u8 <- function(s){ Encoding(s) <- "UTF-8"; s }
SUP2 <- intToUtf8(178)
ex <- function(x) vapply(x, function(v){ if(is.na(v)) return(NA_character_); if(!is.finite(v)) return(as.character(v))
  for(d in 1:17){ z <- sprintf(paste0("%.", d, "g"), v); if(identical(as.numeric(z), v)) return(z) }; sprintf("%.17g", v) }, character(1), USE.NAMES=FALSE)
exact <- function(dt){ dt <- copy(dt); for(cn in names(dt)) if(is.double(dt[[cn]])) set(dt, j=cn, value=ex(dt[[cn]])); dt }
devs <- function(out, w, h, fun, ps=8){
  cairo_pdf(paste0(out, ".pdf"), width=w, height=h, pointsize=ps, bg="white"); fun(); invisible(dev.off())
  png(paste0(out, ".png"), width=w, height=h, units="in", res=600, type="cairo", bg="white", pointsize=ps); fun(); invisible(dev.off()) }

## ---- Figure 1: portability ---------------------------------------------------

fig1_portability <- function(R, out){
  R <- as.data.table(R)[role == "Primary portability panel" & family %in% FAM4]
  stopifnot(nrow(R) == 40L, !anyDuplicated(R[, .(cohort, trait, family)]))
  R[, pct := 100 * incR2]; if(any(R$pct < -1e-8)) stop("negative incremental R2")
  top <- ceiling(max(R$pct) * 2) / 2
  pal <- colorRampPalette(c("#FFFFFF", "#1F3864"))(1001)
  colf <- function(v) pal[1 + round(1000 * pmin(1, pmax(0, v / top)))]
  lab  <- function(v) sprintf("%.1f", pmax(v, 0))
  R[, `:=`(fill=colf(pct), label=lab(pct), txtcol=fifelse(pmax(pct, 0) / top > 0.58, "white", "black"))]
  yr <- c(0, 1, 2, 3.3, 4.3); names(yr) <- COH5
  draw <- function(){
    par(mar=c(0, 0, 0, 0), family="sans", xaxs="i", yaxs="i", lend=1)
    plot.new(); plot.window(xlim=c(0, 7.09), ylim=c(3.0, 0))
    X0 <- c(SBP=0.93, DBP=3.78); CW <- 0.60; CH <- 0.36; Y0 <- 0.42
    for(tr in c("SBP", "DBP")){
      x0 <- X0[[tr]]
      text(x0 - 0.02, 0.2, paste0(if(tr == "SBP") "A" else "B", "   ", if(tr == "SBP") "Systolic blood pressure" else "Diastolic blood pressure"),
           adj=c(0, 0.5), font=2, cex=1.12)
      for(c in COH5) for(j in seq_along(FAM4)){
        r <- R[cohort == c & trait == tr & family == FAM4[j]]
        xl <- x0 + (j - 1) * CW; yt <- Y0 + yr[[c]] * CH
        rect(xl, yt, xl + CW, yt + CH, col=r$fill, border="white", lwd=1.6)
        text(xl + CW / 2, yt + CH / 2, r$label, col=r$txtcol, cex=1.0)
      }
      yb <- Y0 + (yr[[5]] + 1) * CH
      text(x0 + (seq_along(FAM4) - 0.5) * CW, yb + 0.17, c("European", "South\nAsian", "East\nAsian", "Diverse\nancestry"), cex=0.93, adj=c(0.5, 1))
      text(x0 + 2 * CW, yb + 0.62, "PGS development ancestry", cex=1.0)
    }
    x0 <- X0[["SBP"]]
    for(c in COH5) text(x0 - 0.08, Y0 + (yr[[c]] + 0.5) * CH, SHORT[[c]], adj=c(1, 0.5), cex=1.0)
    for(g in list(list("South Asian", 1, 3), list("African", 4, 5))){
      ya <- Y0 + yr[[COH5[g[[2]]]]] * CH + 0.03; yz <- Y0 + (yr[[COH5[g[[3]]]]] + 1) * CH - 0.03
      segments(0.30, ya, 0.30, yz, col="#555555", lwd=0.8); text(0.22, (ya + yz) / 2, g[[1]], srt=90, cex=1.0, col="#333333")
    }

    bx <- 6.40; bw <- 0.11; by0 <- 0.62; by1 <- 2.18
    ys <- seq(by1, by0, length.out=301)
    for(k in 1:300) rect(bx, ys[k], bx + bw, ys[k + 1], col=colf(top * (k - 0.5) / 300), border=NA)
    rect(bx, by0, bx + bw, by1, border="black", lwd=0.5)
    tk <- pretty(c(0, top), n=4); tk <- tk[tk <= top + 1e-9]
    yt <- by1 - (by1 - by0) * tk / top
    segments(bx + bw, yt, bx + bw + 0.04, yt, lwd=0.5); text(bx + bw + 0.07, yt, format(tk), adj=c(0, 0.5), cex=0.88)
    text(bx + bw + 0.36, (by0 + by1) / 2, u8(paste0("Incremental R", SUP2, " (%)")), srt=90, cex=0.98)
  }
  devs(out, 7.09, 3.0, draw)
  V <- R[, .(panel=fifelse(trait == "SBP", "A", "B"), trait, cohort, cohort_label=SHORT[cohort], family, score_id, N, incR2, incR2_pct=pct, label_displayed=label)]
  fwrite(exact(V[order(panel, match(cohort, COH5), match(family, FAM4))]), paste0(out, "_values.tsv"), sep="\t")
  invisible(list(values=V, scale_top=top))
}

fig1_portability(rd(file.path(RESULTS, "pgs_portability_results.tsv")), file.path(FIGURES, "figure1_portability"))

## ---- Figure 2: Mendelian-randomization forest plot ---------------------------

fig2_forest <- function(C, P, out){
  C <- as.data.table(C); P <- as.data.table(P)
  ROWSPEC <- data.table(slot=1:8, kind=c("cohort", "cohort", "cohort", "pool_sub", "cohort", "cohort", "pool_sub", "pool_all"),
                        rkey=c(COH5[1:3], "South Asian", COH5[4:5], "African", "Overall"),
                        label=c(LONG[COH5[1:3]], "South Asian pooled", LONG[COH5[4:5]], "African pooled", "Overall pooled"))
  SECTIONS <- data.table(outcome=c("PTB", "LBW", "SGA", "BWT"),
                         title=c("Preterm birth", "Low birth weight (<2,500 g)", "Small for gestational age (<10th centile)", "Birth weight"),
                         logscale=c(TRUE, TRUE, TRUE, FALSE), ref=c(1, 1, 1, 0),
                         axlab=c(rep("Odds ratio per 10 mmHg", 3), "Grams per 10 mmHg"))
  get_row <- function(tr, oc, sp){
    bin <- oc != "BWT"
    r <- if(sp$kind == "cohort") C[bp_trait == tr & outcome == oc & cohort == sp$rkey] else P[bp_trait == tr & outcome == oc & population_group == sp$rkey]
    stopifnot(nrow(r) == 1L)
    est <- if(bin) r$or_per_10mmHg else r$grams_per_10mmHg; lo <- if(bin) r$or_lo95 else r$grams_lo95; hi <- if(bin) r$or_hi95 else r$grams_hi95
    txt <- if(bin) sprintf("%.2f (%.2f to %.2f)", est, lo, hi) else mn(sprintf("%.0f (%.0f to %.0f)", est, lo, hi))
    het <- if(sp$kind == "pool_all") u8(sprintf("I%s = %.0f%%; Q p %s", SUP2, r$I2, if(r$Q_p < 0.001) "< 0.001" else sprintf("= %.3f", r$Q_p))) else ""
    data.table(bp_trait=tr, outcome=oc, slot=sp$slot, kind=sp$kind, row_key=sp$rkey, row_label=sp$label, estimate=est, ci_lo=lo, ci_hi=hi,
               value_text=txt, heterogeneity_text=het, I2=if(sp$kind == "cohort") NA_real_ else r$I2, Q_p=if(sp$kind == "cohort") NA_real_ else r$Q_p,
               k_cohorts=if(sp$kind == "cohort") 1L else as.integer(r$k_cohorts))
  }
  V <- rbindlist(lapply(SECTIONS$outcome, function(oc) rbindlist(lapply(c("SBP", "DBP"), function(tr)
         rbindlist(lapply(seq_len(nrow(ROWSPEC)), function(j) get_row(tr, oc, as.list(ROWSPEC[j]))))))))

  nice_log <- c(0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100)
  LIM <- rbindlist(lapply(SECTIONS$outcome, function(oc){
    d <- V[outcome == oc]; s <- SECTIONS[outcome == oc]
    if(s$logscale){ lo <- min(d$estimate) / 1.25; hi <- max(d$estimate) * 1.25
      core <- d[kind != "cohort"]; lo <- min(lo, core$ci_lo); hi <- max(hi, core$ci_hi)
      L <- max(nice_log[nice_log <= lo]); H <- min(nice_log[nice_log >= hi]); if(!is.finite(L)) L <- lo; if(!is.finite(H)) H <- hi
      data.table(outcome=oc, lo=L, hi=H) }
    else { core <- d[kind != "cohort"]; lo <- min(core$ci_lo, d$estimate) - 60; hi <- max(core$ci_hi, d$estimate) + 60
      data.table(outcome=oc, lo=floor(lo / 100) * 100, hi=ceiling(hi / 100) * 100) } }))
  SECTIONS <- merge(SECTIONS, LIM, by="outcome", sort=FALSE)
  V <- merge(V, SECTIONS[, .(outcome, lo, hi, logscale, ref)], by="outcome", sort=FALSE)
  V[, `:=`(trunc_lo=ci_lo < lo, trunc_hi=ci_hi > hi)]
  if(any(V$estimate < V$lo | V$estimate > V$hi)) stop("a point estimate falls outside its axis")
  V[, `:=`(drawn_lo=pmax(ci_lo, lo), drawn_hi=pmin(ci_hi, hi))]
  ticks <- function(s){ if(s$logscale){ nice_log[nice_log >= s$lo * (1 - 1e-9) & nice_log <= s$hi * (1 + 1e-9)] }
                        else { t <- pretty(c(s$lo, s$hi), n=5); t[t >= s$lo & t <= s$hi] } }
  X_LAB <- 0.000; X_IND <- 0.013; BLOCK <- c(SBP=0.1450, DBP=0.5875); PLOTW <- 0.2450; VALOFF <- 0.2570
  INK <- "#1B3A5C"; REF <- "grey55"; AXC <- "grey25"; COHF <- "grey35"
  R_HEAD <- 1.45; S_TITLE <- 1.15; S_AXIS <- 2.35; S_GAP <- 1.05; ROWH <- c(rep(1, 7), 1.55)
  sec_y0 <- numeric(4); y <- R_HEAD + 0.55
  for(i in 1:4){ sec_y0[i] <- y; y <- y + S_TITLE + sum(ROWH) + S_AXIS + S_GAP }
  TOTAL <- y - S_GAP + 0.35
  rowy <- function(y0, j) y0 + S_TITLE + c(0, cumsum(ROWH))[j] + 0.5
  mapx <- function(v, s, xl){ S <- SECTIONS[s]
    z <- if(S$logscale) (log(v) - log(S$lo)) / (log(S$hi) - log(S$lo)) else (v - S$lo) / (S$hi - S$lo); xl + z * PLOTW }
  head_at <- function(x, y, dir, col, lw){ ux <- diff(grconvertX(c(0, 0.052), from="inches", to="user")); uy <- diff(grconvertY(c(0, 0.052), from="inches", to="user"))
    a <- 22 * pi / 180; lines(c(x - dir * ux * cos(a), x, x - dir * ux * cos(a)), c(y - uy * sin(a), y, y + uy * sin(a)), col=col, lwd=lw) }
  draw <- function(){
    par(mar=c(0.15, 0.15, 0.15, 0.15), bg="white", family="sans", lend=1, ljoin=1, xaxs="i", yaxs="i")
    plot.new(); plot.window(xlim=c(0, 1), ylim=c(TOTAL, 0))
    for(tr in c("SBP", "DBP")){ xl <- BLOCK[[tr]]
      text(xl + PLOTW / 2, R_HEAD - 0.55, if(tr == "SBP") "Systolic blood pressure" else "Diastolic blood pressure", font=2, cex=1.05)
      text(xl + VALOFF, R_HEAD - 0.55, "Estimate (95% CI)", adj=c(0, 0.5), cex=0.95) }
    for(s in 1:4){ y0 <- sec_y0[s]; oc <- SECTIONS$outcome[s]; S <- SECTIONS[s]
      text(X_LAB, y0, S$title, adj=c(0, 0.5), font=2, cex=1.02)
      for(j in 1:8){ sp <- ROWSPEC[j]; yy <- rowy(y0, j)
        text(if(sp$kind == "cohort") X_IND else X_LAB, yy, sp$label, adj=c(0, 0.5), cex=1, font=if(sp$kind == "pool_all") 2 else 1,
             col=if(sp$kind == "cohort") "grey20" else "black") }
      for(tr in c("SBP", "DBP")){ xl <- BLOCK[[tr]]; xr <- mapx(S$ref, s, xl)
        segments(xr, y0 + S_TITLE + 0.10, xr, y0 + S_TITLE + sum(ROWH) - 0.45, col=REF, lty=2, lwd=0.7)
        D <- V[outcome == oc & bp_trait == tr][order(slot)]
        for(j in 1:8){ d <- D[j]; yy <- rowy(y0, j)
          xa <- mapx(d$drawn_lo, s, xl); xb <- mapx(d$drawn_hi, s, xl); xe <- mapx(d$estimate, s, xl)
          col <- if(d$kind == "cohort") "black" else INK; lw <- if(d$kind == "pool_all") 1.15 else if(d$kind == "pool_sub") 0.95 else 0.8
          segments(xa, yy, xb, yy, col=col, lwd=lw)
          if(d$trunc_lo) head_at(xa, yy, -1, col, lw) else segments(xa, yy - 0.16, xa, yy + 0.16, col=col, lwd=lw)
          if(d$trunc_hi) head_at(xb, yy, +1, col, lw) else segments(xb, yy - 0.16, xb, yy + 0.16, col=col, lwd=lw)
          if(d$kind == "cohort") points(xe, yy, pch=21, bg=COHF, col="black", cex=0.62, lwd=0.6) else {
            hw <- if(d$kind == "pool_all") 0.0062 else 0.0050; hh <- if(d$kind == "pool_all") 0.36 else 0.30
            polygon(c(xe - hw, xe, xe + hw, xe), c(yy, yy - hh, yy, yy + hh), col=if(d$kind == "pool_all") INK else "white", border=INK,
                    lwd=if(d$kind == "pool_all") 1.0 else 0.85) }
          text(xl + VALOFF, yy, d$value_text, adj=c(0, 0.5), cex=0.95, font=if(d$kind == "pool_all") 2 else 1)
          if(nzchar(d$heterogeneity_text)) text(xl + VALOFF, yy + 0.78, d$heterogeneity_text, adj=c(0, 0.5), cex=0.80, col="grey25") }
        ya <- y0 + S_TITLE + sum(ROWH) + 0.35 - 0.55 + 0.55
        segments(xl, ya, xl + PLOTW, ya, col=AXC, lwd=0.7)
        tk <- ticks(S); xt <- mapx(tk, s, xl)
        segments(xt, ya, xt, ya + 0.20, col=AXC, lwd=0.7)
        text(xt, ya + 0.62, if(S$logscale) formatC(tk, format="fg") else mn(formatC(tk, format="d")), cex=0.88)
        text(xl + PLOTW / 2, ya + 1.45, S$axlab, cex=0.92) } } }
  devs(out, 7.09, 9.30, draw)
  setorderv(V, c("outcome", "bp_trait", "slot"))
  fwrite(exact(V[, .(outcome, bp_trait, slot, kind, row_key, row_label, k_cohorts, estimate, ci_lo, ci_hi, axis_lo=lo, axis_hi=hi,
                     axis_scale=fifelse(logscale, "logarithmic", "linear"), reference=ref, truncated_lo=trunc_lo, truncated_hi=trunc_hi,
                     drawn_lo, drawn_hi, value_text, heterogeneity_text, I2, Q_p)]), paste0(out, "_values.tsv"), sep="\t")
  invisible(V)
}

fig2_forest(rd(file.path(RESULTS, "pooling_input.tsv")), rd(file.path(RESULTS, "ancestry_pooled_mr_results.tsv")),
            file.path(FIGURES, "figure2_mr_forest"))

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

## ---- Supplementary Figure S2: participant flow by study site -----------------
local({
T   <- fread(file.path(MOMI_ROOT, "data", "aggregate", "figureS2_participant_flow_by_site.tsv"))
OUT <- file.path(FIGURES, "figureS2_participant_flow_by_site")

SITES <- c("Sylhet", "Karachi", "Matlab", "Pemba", "Lusaka", "GARBH-INi")
HEAD1 <- c(Sylhet = "Sylhet", Karachi = "Karachi", Matlab = "Matlab", Pemba = "Pemba", Lusaka = "Lusaka", `GARBH-INi` = "Gurugram", total = "All sites")
HEAD2 <- c(Sylhet = "AMANHI", Karachi = "AMANHI", Matlab = "PreSSMat", Pemba = "AMANHI", Lusaka = "ZAPPS", `GARBH-INi` = "GARBH-INi", total = "")
COLS  <- c(SITES, "total")
stopifnot(identical(T$kind, c("retained", "excluded", "retained", "excluded", "retained", "excluded", "retained", "excluded", "retained")))
for(i in c(3, 5, 7, 9)) stopifnot(all(unlist(T[i, ..COLS]) == unlist(T[i - 2, ..COLS]) - unlist(T[i - 1, ..COLS])))
cnt <- function(x) formatC(x, format = "d", big.mark = ",")
LABELS <- c("Women with a first recorded study pregnancy", "Excluded: known multiple pregnancy",
            "Women after excluding known multiples", "Excluded: no harmonized PGS data",
            "Women with harmonized PGS data", "Excluded: no valid antenatal blood pressure",
            "Women with valid antenatal blood pressure", "Excluded: genetic quality control and identity resolution",
            "Cleaned genetic-analysis base sample")
W_IN <- 7.09; H_IN <- 3.62
draw <- function(){
  par(mar = c(0, 0, 0, 0), bg = "white", family = "sans", lend = 1, ljoin = 1, xaxs = "i", yaxs = "i")
  plot.new(); plot.window(xlim = c(0, W_IN), ylim = c(H_IN, 0))
  X0 <- 2.52; CW <- (W_IN - X0 - 0.04) / length(COLS); xc <- X0 + CW * (seq_along(COLS) - 0.5)
  BW <- 0.60; BH <- 0.215; YH <- 0.20; Y0 <- 0.62; DY <- 0.355
  yr <- Y0 + DY * (0:8)
  for(j in seq_along(COLS)){
    text(xc[j], YH, HEAD1[[COLS[j]]], font = 2, cex = 0.86)
    if(nzchar(HEAD2[[COLS[j]]])) text(xc[j], YH + 0.17, HEAD2[[COLS[j]]], cex = 0.74, col = "grey30")
  }
  segments(X0 + CW * 6, 0.06, X0 + CW * 6, H_IN - 0.06, col = "grey70", lwd = 0.6)
  for(i in 1:9){
    ret <- T$kind[i] == "retained"
    text(0.04, yr[i], LABELS[i], adj = c(0, 0.5), cex = if(ret) 0.80 else 0.74, font = if(ret) 1 else 3,
         col = if(ret) "black" else "grey30")
    for(j in seq_along(COLS)){
      v <- T[[COLS[j]]][i]; col_ <- COLS[j]
      if(col_ == "GARBH-INi" && i >= 5){
        if(i == 5) text(xc[j], yr[i] + 0.10, "Phenotype data\nonly; no genetic\nanalyses", cex = 0.66, col = "grey35", font = 3)
        next }
      if(ret){
        rect(xc[j] - BW / 2, yr[i] - BH / 2, xc[j] + BW / 2, yr[i] + BH / 2, border = "black", lwd = 0.7,
             col = if(i == 9) "grey92" else "white")
        text(xc[j], yr[i], cnt(v), cex = 0.80, font = if(i == 9 || col_ == "total") 2 else 1)
        if(i < 9 && !(col_ == "GARBH-INi" && i >= 3)){
          arrows(xc[j], yr[i] + BH / 2, xc[j], yr[i + 2] - BH / 2, length = 0.045, angle = 22, lwd = 0.6, col = "grey25")
        }
      } else {
        text(xc[j] + 0.035, yr[i], cnt(v), adj = c(0, 0.5), cex = 0.72, font = 3, col = "grey30")
      }
    }
  }
}
cairo_pdf(paste0(OUT, ".pdf"), width = W_IN, height = H_IN, pointsize = 9, bg = "white"); draw(); invisible(dev.off())
png(paste0(OUT, ".png"), width = W_IN, height = H_IN, units = "in", res = 600, type = "cairo", bg = "white", pointsize = 9); draw(); invisible(dev.off())
cat("Supplementary Figure S2 written:", OUT, "\n")

})

## ---- Supplementary Figure S3: descriptive cross-cohort principal components --

figS3_cross_pca <- function(B, E, out, seed=20260925L){
  B <- as.data.table(B); E <- as.data.table(E)
  COL <- c("AMANHI-Bangladesh"="#0072B2", "AMANHI-Pakistan"="#56B4E9", "GAPPS-Bangladesh"="#009E73", "AMANHI-Pemba"="#E69F00", "GAPPS-Zambia"="#D55E00")
  PCH <- c("low-pass WGS only"=16, "GSA only"=17, "both"=15)
  ex_ <- E[axis == "PC1"]; ey_ <- E[axis == "PC2"]; nb <- ex_$n_bins
  set.seed(seed)
  P <- B[rep(seq_len(.N), count)]
  P[, `:=`(x=ex_$lower + (bin_x - 1 + runif(.N)) * (ex_$upper - ex_$lower) / nb, y=ey_$lower + (bin_y - 1 + runif(.N)) * (ey_$upper - ey_$lower) / nb)]
  P <- P[sample(.N)]
  draw <- function(){
    par(mar=c(4.1, 4.4, 0.8, 11.5), mgp=c(2.5, 0.65, 0), las=1, family="sans", cex.axis=0.95, tcl=-0.3)
    plot(NA, xlim=c(ex_$lower, ex_$upper), ylim=c(ey_$lower, ey_$upper), xaxs="i", yaxs="i", bty="l", xaxt="n", yaxt="n",
         xlab=sprintf("PC1 (%.1f%% of variance)", ex_$pct_variance), ylab=sprintf("PC2 (%.1f%% of variance)", ey_$pct_variance))
    for(s in 1:2) axis(s, at=axTicks(s), labels=mn(format(axTicks(s), trim=TRUE)))
    points(P$x, P$y, pch=PCH[P$technology], cex=0.42, col=adjustcolor(COL[P$cohort], 0.6))
    u <- par("usr"); par(xpd=NA)
    legend(u[2] + 0.03 * diff(u[1:2]), u[4], title="Cohort", title.adj=0, legend=SHORT[COH5], pch=16, col=COL[COH5], bty="n", cex=0.95, pt.cex=1.1)
    legend(u[2] + 0.03 * diff(u[1:2]), u[4] - 0.50 * diff(u[3:4]), title="Genotyping technology", title.adj=0,
           legend=c("Low-pass WGS only", "GSA only", "Both"), pch=PCH, col="grey30", bty="n", cex=0.95, pt.cex=1.0)
  }
  devs(out, 7.09, 4.6, draw, ps=9)
  S <- B[, .(women=sum(count)), by=.(cohort, technology)][order(match(cohort, COH5), technology)]
  fwrite(S, paste0(out, "_values.tsv"), sep="\t")
  invisible(S)
}

XB <- file.path(config$pca_dir, "aggregate", "cross", "agg_cross_bins.tsv")
XE <- file.path(config$pca_dir, "aggregate", "cross", "agg_cross_bin_edges.tsv")
if(file.exists(XB) && file.exists(XE)) figS3_cross_pca(rd(XB), rd(XE), file.path(FIGURES, "figureS3_cross_cohort_pca")) else
  cat("06: cross-cohort PCA output not found in config$pca_dir; Supplementary Figure S3 not drawn\n")

## ---- Supplementary Figure S4: exposure definitions ---------------------------

figS4_exposure <- function(P, out){
  P <- as.data.table(P)
  DEFS <- c("mean", "ge20", "resid", "ge20last")
  DEFLAB <- u8(c(mean="Overall mean antenatal BP", ge20=paste0("Mean at ", intToUtf8(8805), "20 weeks"), resid="Gestational-age-standardised",
                 ge20last=paste0("Last reading at ", intToUtf8(8805), "20 weeks")))
  names(DEFLAB) <- DEFS
  DEFPCH <- c(mean=16, ge20=17, resid=1, ge20last=5); DEFCOL <- c(mean="#1B3A5C", ge20="grey45", resid="#1B3A5C", ge20last="grey45")
  DOFF <- c(mean=-0.30, ge20=-0.10, resid=0.10, ge20last=0.30)
  ROWS <- data.table(slot=1:8, panel=c(rep("A", 6), rep("B", 2)), outcome=c("PTB", "PTB", "LBW", "LBW", "SGA", "SGA", "BWT", "BWT"),
                     bp_trait=rep(c("SBP", "DBP"), 4), outtext=c("Preterm birth", "", "Low birth weight", "", "Small for gestational age", "", "Birth weight", ""),
                     trtext=rep(c("Systolic", "Diastolic"), 4))
  V <- rbindlist(lapply(seq_len(nrow(ROWS)), function(i) rbindlist(lapply(DEFS, function(d){
    r <- P[exposure_definition == d & bp_trait == ROWS$bp_trait[i] & outcome == ROWS$outcome[i]]; stopifnot(nrow(r) == 1L)
    bin <- ROWS$outcome[i] != "BWT"
    data.table(slot=ROWS$slot[i], panel=ROWS$panel[i], outcome=ROWS$outcome[i], bp_trait=ROWS$bp_trait[i], exposure_definition=d,
               estimate=if(bin) r$or_per_10mmHg else r$grams_per_10mmHg, ci_lo=if(bin) r$or_lo95 else r$grams_lo95, ci_hi=if(bin) r$or_hi95 else r$grams_hi95,
               k_cohorts=r$k_cohorts) }))))
  stopifnot(all(is.finite(V$estimate)), all(V$k_cohorts == 5L))
  lims <- function(lo, hi, logscale, ref, minspan){ tr <- function(v) if(logscale) log(v) else v; inv <- function(v) if(logscale) exp(v) else v
    full <- c(min(tr(lo)), max(tr(hi))); bulk <- c(as.numeric(quantile(tr(lo), 0.10)), as.numeric(quantile(tr(hi), 0.90)))
    use <- if(diff(full) > 2.5 * diff(bulk)) bulk else full; pad <- 0.06 * diff(use); o <- c(use[1] - pad, use[2] + pad)
    inv(c(min(o[1], tr(ref) - minspan / 2), max(o[2], tr(ref) + minspan / 2))) }
  A <- lims(V[panel == "A"]$ci_lo, V[panel == "A"]$ci_hi, TRUE, 1, log(1.6)); Bl <- lims(V[panel == "B"]$ci_lo, V[panel == "B"]$ci_hi, FALSE, 0, 120)
  nice_log <- function(lo, hi){ k <- seq(floor(log10(lo)), ceiling(log10(hi))); t <- sort(unique(as.vector(outer(c(1, 1.5, 2, 3, 5, 7), 10^k))))
    t <- t[t >= lo & t <= hi]; if(!(1 %in% t) && lo <= 1 && hi >= 1) t <- sort(c(t, 1))
    while(length(t) > 7){ keep <- t == 1 | seq_along(t) %% 2 == 1; if(all(keep)) break; t <- t[keep] }; t }
  A_TICK <- nice_log(A[1], A[2]); B_TICK <- pretty(Bl, n=5); B_TICK <- B_TICK[B_TICK >= Bl[1] & B_TICK <= Bl[2]]
  V[, `:=`(trunc_lo=fifelse(panel == "A", ci_lo < A[1], ci_lo < Bl[1]), trunc_hi=fifelse(panel == "A", ci_hi > A[2], ci_hi > Bl[2]))]
  V[, `:=`(drawn_lo=fifelse(panel == "A", pmax(ci_lo, A[1]), pmax(ci_lo, Bl[1])), drawn_hi=fifelse(panel == "A", pmin(ci_hi, A[2]), pmin(ci_hi, Bl[2])))]
  if(any(V[panel == "A"]$estimate < A[1] | V[panel == "A"]$estimate > A[2]) || any(V[panel == "B"]$estimate < Bl[1] | V[panel == "B"]$estimate > Bl[2]))
    stop("a pooled point estimate falls outside its axis")
  X_OUT <- 0.000; X_TR <- 0.215; X_PL <- 0.300; X_PR <- 0.985
  mapx <- function(v, pan) if(pan == "A") X_PL + (log(v) - log(A[1])) / (log(A[2]) - log(A[1])) * (X_PR - X_PL) else X_PL + (v - Bl[1]) / (Bl[2] - Bl[1]) * (X_PR - X_PL)
  REF <- "grey55"; AXC <- "grey25"
  head_at <- function(x, y, dir, col, lw){ ux <- diff(grconvertX(c(0, 0.045), from="inches", to="user")); uy <- diff(grconvertY(c(0, 0.045), from="inches", to="user"))
    a <- 22 * pi / 180; lines(c(x - dir * ux * cos(a), x, x - dir * ux * cos(a)), c(y - uy * sin(a), y, y + uy * sin(a)), col=col, lwd=lw) }
  panel_ <- function(pan, lab, ticks, tlabs, ref, axtitle){
    D <- V[panel == pan]; slots <- sort(unique(D$slot)); nr <- length(slots)
    par(mar=c(2.6, 0.2, 1.0, 0.2), xaxs="i", yaxs="i"); plot.new(); plot.window(xlim=c(0, 1), ylim=c(nr + 0.75, 0.15))
    text(X_OUT, 0.42, lab, adj=c(0, 0.5), font=2, cex=1.15, xpd=NA)
    xr <- mapx(ref, pan); segments(xr, 0.72, xr, nr + 0.5, col=REF, lty=2, lwd=0.8)
    for(j in seq_along(slots)){ s <- slots[j]; rw <- ROWS[slot == s]; y <- j
      if(nzchar(rw$outtext)) text(X_OUT, y, rw$outtext, adj=c(0, 0.5), xpd=NA)
      text(X_TR, y, rw$trtext, adj=c(0, 0.5), xpd=NA)
      for(d in DEFS){ r <- D[slot == s & exposure_definition == d]; yy <- y + DOFF[[d]]
        xa <- mapx(r$drawn_lo, pan); xb <- mapx(r$drawn_hi, pan); xe <- mapx(r$estimate, pan)
        segments(xa, yy, xb, yy, col=DEFCOL[[d]], lwd=0.9)
        if(r$trunc_lo) head_at(xa, yy, -1, DEFCOL[[d]], 0.9) else segments(xa, yy - 0.10, xa, yy + 0.10, col=DEFCOL[[d]], lwd=0.9)
        if(r$trunc_hi) head_at(xb, yy, +1, DEFCOL[[d]], 0.9) else segments(xb, yy - 0.10, xb, yy + 0.10, col=DEFCOL[[d]], lwd=0.9)
        points(xe, yy, pch=DEFPCH[[d]], col=DEFCOL[[d]], bg="white", cex=0.66, lwd=0.9) } }
    ya <- nr + 0.55; segments(X_PL, ya, X_PR, ya, col=AXC, lwd=0.8, xpd=NA); xt <- mapx(ticks, pan)
    segments(xt, ya, xt, ya + 0.13, col=AXC, lwd=0.8, xpd=NA); text(xt, ya + 0.40, tlabs, cex=0.88, xpd=NA)
    text((X_PL + X_PR) / 2, ya + 0.95, axtitle, cex=0.95, xpd=NA) }
  draw <- function(){ par(bg="white", family="sans", lend=1, ljoin=1); layout(matrix(1:3, ncol=1), heights=c(6.6, 2.9, 1.05))
    panel_("A", "A", A_TICK, formatC(A_TICK, format="fg"), 1, "Odds ratio per 10 mmHg")
    panel_("B", "B", B_TICK, mn(formatC(B_TICK, format="d")), 0, "Grams per 10 mmHg")
    par(mar=c(0.1, 0.2, 0.1, 0.2)); plot.new(); plot.window(xlim=c(0, 1), ylim=c(0, 1)); xs <- c(0.02, 0.27, 0.52, 0.74)
    for(i in seq_along(DEFS)){ d <- DEFS[i]; points(xs[i], 0.55, pch=DEFPCH[[d]], col=DEFCOL[[d]], bg="white", cex=0.66, lwd=0.9)
      text(xs[i] + 0.018, 0.55, DEFLAB[[d]], adj=c(0, 0.5), cex=0.92) } }
  devs(out, 7.09, 5.10, draw, ps=9)
  fwrite(exact(V), paste0(out, "_values.tsv"), sep="\t"); invisible(V)
}

figS4_exposure(rd(file.path(RESULTS, "exposure_four_definitions_pooled.tsv")), file.path(FIGURES, "figureS4_exposure_definitions"))
cat("06: figures and tables written\n")
