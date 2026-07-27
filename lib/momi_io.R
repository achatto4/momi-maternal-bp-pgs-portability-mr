#!/usr/bin/env Rscript
# ============================================================
# momi_io.R — shared I/O, paths, run-context, and MANIFEST logging for the
# synchronized MOMI paper pipeline. Every deliverable module sources this so that:
#   * paths are computed one way (no hard-coded absolute paths),
#   * every table/figure is written through the same writer,
#   * every module appends ONE row to results/current/manifest.tsv with its
#     status (OK / FAIL / SKIP), N, key numbers, inputs, outputs, runtime, git SHA.
# That manifest + the per-module logs are what let us trace any error later.
#
#   usage (top of every module):
#     PIPE <- Sys.getenv("MOMI_PIPE", unset=".")
#     source(file.path(PIPE,"lib/momi_io.R"))
#     source(file.path(PIPE,"lib/momi_config.R"))
# ============================================================
suppressMessages(library(data.table))

## ---- CLI arg parser (shared by all modules) ----
momi_arg <- function(flag, default=NULL){
  a <- commandArgs(TRUE); i <- which(a==flag)
  if(length(i)) a[i+1] else default
}

## ---- locate the pipeline root ----
# priority: env MOMI_PIPE > --pipe arg > current dir. Must contain lib/momi_config.R.
momi_pipe <- function(){
  p <- Sys.getenv("MOMI_PIPE", unset="")
  if(!nzchar(p)) p <- momi_arg("--pipe", ".")
  p <- normalizePath(p, mustWork=FALSE)
  if(!file.exists(file.path(p,"lib","momi_config.R")))
    stop(sprintf("momi_pipe(): '%s' is not the pipeline root (no lib/momi_config.R). Set MOMI_PIPE.", p))
  p
}

## ---- canonical directory layout under results/current/ ----
# intermediates/  frozen RDS that every deliverable reads (the synchronization point)
# tables/         <ID>.tsv, one per table deliverable
# figures/        <ID>.png/.pdf, one per figure deliverable
# logs/           <ID>.log, per-module stdout+stderr (written by the driver)
# manifest.tsv    one row per deliverable run
momi_paths <- function(pipe=momi_pipe(), root=NULL){
  if(is.null(root)) root <- Sys.getenv("MOMI_RESULTS", unset=file.path(pipe,"results","current"))
  P <- list(
    pipe          = pipe,
    root          = root,
    intermediates = file.path(root,"intermediates"),
    tables        = file.path(root,"tables"),
    figures       = file.path(root,"figures"),
    logs          = file.path(root,"logs"),
    manifest      = file.path(root,"manifest.tsv"),
    provenance    = file.path(root,"PROVENANCE.txt"))
  for(d in c(P$root,P$intermediates,P$tables,P$figures,P$logs))
    dir.create(d, showWarnings=FALSE, recursive=TRUE)
  P
}

## ---- run context (reproducibility stamp) ----
momi_git_sha <- function(pipe=momi_pipe()){
  s <- tryCatch(system(sprintf("git -C '%s' rev-parse --short HEAD",pipe),
                       intern=TRUE, ignore.stderr=TRUE), error=function(e) NA_character_)
  if(length(s)) s[1] else NA_character_
}
momi_md5 <- function(path){
  if(is.null(path)||is.na(path)||!file.exists(path)) return(NA_character_)
  tryCatch(tools::md5sum(path)[[1]], error=function(e) NA_character_)
}

## ---- intermediates: read/write the frozen RDS everyone shares ----
momi_save_intermediate <- function(obj, name, P=momi_paths()){
  f <- file.path(P$intermediates, paste0(name,".rds")); saveRDS(obj, f); f
}
momi_read_intermediate <- function(name, P=momi_paths()){
  f <- file.path(P$intermediates, paste0(name,".rds"))
  if(!file.exists(f)) stop(sprintf("intermediate '%s' not built yet (%s). Run its build step first.", name, f))
  readRDS(f)
}
momi_has_intermediate <- function(name, P=momi_paths())
  file.exists(file.path(P$intermediates, paste0(name,".rds")))

## ---- table writer ----
momi_write_table <- function(dt, id, P=momi_paths()){
  f <- file.path(P$tables, paste0(id,".tsv")); fwrite(dt, f, sep="\t"); f
}

## ---- figure writer (all figures are ggplot; write png + pdf) ----
momi_save_fig <- function(plot, id, width=7, height=5, P=momi_paths()){
  suppressMessages(require(ggplot2))
  png <- file.path(P$figures, paste0(id,".png")); pdf <- file.path(P$figures, paste0(id,".pdf"))
  ggplot2::ggsave(png, plot, width=width, height=height, dpi=300)
  ggplot2::ggsave(pdf, plot, width=width, height=height)
  c(png, pdf)
}

## ---- input/output fingerprints: how a change UPSTREAM rolls down ----
# Every deliverable declares the intermediates it reads (its `inputs` string). We hash those
# files at run time and store the digests in the manifest. A later staleness check
# (bin/check_stale.R) recomputes them: if the digest recorded when a deliverable last ran no
# longer matches the intermediate on disk, that deliverable is STALE -- something upstream
# was rebuilt and this item has not caught up. Tokens that are not intermediates (raw EPI
# path, sscore dir, ref/ files) hash as "-" and are ignored.
momi_fingerprint <- function(tokens, P, dir=NULL){
  toks <- trimws(unlist(strsplit(paste(tokens, collapse=";"), ";")))
  toks <- toks[nzchar(toks)]
  if(!length(toks)) return("")
  h <- vapply(toks, function(t){
    cand <- c(file.path(P$intermediates, paste0(t, ".rds")),
              if(!is.null(dir)) file.path(dir, t) else NULL)
    f <- cand[file.exists(cand)]
    if(length(f)) substr(momi_md5(f[1]), 1, 8) else NA_character_
  }, "")
  paste(sprintf("%s:%s", toks, ifelse(is.na(h), "-", h)), collapse=";")
}

## ---- manifest append (the traceability record) ----
# key = a short human-readable summary string of the headline numbers for this item.
momi_manifest_append <- function(id, script, status, n=NA, key="", inputs="",
                                 outputs="", seconds=NA, msg="", P=momi_paths()){
  # NB: build via a named list, NOT data.table(key=...) — `key` is a reserved
  # data.table() argument (sort key) and would be misread as column names.
  row <- as.data.table(list(
    ts       = format(Sys.time(),"%Y-%m-%d %H:%M:%S"),
    id       = id,
    script   = script,
    status   = status,               # OK | FAIL | SKIP
    n        = if(is.null(n)||length(n)==0) NA else n,
    key      = key,
    inputs   = inputs,
    outputs  = outputs,
    # digests of the intermediates this item READ, and of the files it WROTE. These are what
    # make an upstream change roll down: check_stale.R compares in_md5 against the files on
    # disk now, so rebuilding B01 immediately marks every downstream item stale.
    in_md5   = momi_fingerprint(inputs,  P),
    out_md5  = momi_fingerprint(outputs, P, dir=NULL),
    seconds  = round(as.numeric(seconds),1),
    git      = momi_git_sha(P$pipe),
    message  = gsub("[\r\n\t]+"," ",msg)))
  fwrite(row, P$manifest, sep="\t", append=file.exists(P$manifest))
  invisible(row)
}

## ---- deliverable wrapper: run one module body with timing + manifest logging ----
# Every deliverable calls momi_deliverable("<ID>", function(ctx){ ... }, script=...)
# The body returns a list(n=, key=, outputs=) which is logged. Errors are caught,
# logged as FAIL with the message, and (by default) re-raised so the driver sees a
# non-zero exit — set stop_on_error=FALSE for optional/Stage-7 items to log SKIP/FAIL
# without aborting the whole run.
momi_deliverable <- function(id, body, script=NA_character_, inputs="",
                             stop_on_error=TRUE, P=momi_paths()){
  t0 <- Sys.time()
  cat(sprintf("\n===== [%s] %s =====\n", id, script))
  res <- tryCatch(body(list(P=P)), error=function(e) structure(list(msg=conditionMessage(e)), class="momi_err"))
  secs <- as.numeric(difftime(Sys.time(), t0, units="secs"))
  if(inherits(res,"momi_err")){
    momi_manifest_append(id, script, "FAIL", inputs=inputs, seconds=secs, msg=res$msg, P=P)
    cat(sprintf("----- [%s] FAIL (%.1fs): %s\n", id, secs, res$msg))
    if(stop_on_error) stop(sprintf("[%s] failed: %s", id, res$msg))
    return(invisible(res))
  }
  if(is.list(res) && isTRUE(res$skip)){
    momi_manifest_append(id, script, "SKIP", inputs=inputs, seconds=secs, msg=res$reason %||% "", P=P)
    cat(sprintf("----- [%s] SKIP (%.1fs): %s\n", id, secs, res$reason %||% ""))
    return(invisible(res))
  }
  n   <- if(is.list(res)) res$n   else NA
  key <- if(is.list(res)) res$key else ""
  out <- if(is.list(res)) paste(res$outputs, collapse=";") else ""
  momi_manifest_append(id, script, "OK", n=n, key=key, inputs=inputs,
                       outputs=out, seconds=secs, P=P)
  cat(sprintf("----- [%s] OK (%.1fs) n=%s  %s\n", id, secs, as.character(n), key))
  invisible(res)
}

`%||%` <- function(a,b) if(is.null(a)||length(a)==0||(length(a)==1&&is.na(a))) b else a

invisible(TRUE)
