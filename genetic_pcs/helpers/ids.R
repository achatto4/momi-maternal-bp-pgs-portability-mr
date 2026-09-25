suppressMessages(library(data.table))
normid <- function(x) { x <- as.character(x); if(!length(x)) return(character(0))
  as.character(fread(text=paste(c("IID", x), collapse="\n"), sep="\t", header=TRUE)$IID) }
rd_psam <- function(f){ x <- fread(f, colClasses="character"); setnames(x, 1, sub("^#", "", names(x)[1]))
  if(!"IID" %in% names(x)) setnames(x, 1, "IID"); x[, nid := normid(IID)]; x }
wr_keep <- function(psam, nids, f){ cols <- intersect(c("FID", "IID"), names(psam)); y <- psam[nid %chin% nids, ..cols]
  setnames(y, 1, paste0("#", names(y)[1])); fwrite(y, f, sep="\t"); nrow(y) }
rd_tab <- function(f, idcols){ idcols <- intersect(idcols, names(fread(f, nrows=0))); x <- fread(f, colClasses=list(character=idcols))
  setnames(x, 1, sub("^#", "", names(x)[1]))
  for(c in sub("^#", "", idcols)) if(c %in% names(x)) set(x, j=c, value=normid(x[[c]])); x }
rd_kin0 <- function(f){ h <- names(fread(f, nrows=0)); idc <- intersect(h, c("#FID1", "FID1", "#IID1", "IID1", "FID2", "IID2"))
  x <- fread(f, colClasses=list(character=idc)); names(x) <- sub("^#", "", names(x))
  if(!nrow(x)) return(data.table(ID1=character(0), ID2=character(0), NSNP=numeric(0), IBS0_FRAC=numeric(0), KINSHIP=numeric(0)))
  x[, IBS0 := as.numeric(IBS0)]; if(any(x$IBS0 > 1, na.rm=TRUE)) x[, IBS0 := IBS0 / as.numeric(NSNP)]
  x[, .(ID1=normid(IID1), ID2=normid(IID2), NSNP=as.numeric(NSNP), IBS0_FRAC=IBS0, KINSHIP=as.numeric(KINSHIP))] }
