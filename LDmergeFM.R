# Load locus number
LOCUS <- commandArgs(trailingOnly = TRUE)

# Load required libraries (updated CRAN versions as of 01/09/18)
suppressMessages(library("methods"))
suppressMessages(library("readr"))
suppressMessages(library("dplyr"))
suppressMessages(library("purrr"))
suppressMessages(library("reshape2"))
suppressMessages(library("Matrix"))

# Load list of SNPs in locus
SNP <-
  unlist(read.table(
    paste0(LOCUS, ".ref"),
    header = F,
    colClasses = c("character", "NULL")
  ))

# Load correlation files into a list
corfile <- Sys.glob(paste0("*", LOCUS, ".cor.gz"))
samplevec <-
  unlist(strsplit(corfile, paste0("_", LOCUS, ".cor.gz")))
LDlist <-
  lapply(corfile,
         read_table2,
         col_types = cols_only(
           RSID1 = "c",
           RSID2 = "c",
           correlation = "d"
         ))
names(LDlist) <- samplevec

# Expand correlation items to cover all pairwise SNP combinations
expand.cor <- function(cortable, snplist) {
  snpvalid <- unique(c(cortable$RSID1, cortable$RSID2))
  corlist.element <-
    data.frame(RSID1 = SNP,
               RSID2 = SNP,
               stringsAsFactors = F)
  corlist.element <-
    full_join(corlist.element, cortable, by = c("RSID1", "RSID2"))
  corlist.element$correlation[is.na(corlist.element$correlation) &
                                (corlist.element$RSID1 %in% snpvalid &
                                   corlist.element$RSID2 %in% snpvalid)] <-
    0
  return(corlist.element)
}
LDlist <- lapply(LDlist, expand.cor, snplist = SNP)

# Load fam files into a list
famfile <- Sys.glob(paste0("*", LOCUS, ".fam"))
famlist <-
  lapply(
    famfile,
    read.table,
    header = F,
    na.strings = c("NA", "-9"),
    colClasses = c("NULL", "NULL", "NULL", "NULL", "NULL", "integer")
  )
names(famlist) <- samplevec
famlist <-
  lapply(famlist, function(x)
    x[!is.na(x)]) # This transforms each element into a phenotype vector

# Calculate effective sample sizes
neffvec <- samplevec
for (sample in samplevec) {
  ntot <- length(famlist[[sample]])
  ncas <- sum(famlist[[sample]] == 2)
  neff <- 4 / (1 / ncas + 1 / (ntot - ncas)) # METAL formula
  # neff <- ntot*(ncas/ntot)*(1-(ncas/ntot)) # Matti's formula
  neffvec[neffvec == sample] <- neff
}

# Create LD data frame with missing diagonal
LDframe <- data.frame(RSID1 = SNP,
                      RSID2 = SNP,
                      stringsAsFactors = F)
LDframe <-
  full_join(LDframe,
            reduce(LDlist, full_join, by = c("RSID1", "RSID2")),
            by = c("RSID1", "RSID2"))
rm(LDlist) # Cleanup large object
gc() # Reallocate memory

# Calculate weighted correlation average
LDframe$wcor <-
  apply(LDframe[, c(-1, -2)],
        1,
        weighted.mean,
        w = as.numeric(neffvec),
        na.rm = T)
LDframe$wcor[LDframe$RSID1 == LDframe$RSID2] <- 1
LDframe <- subset(LDframe, select = c("RSID1", "RSID2", "wcor"))

# Transform into LD correlation format
LDmatrix <-
  acast(LDframe, RSID1 ~ RSID2, value.var = "wcor", drop = F, fill = 0) # Fill might not be necessary but is a safeguard for a FINEMAP error caused if NAs are present
rm(LDframe) # Cleanup large object
gc() # Reallocate memory
LDmatrix <- LDmatrix[SNP, SNP] # Keep SNP order, creating the upper triangle
LDmatrix <- forceSymmetric(LDmatrix, "U")
LDmatrix <-
  as.matrix(LDmatrix) # Return actual correlation values in writable format

# Export matrix with LD format
write.table(
  LDmatrix,
  file = paste0(LOCUS, ".ld"),
  col.names = F,
  row.names = F
)

# Log files used in computing correlations
options(max.print = 1000000)
sink(paste0(LOCUS, ".samples.log"))
for (obj in ls(pattern = "corfile")) {
  cat(unname(get(obj)), sep = "\n")
}
sink()
sink(paste0(LOCUS, ".snp.log"))
for (obj in ls(pattern = "SNP")) {
  cat(unname(get(obj)), sep = "\n")
}
sink()
