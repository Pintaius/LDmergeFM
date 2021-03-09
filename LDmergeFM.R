# Load arguments
script.flags <- commandArgs(trailingOnly = TRUE)
script.flags.n <- length(script.flags)
LOCUS <- script.flags[1]
cor.format <-
  if (script.flags.n < 2) {
    "LDSTORE"
  } else {
    script.flags[2]
  }
neff.formula <-
  if (script.flags.n < 3) {
    "Willer2010"
  } else {
    script.flags[3]
  }

# Load required libraries (CRAN versions after 01/09/18 will do)
suppressMessages(library("methods"))
suppressMessages(library("readr"))
suppressMessages(library("dplyr"))
suppressMessages(library("purrr"))
suppressMessages(library("reshape2"))
suppressMessages(library("Matrix"))

# Assertions
if (script.flags.n > 3) {
  stop("Too many arguments for script to understand.")
}
if (!file.exists(paste0(LOCUS, ".ref"))) {
  stop("Effect allele reference file not found, check locus name.")
}
if (!(neff.formula %in% c("Willer2010", "Vukcevic2011"))) {
  stop("Invalid option for sample size formula.")
}
if (!(cor.format %in% c("LDSTORE", "PLINK"))) {
  stop("Invalid option for correlation matrix format.")
}

# Load list of SNPs in locus
SNP <-
  unlist(read.table(
    paste0(LOCUS, ".ref"),
    header = F,
    colClasses = c("character", "NULL")
  ))

# Load correlation files into a list
if (cor.format=="LDSTORE"){
corfile <- Sys.glob(paste0("*", LOCUS, ".cor.gz"))
if (length(corfile) < 2) {
  stop("Fewer than two correlation matrices found, check locus name.")
}
samplevec <-
  unlist(strsplit(corfile, paste0("_", LOCUS, ".cor.gz")))
LDlist <-
  lapply(corfile,
         read_table2,
         col_types = cols_only(
           RSID1 = "c",
           RSID2 = "c",
           correlation = "d"
         ))}
if (cor.format=="PLINK") {
  corfile <- Sys.glob(paste0("*", LOCUS, ".ld.gz"))
  if (length(corfile) < 2) {
    stop("Fewer than two correlation matrices found, check locus name.")
  }
  samplevec <-
    unlist(strsplit(corfile, paste0("_", LOCUS, ".ld.gz")))
  LDlist <-
    suppressWarnings(lapply(
      corfile,
      read_table2,
      col_types = cols_only(
        SNP_A = "c",
        SNP_B = "c",
        R = "d"
      )
    ))
  LDlist <-
    LDlist %>% map( ~ rename(
      .x,
      RSID1 = SNP_A,
      RSID2 = SNP_B,
      correlation = R
    ))
}
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
if (length(famfile) < 2) {
  stop("Fewer than two .fam files found, check locus name.")
}
famlist <-
  lapply(
    famfile,
    read.table,
    header = F,
    na.strings = c("NA", "-9"),
    colClasses = c("NULL", "NULL", "NULL", "NULL", "NULL", "integer")
  )
names(famlist) <- samplevec
famlist.pheno <-
  lapply(famlist, function(x)
    x[!is.na(x)]) # This transforms each element into a phenotype vector

# Calculate effective sample sizes
neffvec <- samplevec
for (sample in samplevec) {
  ntot <- dim(famlist[[sample]])[1] # Takes into account missing phenotypes if used to calculate R2
  ncas <- sum(famlist.pheno[[sample]] == 2)
  if(ncas==0){ncas<-round(ntot/2)}
  if(neff.formula=="Willer2010") {neff <- 4 / (1 / ncas + 1 / (ntot - ncas))}
  if(neff.formula=="Vukcevic2011") {neff <- ntot*(ncas/ntot)*(1-(ncas/ntot))}
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

# Generate heatmap
heat.pal <- colorRampPalette(c("#d73027", "#fc8d59", "#fee090","#ffffbf","#e0f3f8","#91bfdb","#4575b4"))(256)
png(filename=paste0(LOCUS, ".heatmap.png"),type="cairo",width=4800,height=3200,pointsize=16,res=300)
heatmap(LDmatrix, Colv = NA, Rowv = NA, col = col, scale = "none",symm=T,labRow = F,labCol = F)
dev.off()

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
