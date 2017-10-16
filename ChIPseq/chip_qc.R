args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=FALSE)
qc_experiment_csv = args[1]
source("https://bioconductor.org/biocLite.R")
biocLite("ChIPQC")
library(ChIPQC)
samples = read.csv(qc_experiment_csv,header=T)
# for each line in samples data frame, run individual QC and append QCd sample to list
# make the full QC object after.

