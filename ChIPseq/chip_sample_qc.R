args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=FALSE)
bam_file = as.character(args[1])
peaks_file = as.character(args[2])
name = as.character(args[3])
#qc_experiment_csv = args[1]
source("https://bioconductor.org/biocLite.R")
biocLite("ChIPQC")
library(ChIPQC)
sample = ChIPQCsamples(bam_file,peaks=peaks_file)
sample_list = list(sample)
#samples = read.csv(qc_experiment_csv,header=T)
names(sample_list) = c(name)
save(sample_list,file=paste("RData/",name,"_sample_qc.RData",sep=''))
# for each line in samples data frame, run individual QC and append QCd sample to list
# make the full QC object after.

