options(stringsAsFactors = F)
args <- commandArgs(trailingOnly = TRUE)
genes_file = args[1]
canc_data = args[2]
data_out_dir = args[3]
n_cores = as.numeric(args[4])
workdir = '/sc/orga/work/holtj02/Other_Projects/Reva_SurvAnalysis/TanyaHypothesis_Data'
setwd(workdir)
plots_out = paste(data_out_dir,canc_data,sep='/')
#plots_out = paste('/sc/orga/work/holtj02/Other_Projects/Reva_SurvAnalysis/TanyaHypothesis_Data/Filtered_Plots',canc_data,sep='/')
dir.create(plots_out)
out_dir = 'r_cnvs'

make_class = function(row,comp) {
    bool_test = FALSE
    if(comp=="Pval_D1AD1B"){
        bool_test = (row[1] == -1 && row[2] == -1)
    } else if(comp=="Pval_D1AD2B"){
        bool_test = (row[1] == -1 && row[2] == -2)
    } else if(comp=="Pval_D1AD3B"){
        bool_test = (row[1] == -1 && row[2] < 0)
    } else if(comp=="Pval_D2AD1B"){
        bool_test = (row[1] == -2 && row[2] == -1)
    } else if(comp=="Pval_D3AD1B"){
        bool_test = (row[1] < 0 && row[2] == -1)
    } else if(comp=="Pval_D2AD2B"){
        bool_test = (row[1] == -2 && row[2] == -2)
    } else if(comp=="Pval_D2AD3B"){
        bool_test = (row[1] == -2 && row[2] < 0)
    } else if(comp=="Pval_D3AD2B"){
        bool_test = (row[1] < 0 && row[2] == -2)
    } else if(comp=="Pval_D3AD3B"){
        bool_test = (row[1] < 0 && row[2] < 0)
    }
    if(bool_test){
        return(1)
    } else {
        return(0)
    }
}

#### anti-correllation Analysis
library(parallel)
library(survival)
library(foreach)
library(doParallel)

## Read in the files and make parallilzation lists for processing.
all_data = read.table(paste(out_dir,canc_data,sep='/'),sep='\t',header=T,check.names=F)
genes_of_interest = read.csv(genes_file,header=T)
combinations = as.data.frame(cbind(apply(genes_of_interest,2,function(x) match(x,colnames(all_data)))))
#sig_cutoff = 0.05 / dim(genes_of_interest)[1]
sig_cutoff = 0.05
chunk_size = ceiling(dim(genes_of_interest)[1]/n_cores)
start = 1
i = 1
l_to_apply = list()
while(start < dim(genes_of_interest)[1]){
    l_to_apply[[i]] = list(genes_of_interest[start:(start+chunk_size-1),],combinations[start:(start+chunk_size-1),])
    start = start + chunk_size
    i = i + 1
}

## Create the cluster
if(n_cores > (detectCores()-1)){
    print(paste(n_cores,"cores is greater than the available number of cores",detectCores()-1))
    print(paste("changing n_cores from",n_cores,"to",detectCores()-1))
    n_cores = detectCores()-1
} else {
    debug_file = paste(plots_out,"debug.txt",sep="/")
    cl = makeCluster(n_cores,outfile = debug_file)
    registerDoParallel(cl)

    ## Run the parallelization
    tests = c("Pval_D1AD1B","Pval_D1AD2B","Pval_D1AD3B","Pval_D2AD1B","Pval_D3AD1B","Pval_D2AD2B","Pval_D2AD3B","Pval_D3AD2B","Pval_D3AD3B")
    n_tot = dim(all_data)[1]

    foreach(x=l_to_apply) %dopar% {
        pid = Sys.getpid()
        row_start = rownames(x[[1]])[1]
        out_file = paste(plots_out,paste(canc_data,"pid",as.character(pid),"row_start",row_start,paste(strsplit(basename(genes_file),split='.',fixed=T)[[1]][-4],collapse="."),"anti-corr.csv",sep='_'),sep='/')
        cat(paste("Cancer_type","GeneA","GeneB","Ntot","%D1A","%D2A","%D3A","%D1B","%D2B","%D3B","Pval_D1AD1B","Pval_D1AD2B","Pval_D1AD3B","Pval_D2AD1B","Pval_D3AD1B","Pval_D2AD2B","Pval_D2AD3B","Pval_D3AD2B",paste("Pval_D3AD3B",sep=''),sep=','),append=T,file=out_file)
        #sink(paste(plots_out,paste(canc_data,"chunk",pid,"row_start",row_start,strsplit(basename(genes_files),split='.',fixed=T)[[1]][1],"anti-corr.csv",sep='_'),sep='/'))
        for(i in 1:dim(x[[1]])){
            if (is.na(x[[1]][i,1])){ # check if the data is NA ... this is easier than filtering the list of listed data frames
                next()
            }
            n_d1ad1b = length(which((all_data[,x[[2]][i,1]] == -1) & (all_data[,x[[2]][i,2]] == -1)))
            n_d1a = length(which(all_data[,x[[2]][i,1]]==-1))
            n_d1b = length(which(all_data[,x[[2]][i,2]]==-1))
            f_test_d1ad1b = fisher.test(matrix(as.numeric(c(n_d1ad1b,(n_d1a-n_d1ad1b),(n_d1b-n_d1ad1b),(n_tot-(n_d1a+n_d1b-n_d1ad1b)))),nrow=2),alternative='l')
            n_d1ad2b = length(which((all_data[,x[[2]][i,1]] == -1) & (all_data[,x[[2]][i,2]] == -2)))
            n_d2b = length(which(all_data[,x[[2]][i,2]]==-2))
            f_test_d1ad2b = fisher.test(matrix(as.numeric(c(n_d1ad2b,(n_d1a-n_d1ad2b),(n_d2b-n_d1ad2b),(n_tot-(n_d1a+n_d2b-n_d1ad2b)))),nrow=2),alternative='l')
            n_d1ad3b = length(which((all_data[,x[[2]][i,1]] == -1) & (all_data[,x[[2]][i,2]] <= -1)))
            n_d3b = n_d1b + n_d2b
            f_test_d1ad3b = fisher.test(matrix(as.numeric(c(n_d1ad3b,(n_d1a-n_d1ad3b),(n_d3b-n_d1ad3b),(n_tot-(n_d1a+n_d3b-n_d1ad3b)))),nrow=2),alternative='l')
            n_d2ad1b = length(which((all_data[,x[[2]][i,1]] == -2) & (all_data[,x[[2]][i,2]] == -1)))
            n_d2a = length(which(all_data[,x[[2]][i,1]]==-2))
            f_test_d2ad1b = fisher.test(matrix(as.numeric(c(n_d2ad1b,(n_d2a-n_d2ad1b),(n_d1b-n_d2ad1b),(n_tot-(n_d2a+n_d1b-n_d2ad1b)))),nrow=2),alternative='l')
            n_d3ad1b = length(which((all_data[,x[[2]][i,1]] <= -1) & (all_data[,x[[2]][i,2]] == -1)))
            n_d3a = n_d1a + n_d2a
            f_test_d3ad1b = fisher.test(matrix(as.numeric(c(n_d3ad1b,(n_d3a-n_d3ad1b),(n_d1b-n_d3ad1b),(n_tot-(n_d3a+n_d1b-n_d3ad1b)))),nrow=2),alternative='l')
            n_d2ad2b = length(which((all_data[,x[[2]][i,1]] == -2) & (all_data[,x[[2]][i,2]] == -2)))
            f_test_d2ad2b = fisher.test(matrix(as.numeric(c(n_d2ad2b,(n_d2a-n_d2ad2b),(n_d2b-n_d2ad2b),(n_tot-(n_d2a+n_d2b-n_d2ad2b)))),nrow=2),alternative='l')
            n_d2ad3b = length(which((all_data[,x[[2]][i,1]] == -2) & (all_data[,x[[2]][i,2]] <= -1)))
            f_test_d2ad3b = fisher.test(matrix(as.numeric(c(n_d2ad3b,(n_d2a-n_d2ad3b),(n_d3b-n_d2ad3b),(n_tot-(n_d2a+n_d3b-n_d2ad3b)))),nrow=2),alternative='l')
            n_d3ad2b = length(which((all_data[,x[[2]][i,1]] <= -1) & (all_data[,x[[2]][i,2]] == -2)))
            f_test_d3ad2b = fisher.test(matrix(as.numeric(c(n_d3ad2b,(n_d3a-n_d3ad2b),(n_d2b-n_d3ad2b),(n_tot-(n_d3a+n_d2b-n_d3ad2b)))),nrow=2),alternative='l')
            n_d3ad3b = length(which((all_data[,x[[2]][i,1]] <= -1) & (all_data[,x[[2]][i,2]] <= -1)))
            f_test_d3ad3b = fisher.test(matrix(as.numeric(c(n_d3ad3b,(n_d3a-n_d3ad3b),(n_d3b-n_d3ad3b),(n_tot-(n_d3a+n_d3b-n_d3ad3b)))),nrow=2),alternative='l')
            #if(f_test_d2ad2b$p.value < sig_cutoff){
            cat(paste(paste(x[[1]][i,1],x[[1]][i,2],n_tot,(n_d1a/n_tot)*100,(n_d2a/n_tot)*100,(n_d3a/n_tot)*100,(n_d1b/n_tot)*100,(n_d2b/n_tot)*100,(n_d3b/n_tot)*100,f_test_d1ad1b$p.value,f_test_d1ad2b$p.value,f_test_d1ad3b$p.value,f_test_d2ad1b$p.value,f_test_d3ad1b$p.value,f_test_d2ad2b$p.value,f_test_d2ad3b$p.value,f_test_d3ad2b$p.value,f_test_d3ad3b$p.value,sep=','),'\n',sep=''),append=T,file=out_file)
        }
    }
    stopCluster(cl)
}
if (FALSE) {
cat(paste("Cancer_type","GeneA","GeneB","Ntot","%D1A","%D2A","%D3A","%D1B","%D2B","%D3B","Pval_D1AD1B","Pval_D1AD2B","Pval_D1AD3B","Pval_D2AD1B","Pval_D3AD1B","Pval_D2AD2B","Pval_D2AD3B","Pval_D3AD2B","Pval_D3AD3B","AvD1AD1B","AvNoD1AD1B","Chisq_surv_D1AD1B","AvD1AD2B","AvNoD1AD2B","Chisq_surv_D1AD2B","AvD1AD3B","AvNoD1AD3B","Chisq_surv_D1AD3B","AvD2AD1B","AvNoD2AD1B","Chisq_surv_D2AD1B","AvD3AD1B","AvNoD3AD1B","Chisq_surv_D3AD1B","AvD2AD2B","AvNoD2AD2B","Chisq_surv_D2AD2B","AvD2AD3B","AvNoD2AD3B","Chisq_surv_D2AD3B","AvD3AD2B","AvNoD3AD2B","Chisq_surv_D3AD2B","AvD3AD3B","AvNoD3AD3B",paste("Chisq_surv_D3AD3B","\n",sep=''),sep=','))
}
