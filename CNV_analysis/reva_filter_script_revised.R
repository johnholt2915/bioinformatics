options(stringsAsFactors=F)
args <- commandArgs(trailingOnly = TRUE)
out_dir = 'r_cnvs'
workdir = '/sc/orga/work/holtj02/Other_Projects/Reva_SurvAnalysis/TanyaHypothesis_Data'
canc_data = args[1]
data_out_dir = args[2]
dir.create(data_out_dir)
dir.create(paste(data_out_dir,canc_data,sep='/'))
n_deletions = as.numeric(args[3]) ## values must be either c(-1, -2, -3)
number_deletions = as.numeric(args[4]) ## the higher this value, the more stringent ... recommend to set to 0.02
max_expected_score = as.numeric(args[5]) ## the lower this value, the more stringent ... recommended to set this value to 0.8
max_overlap_thresh = as.numeric(args[6]) ## the lower the value the more stringent ... recommended to set this value to 0.05

## Functions ##
exp_score = function(inds,sub_dat,n_tot){
    which_a = which(sub_dat[,inds[1]]==-2)
    which_b = which(sub_dat[,inds[2]]==-2)
    n_overlap = length(which(!is.na(match(which_a,which_b))))
    n_total_samples = length(which_a) + length(which_b) - n_overlap
    return(c((length(which_a)*length(which_b))/n_tot,n_overlap/n_total_samples))
}
###############

## load the cnv data and determine the number of deletions (single or double) per gene
all_data = read.table(paste(out_dir,canc_data,sep='/'),sep='\t',header=T,check.names=F)
cat(paste("loaded",paste(out_dir,canc_data,sep='/'),"\n"))
ntot = dim(all_data)[1]
if(n_deletions == -3){
    out = apply(all_data[,4:24779],2,function(x) return(length(which(x<0))/ntot))
} else {
    out = apply(all_data[,4:24779],2,function(x) return(length(which(x==n_deletions))/ntot))
}
cat(paste(ntot,"total samples\n"))
## Filter the genes by those having more than number_deletions
sub_data = all_data[,c(1,2,3,which(out>number_deletions)+3)]
cat(paste("filtered all genes by those with fewer than",number_deletions,"samples with",n_deletions,"deletions\n"))
cat(paste(dim(sub_data)[2]-3,"genes remaining\n"))

## find all pairwise combinations of gene pairs for applying further filtering steps
combos = combn(c(4:dim(sub_data)[2]),2)
#print(paste("generated all pair-wise combinations of genes passing the first filter"))
cat(paste(dim(combos)[2],"gene pairs remaining\n"))

## Calculate the Expected score ((Number of samples with n_deletions of GeneA) * (Number of samples with n_deletions of GeneB)) / (Number of total samples)
## Calculate the percentage of all samples with n_deletions deletions which overlap between both genes of a gene_pair
if(dim(combos)[2] > 400000){
    library(parallel)
    library(survival)
    library(foreach)
    library(doParallel)
    n_cores = ceiling(dim(combos)[2] / 200000)
    if (detectCores() < n_cores){
        n_cores = detectCores()-2
    }
    chunk_size = ceiling(dim(combos)[2] / n_cores)
    l_to_apply = list()
    start = 1
    count = 1
    while(start < dim(combos)[2]){
        l_to_apply[[count]] = combos[,start:(start+chunk_size-1)]
        count = count + 1
        start = start + chunk_size
    }
    debug_file = paste("debug.log.txt",sep="/")
    cl = makeCluster(n_cores,outfile = debug_file)
    registerDoParallel(cl)
    out_2nd = foreach(x=l_to_apply) %dopar% {
        return(apply(x,2,function(x) exp_score(x,sub_data,ntot)))
    }
    out_2nd = do.call("cbind",out_2nd)
    stopCluster(cl)
} else {
    out_2nd = apply(combos,2,function(x) exp_score(x,sub_data,ntot))
}
out_2nd = as.data.frame(t(out_2nd))
genes_to_keep = t(combos[,which(out_2nd$V1 < max_expected_score & out_2nd$V2 < max_overlap_thresh)])
cat(paste(dim(genes_to_keep)[1],"gene pairs remaining\n"))
genes_to_keep = apply(genes_to_keep,c(1,2),function(x) colnames(sub_data)[x])
write.csv(genes_to_keep,row.names=F,file=paste(data_out_dir,"/",canc_data,"/filtered_gene_pairs_",canc_data,"_ndel_",n_deletions,"_numdel_",number_deletions,"_mes_",max_expected_score,"_mot_",max_overlap_thresh,".csv",sep=""))
cat(paste("wrote output to",paste(data_out_dir,"/",canc_data,"/filtered_gene_pairs_",canc_data,"_ndel_",n_deletions,"_numdel_",number_deletions,"_mes_",max_expected_score,"_mot_",max_overlap_thresh,".csv\n",sep="")))


