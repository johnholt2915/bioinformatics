options(stringsAsFactors=F)
library(parallel)
library(foreach)
library(doParallel)
args <- commandArgs(trailingOnly = TRUE)
workdir = '/sc/orga/work/holtj02/Other_Projects/Reva_SurvAnalysis/TanyaHypothesis_Data/r_cnvs'
setwd(workdir)
gp_file = '/sc/orga/projects/Signatures/Cocorrelation_Analysis/johns_gene_pair_lists/ndel_-2_percdel_0.02_mes_1.2_mot_0.06_coordinates.tsv'
gp_anticorr_pval_dir = '/sc/orga/projects/Signatures/Cocorrelation_Analysis/johns_gene_pair_lists/stats_out_revised'
genes_info = '/sc/orga/projects/Signatures/Cocorrelation_Analysis/johns_gene_pair_lists/individual_gene_histogram_chromosome_simplified.tsv'

gps = read.table(gp_file,header=F)
gps_keys = as.character(apply(gps[,c(1,2)],1,function(x) paste(x,collapse='|')))

n_cores = 60
chunk_size = ceiling(dim(gps)[1] / n_cores)

ind_genes = read.table(genes_info,header=F)
gene_info = as.data.frame(cbind(ind_genes[match(gps[,1],ind_genes[,1]),c(2,3,4,5)],ind_genes[match(gps[,2],ind_genes[,1]),c(2,3,4,5)]))
colnames(gene_info) = c("APaired-count","ACancer-count","APaired-Cancer-product","AChromosome","BPaired-count","BCancer-count","BPaired-Cancer-product","BChromosome")

gps = as.data.frame(cbind(gps,gene_info))

cancs = list.files('.')
abbr_cancs = as.character(sapply(cancs,function(x) strsplit(x,split='.',fixed=T)[[1]][1]))

get_n_patients = function(canc_data,gene_a,gene_b,surv=FALSE){
    col_a = match(gene_a,colnames(canc_data))
    col_b = match(gene_b,colnames(canc_data))
    count_a = which(canc_data[,col_a] < 0)
    count_b = which(canc_data[,col_b] < 0)
    count_ovr = length(which(!is.na(match(count_a,count_b))))
    tot = length(count_a) + length(count_b) - count_ovr
    return(c(length(count_a),length(count_b),count_ovr,tot))
}

calc_score = function(n_tot,p_value,n_patients){
    return((1/p_value) * (n_patients/n_tot))
}

header = c("GeneA","GeneB","Frequency-Observed","Cancer-count","APaired-count","ACancer-count","APaired-Cancer-product","AChromosome","BPaired-count","BCancer-count","BPaired-Cancer-product","BChromosome")
additional_header = c("Cancer","Anti-corr-pvalue","nda","ndb","nboth","nwithd","scores","n_patients")
weights = list()
for(i in 2:length(cancs)){
    print(abbr_cancs[i])
    print("Reading CNV data...")
    canc_data = read.table(cancs[i],sep='\t',header=T,check.names=F)
    n_tot = dim(canc_data)[1]
    canc_pval_data = read.table(paste(gp_anticorr_pval_dir,cancs[i],paste(abbr_cancs[i],'_ndel_-2_percdel_0.02_mes_1.2_mot_0.6.anti-corr.tsv',sep=''),sep='/'),header=F)
    print("Retrieving gene pair prevalence metrics...")
    l_to_apply = list()
    start = 1
    count = 1
    while(start < dim(gps)[1]) {
        if(start+chunk_size-1 <= dim(gps)[1]){
            l_to_apply[[count]] = canc_pval_data[start:(start+chunk_size-1),]
        } else {
            l_to_apply[[count]] = canc_pval_data[start:dim(gps)[1],]
        }
        count = count + 1
        start = start + chunk_size
    }
    cl = makeCluster(n_cores)
    registerDoParallel(cl)
    n_patients = foreach(x=l_to_apply) %dopar% {
        return(apply(x,1,function(y) get_n_patients(canc_data,y[2],y[3])))
    }
    n_patients = do.call("cbind",n_patients)
    n_patients = t(n_patients)
    canc_pval_data = as.data.frame(cbind(canc_pval_data,n_patients))
    print("Calculating genepair scores...")
    l_to_apply = list()
    start = 1
    count = 1
    while(start < dim(gps)[1]) {
        if(start+chunk_size-1 <= dim(gps)[1]){
            l_to_apply[[count]] = canc_pval_data[start:(start+chunk_size-1),]
        } else {
            l_to_apply[[count]] = canc_pval_data[start:dim(gps)[1],]
        }
        count = count + 1
        start = start + chunk_size
    }
    scores = foreach(x=l_to_apply) %dopar% {
        return(apply(x,1,function(y) calc_score(n_tot,as.numeric(y[4]),as.numeric(y[8]))))
    }
    scores = unlist(scores)
    stopCluster(cl)
    canc_pval_data = as.data.frame(cbind(canc_pval_data,scores))
    print("Reordering the genepairs...")
    pval_keys = as.character(apply(canc_pval_data[,c(2,3)],1,function(x) paste(x,collapse='|')))
    tmp = canc_pval_data[match(gps_keys,pval_keys),]
    tmp[,1] = rep(abbr_cancs[i],dim(tmp)[1])
    tmp = as.data.frame(cbind(tmp,rep(n_tot,dim(tmp)[1])))
    weights[[i]] = tmp[,c(1,4:dim(tmp)[2])]
}

for(i in 3:length(weights)){
    gps = as.data.frame(cbind(gps,weights[[i]]))
    header = c(header,additional_header)
}

ran = range(13,132,dim(weights[[1]])[2])
for(i in 1:length(ran)){
    colnames(gps)[ran[i]:(ran[i]+7)] = paste(rep(abbr_cancs[i],8),rep("_",8),colnames(gps)[ran[i]:(ran[i]+7)],sep='')
}

write.table(gps,file='/sc/orga/projects/Signatures/Cocorrelation_Analysis/johns_gene_pair_lists/merged_data.tsv',row.names=F,quote=F,sep='\t')

inds = seq(13,132,8)
fixes = list()
for(i in 1:length(inds)){
	fixes[[i]] = gps[,c(1,2,3,4,inds[i]:(inds[i]+7))]
	geom_ave = as.numeric(apply(fixes[[i]],1,function(x) exp((-1/as.numeric(x[3]))*sum(log(as.numeric(rep(x[6],x[3])),base=exp(1))))))
	fixes[[i]] = as.data.frame(cbind(fixes[[i]],geom_ave))
	big_m = as.numeric(apply(fixes[[i]],1,function(x) sum(as.numeric(x[10])*(1/as.numeric(rep(x[6],x[3]))))/sum((1/as.numeric(rep(x[6],x[3]))))))
	fixes[[i]] = as.data.frame(cbind(fixes[[i]],big_m))
}

for(i in 1:length(fixes)){
	tmp = colnames(fixes[[i]])
	canc = strsplit(tmp[5],split='_')[[1]][1]
	print(canc)
	print(paste(canc,tmp[13],sep='_'))
	tmp[13] = paste(canc,tmp[13],sep='_')
	print(paste(canc,tmp[14],sep='_'))
	tmp[14] = paste(canc,tmp[14],sep='_')
	colnames(fixes[[i]]) = tmp
}

fixes_gps = fixes[[1]]
for(i in fixes[2:length(fixes)]){
	fixes_gps = as.data.frame(cbind(fixes_gps,i[5:dim(i)[2]]))
}
write.table(fixes_gps,file='/sc/orga/projects/Signatures/Cocorrelation_Analysis/johns_gene_pair_lists/fixes_merged_data.tsv',row.names=F,quote=F,sep='\t')

sums = apply(fixes_gps,1,function(x) sum(as.numeric(x[seq(14,154,10)])))
fixes_gps = as.data.frame(cbind(fixes_gps,sums))

write.table(fixes_gps,file='/sc/orga/projects/Signatures/Cocorrelation_Analysis/johns_gene_pair_lists/fixes_merged_data.tsv',row.names=F,quote=F,sep='\t')


