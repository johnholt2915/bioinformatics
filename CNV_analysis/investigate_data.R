# go to data directory
setwd('/Users/holtj02/Desktop/Projects/Reva_Survivorship/TanyaHypothesis_Data')
cnv_dir = 'cnvs'
survivorship_dir = 'clean_clinical'
out_dir = 'r_cnvs'
plots_out = '/Users/holtj02/Desktop/Projects/Reva_Survivorship/TanyaHypothesis_Data/RB1_KDM5_plots'
#checking if cnv and survivorship data are all paired (looks like some might not be)
View(cbind(list.files(cnv_dir),c(list.files(survivorship_dir),"")))

"""
matches:
blca.cnv  — clinical_blca.new
Brca.cnv — clinical_brca.new
Coadread.cnv — clinical_coadread.new
Esca.cnv — clinical_esca.new
Kirc.cnv — clinical_kirc.new
Kirp.cnv — clinical_kirp.new
Lgg.cnv — clinical_lgg.new
Lihc.cnv — clinical_lihc.new
Ov.cnv — clinical_ov.new
Paad.cnv — clinical_paad.new
Prad.cnv — clinical_prad.new
Sarc.cnv — clinical_sarc.new
Skcm.cnv — clinical_skcm.new
Stad.cnv — clinical_stad.new
Ucec.cnv — clinincal_ucec.new

mismatches:
Gbm.cnv — clinical_ERpos.new
Hnsc.cnv — clinical_her2+.new
Laml.cnv — clinical_thca.new
Luad.cnv — clinical_tnbc.new
Lusc.cnv — ?
"""
# isolate cnv and survivorship data which are paired (this may change when Boris provides more info)
cancers_cnvs = as.vector(sapply(list.files(cnv_dir),function(x) strsplit(x,split='.',fixed = TRUE)[[1]][1]))
cancers_surv = as.vector(sapply(list.files(survivorship_dir),function(x) strsplit(strsplit(x,split='.',fixed=T)[[1]][1],split='_',fixed=T)[[1]][2]))
paired_cnvs = list.files(cnv_dir)[which(!is.na(match(cancers_cnvs,cancers_surv)))]
paired_survs = list.files(survivorship_dir)[which(!is.na(match(cancers_surv,cancers_cnvs)))]

# now lets try loading the cnv/survivorship datasets and merging them for easier loading.
#options(stringsAsFactors = F)
#for(ind in 1:length(paired_cnvs)){
  print(paired_cnvs[ind])
  print(paired_survs[ind])
  print("Reading Data")
  cnv = data.frame(t(read.table(paste(cnv_dir,paired_cnvs[ind],sep='/'), sep="\t", header = T,check.names = F,row.names = 'Gene')))
  surv = read.table(paste(survivorship_dir,paired_survs[ind],sep='/'),sep='\t',header=T)
  print("Formatting Data")
  surv$bcr_patient_barcode = toupper(surv$bcr_patient_barcode)
  cnv = data.frame(cbind(toupper(rownames(cnv)),cnv))
  colnames(cnv) = c("bcr_patient_barcode",colnames(cnv)[2:dim(cnv)[2]])
  all_data = merge(surv,cnv,by="bcr_patient_barcode")
  print("Writing Data")
  write.table(all_data,file = paste(out_dir,paste(paired_cnvs[ind],paired_survs[ind],sep='_'),sep='/'),sep='\t')
}

#### Survival Analysis
library(survival)
plots_out = 'RB1_KDM5_plots'
make_class = function(row) {
  if(row[1] == -2 && row[2] < 0){
    return(1)
  } else {
    return(0)
  }
}
sink(paste(plots_out,"RB1_KDM5_stats_out.txt",sep='/'))
for(ind in 1:length(paired_cnvs)){
  cat(paste("***********",paired_cnvs[ind],paired_survs[ind],"************\n",sep=' '))
  all_data = read.table(paste(out_dir,paste(paired_cnvs[ind],paired_survs[ind],sep='_'),sep='/'),sep='\t',header=T,check.names = F)
  genes_of_interest = c("RB1","KDM5A","KDM5B","KDM5C")
  genes_of_interest = t(combn(genes_of_interest,2))
  genes_of_interest = genes_of_interest[which(genes_of_interest[,1]=="RB1"),]
  # select only combinations which include RB1
  combinations = apply(genes_of_interest,2,function(x) match(x,colnames(all_data)))
  for(i in 1:dim(combinations)[1]){
    cat("===============================================================\n")
    cat(paste("Kaplan-Meier Test between ",genes_of_interest[i,1]," and ",genes_of_interest[i,2],"\n",sep = ""))
    cat("===============================================================\n")
    surv_input = data.frame(cbind(all_data[,1:3],apply(all_data[,combinations[i,]],1,function(x) make_class(x))))
    colnames(surv_input)[4] = "Class"
    if(length(which(surv_input$Class==1))==0 || length(which(surv_input$Class==0))==0){
      print(paste("only one class for",paired_cnvs[ind],paired_survs[ind],i))
      next()
    }
    mfitA=survfit(Surv(survival_days, vital_status)  ~Class, data=surv_input)
    mdiff=survdiff(Surv(survival_days, vital_status)  ~Class, data=surv_input)
    print(mdiff)
    pdf(paste(plots_out,paste(paired_cnvs[ind],paired_survs[ind],genes_of_interest[i,1],genes_of_interest[i,2],"pdf",sep='.'),sep='/'),height=8.5,width=11)
    plot(mfitA, conf.int=FALSE,col=c("green", "red"),lty=1:1,lwd=2.5,xlab="Days",ylab="Survival %", cex.lab=1.5,cex=1.5,cex.main=1.5)
    title(paste(paired_survs[ind],paste("Chisq:",mdiff$chisq,sep=" "),sep=' ; '),cex=1.5)
    legend('topright', c(paste("RB1==-2 and ",genes_of_interest[i,2],"<0 (",mfitA$n[2],")",sep = ""),paste("RB1>-2 or ",genes_of_interest[i,2],">=0 (",mfitA$n[1],")",sep='')), lty=c(1), lwd=2.5, cex=1.5, col=c("red", "green"))
    dev.off()
  }
}
sink()

## same analysis for all other genes and RB1
plots_out = '/Users/holtj02/Desktop/Projects/Reva_Survivorship/TanyaHypothesis_Data/RB1_All_plots'
sink(paste(plots_out,"RB1_All_stats_out.txt",sep='/'))
for(ind in 1:length(paired_cnvs)){
  cat(paste("***********",paired_cnvs[ind],paired_survs[ind],"************\n",sep=' '))
  all_data = read.table(paste(out_dir,paste(paired_cnvs[ind],paired_survs[ind],sep='_'),sep='/'),sep='\t',header=T,check.names = F)
  genes_of_interest = cbind(rep("RB1",24776),colnames(all_data)[4:24779])
  combinations = cbind(rep(which(colnames(all_data)=="RB1"),24776),c(4:24779))
  for(i in 1:dim(combinations)[1]){
    cat("===============================================================\n")
    cat(paste("Kaplan-Meier Test between ",genes_of_interest[i,1]," and ",genes_of_interest[i,2],"\n",sep = ""))
    cat("===============================================================\n")
    surv_input = data.frame(cbind(all_data[,1:3],apply(all_data[,combinations[i,]],1,function(x) make_class(x))))
    colnames(surv_input)[4] = "Class"
    if(length(which(surv_input$Class==1))==0 || length(which(surv_input$Class==0))==0){
      print(paste("only one class for",paired_cnvs[ind],paired_survs[ind],i))
      next()
    }
    mfitA=survfit(Surv(survival_days, vital_status)  ~Class, data=surv_input)
    mdiff=survdiff(Surv(survival_days, vital_status)  ~Class, data=surv_input)
    print(mdiff)
    pdf(paste(plots_out,paste(paired_cnvs[ind],paired_survs[ind],genes_of_interest[i,1],genes_of_interest[i,2],"pdf",sep='.'),sep='/'),height=8.5,width=11)
    plot(mfitA, conf.int=FALSE,col=c("green", "red"),lty=1:1,lwd=2.5,xlab="Days",ylab="Survival %", cex.lab=1.5,cex=1.5,cex.main=1.5)
    title(paste(paired_survs[ind],paste("Chisq:",mdiff$chisq,sep=" "),sep=' ; '),cex=1.5)
    legend('topright', c(paste("RB1==-2 and ",genes_of_interest[i,2],"<0 (",mfitA$n[2],")",sep = ""),paste("RB1>-2 or ",genes_of_interest[i,2],">=0 (",mfitA$n[1],")",sep='')), lty=c(1), lwd=2.5, cex=1.5, col=c("red", "green"))
    dev.off()
  }
}
sink()

#surv_test = function(row,dat){
  surv_input = data.frame(cbind(dat[,1:3],apply(dat[,row],1,function(y) make_class(y))))
  colnames(surv_input)[4] = "Class"
  if(length(which(surv_input$Class==1))==0 || length(which(surv_input$Class==0))==0){
    return(paste("only one class"))
    next()
  }
  mfitA=survfit(Surv(survival_days, vital_status)  ~Class, data=surv_input)
  mdiff=survdiff(Surv(survival_days, vital_status)  ~Class, data=surv_input)
  return(mdiff)
}

# WGCNA Analysis
library(WGCNA)
start = ''
start_surv = ''
for(file in data_files){
  dat = read.table(paste('../r_cnvs/',file,sep=''),header=T,check.names = F)
  survData = as.data.frame(dat[,1:3])
  survData = as.data.frame(cbind(survData,rep(file,dim(survData)[1])))
  colnames(survData) = c(colnames(survData)[1:3],"cancer_type")
  rownames(dat) = dat$bcr_patient_barcode
  datExpr = as.data.frame(dat[,-c(1:3)])
  if(start == ''){
    start = datExpr
    start_surv = survData
  } else {
    start = as.data.frame(rbind(start,datExpr))
    start_surv = as.data.frame(rbind(start_surv,survData))
  }
  print(paste("loaded",file,sep=' '))
  print(paste("CNV data size",dim(start),sep=" "))
  print(paste("Survival data size",dim(start_surv),sep=" "))
}
datExprC = start
datTraits = start_surv[,2:dim(start_surv)[2]]
rownames(datTraits) = start_surv$bcr_patient_barcode
#saved the data sets at this point

# determine soft threshold power 
sft = pickSoftThreshold(t(datExprC),powers = c(1:10,rep(12,20,2)),verbose = 5)
sft$powerEstimate

# detect modules
bwnet = blockwiseModules(t(datExprC),maxBlockSize=2000,power=sft$powerEstimate,TOMType="signed",saveTOMS=TRUE,saveTOMFileBase="../r_data/bw_modules",networkType="signed",verbose=3)

# Module Characterization
survivor_data_list = list()
for(i in 0:18){
  print(paste("********* ",i," *********",sep=''))
  print(sum(datTraits[which(moduleLabels==i),1])/length(which(moduleLabels==1)))
  # the number of patients represented from each cancer for each module
  survivor_data_list[[i+1]] = list(table(datTraits[which(moduleLabels==i),3]))
  # the number of patients from each cancer within each module divided by the total number of patients for each cancer (respective to each cancer)
  survivor_data_list[[i+1]][[2]] = survivor_data_list[[i+1]][[1]]/table(datTraits$cancer_type)
  # which cancer has the largest percentage of its patients represented in the module
  survivor_data_list[[i+1]][[3]] = names(which(survivor_data_list[[i+1]][[2]]==max(survivor_data_list[[i+1]][[2]])))
  # the rate of individuals living within the module
  survivor_data_list[[i+1]][[4]] = 1-(sum(datTraits[which(moduleLabels==i),1])/length(which(moduleLabels==i)))
  # living rate within the module minus the living rate for all data 
  survivor_data_list[[i+1]][[5]] = survivor_data_list[[i+1]][[4]] - (1-(sum(datTraits$vital_status)/dim(datTraits)[1]))
  # 
  survivor_data_list[[i+1]][[6]] = sapply(names(survivor_data_list[[i+1]][[1]]),function(x) (1-(sum(datTraits[which(datTraits[which(moduleLabels==i),3]==x),1])/length(datTraits[which(datTraits[which(moduleLabels==i),3]==x),1]))))
  # Rate of individuals living for each cancer
  survivor_data_list[[i+1]][[7]] = sapply(names(survivor_data_list[[i+1]][[1]]),function(x) (1-(sum(datTraits[which(datTraits[,3]==x),1])/length(datTraits[which(datTraits[,3]==x),1]))))
  # 
  survivor_data_list[[i+1]][[8]] = survivor_data_list[[i+1]][[6]] - survivor_data_list[[i+1]][[7]]
  names(survivor_data_list[[i+1]]) = c("Tumors_Represented","Norm_By_N_Patients_Representation","Best_Described_Cancer","Living_Rate","Distictive_Living_Rate","Per_Cancer_Module_Living_Rate","Per_Cancer_Global_Living_Rate","Per_Cancer_Distictive_Living_Rate")
}

# which genes have greatest difference between the set of living and dead patients.
dead = dat[which(datTraits$vital_status == 1),]
alive = dat[-which(datTraits$vital_status == 1),]
gene_differences = colSums(dead)/dim(dead)[1] - colSums(alive)/dim(alive)[1]
gene_differences[order(abs(gene_differences),decreasing = TRUE)]

# what are the best predictors for module membership (i.e. what genes play the most significant role in
# the differences between two modules)
# test this on all combinations of modules.
datExprC_labeled = as.data.frame(cbind(moduleLabels,datExprC))
colnames(datExprC_labeled)[1] = "mod_labels"
combinations = combn(unique(moduleLabels),2)
data_files = list.files('r_cnvs')
for(i in 1:dim(combinations)[2]){
  two_mod_comp_data = datExprC_labeled[which(!is.na(match(moduleLabels,combinations[,i]))),]
  two_mod_comp_traits = datTraits[which(!is.na(match(moduleLabels,combinations[,i]))),]
  for(j in 1:15){
    #contains all patients in both modules for a particular cancer
    canc_two_mod_comp_data = two_mod_comp_data[which(two_mod_comp_traits[,3]==data_files[j]),]
    canc_two_mod_comp_traits = two_mod_comp_traits[which(two_mod_comp_traits[,3]==data_files[j]),]
    # subset 
    canc_one_data = canc_two_mod_comp_data[which(canc_two_mod_comp_data[,1]==combinations[1,i]),]
    canc_two_data = canc_two_mod_comp_data[which(canc_two_mod_comp_data[,1]==combinations[2,i]),]
    canc_one_traits = canc_two_mod_comp_traits[which(canc_two_mod_comp_data[,1]==combinations[1,i]),]
    canc_two_traits = canc_two_mod_comp_traits[which(canc_two_mod_comp_data[,1]==combinations[2,i]),]
    # inner module difference between dead and living genes
    max(abs((colSums(canc_one_data[which(canc_one_traits[,1]==1),2:dim(canc_one_data)[2]])/dim(canc_one_data[which(canc_one_traits[,1]==1),])[1] - colSums(canc_one_data[which(canc_one_traits[,1]==0),2:dim(canc_one_data)[2]])/dim(canc_one_data[which(canc_one_traits[,1]==0),])[1]) - (colSums(canc_two_data[which(canc_two_traits[,1]==1),2:dim(canc_two_data)[2]])/dim(canc_two_data[which(canc_two_traits[,1]==1),])[1] - colSums(canc_two_data[which(canc_two_traits[,1]==0),2:dim(canc_two_data)[2]])/dim(canc_two_data[which(canc_two_traits[,1]==0),])[1])))
  }
  #mod_model = glm(mod_labels ~ .,data = two_mod_comp_data,family = binomial)
}

# considering the full data for a particular cancer, rank genes in decreasing order that have the greatest differences between living and dead patients
canc_list = list()
i=1
canc_data = datExprC[which(datTraits[,3]==data_files[i]),]
canc_traits = datTraits[which(datTraits[,3]==data_files[i]),]
canc_list[[1]] = as.data.frame(cbind(order(abs(colSums(canc_data[which(canc_traits[,1]==1),])/length(which(canc_traits[,1]==1)) - colSums(canc_data[which(canc_traits[,1]==0),])/length(which(canc_traits[,1]==0))),decreasing = TRUE)))
colnames(canc_list[[1]])[i] = strsplit(data_files[i],split = '.',fixed = T)[[1]][1]
canc_list[[2]] = as.data.frame(cbind(colnames(canc_data)[order(abs(colSums(canc_data[which(canc_traits[,1]==1),])/length(which(canc_traits[,1]==1)) - colSums(canc_data[which(canc_traits[,1]==0),])/length(which(canc_traits[,1]==0))),decreasing = TRUE)]))
colnames(canc_list[[2]])[i] = strsplit(data_files[i],split = '.',fixed = T)[[1]][1]
for(i in 2:length(data_files)){
  cat(paste(data_files[i],"\n",sep=""))
  canc_data = datExprC[which(datTraits[,3]==data_files[i]),]
  canc_traits = datTraits[which(datTraits[,3]==data_files[i]),]
  canc_list[[1]] = as.data.frame(cbind(canc_list[[1]],order(abs(colSums(canc_data[which(canc_traits[,1]==1),])/length(which(canc_traits[,1]==1)) - colSums(canc_data[which(canc_traits[,1]==0),])/length(which(canc_traits[,1]==0))),decreasing = TRUE)))
  colnames(canc_list[[1]])[i] = strsplit(data_files[i],split = '.',fixed = T)[[1]][1]
  canc_list[[2]] = as.data.frame(cbind(canc_list[[2]],colnames(canc_data)[order(abs(colSums(canc_data[which(canc_traits[,1]==1),])/length(which(canc_traits[,1]==1)) - colSums(canc_data[which(canc_traits[,1]==0),])/length(which(canc_traits[,1]==0))),decreasing = TRUE)]))
  colnames(canc_list[[2]])[i] = strsplit(data_files[i],split = '.',fixed = T)[[1]][1]
}

# repeat the same analysis but compare genes by correlating the frequencies of the 5 cnv values (-2,-1,0,1,2)
apply(canc_data,2,function(l) cor_per_gene(l,canc_traits,base))
cor_per_gene = function(column,trait_data,base){
  tab_1 = table(column[which(trait_data[,1]==1)])
  tab_2 = table(column[which(trait_data[,1]==0)])
  if(length(tab_1) != 5){
    for(j in which(is.na(match(base,names(tab_1))))){
      tab_1 = append(tab_1,0,j-1)
    }
  }
  if(length(tab_2) != 5){
    for(j in which(is.na(match(base,names(tab_2))))){
      tab_2 = append(tab_2,0,j-1)
    }
  }
  return(1-cor(x=tab_1,y=tab_2))
}

canc_list_cor = list()
i=1
canc_data = datExprC[which(datTraits[,3]==data_files[i]),]
canc_traits = datTraits[which(datTraits[,3]==data_files[i]),]
canc_list_cor[[1]] = as.data.frame(cbind(order(apply(canc_data,2,function(l) cor_per_gene(l,canc_traits,base)),decreasing = T)))
colnames(canc_list_cor[[1]])[i] = strsplit(data_files[i],split = '.',fixed = T)[[1]][1]
canc_list_cor[[2]] = as.data.frame(cbind(colnames(canc_data)[canc_list_cor[[1]][,i]]))
colnames(canc_list_cor[[2]])[i] = strsplit(data_files[i],split = '.',fixed = T)[[1]][1]
for(i in 2:length(data_files)){
  cat(paste(data_files[i],"\n",sep=""))
  canc_data = datExprC[which(datTraits[,3]==data_files[i]),]
  canc_traits = datTraits[which(datTraits[,3]==data_files[i]),]
  canc_list_cor[[1]] = as.data.frame(cbind(canc_list_cor[[1]],order(apply(canc_data,2,function(l) cor_per_gene(l,canc_traits,base)),decreasing = T)))
  colnames(canc_list_cor[[1]])[i] = strsplit(data_files[i],split = '.',fixed = T)[[1]][1]
  canc_list_cor[[2]] = as.data.frame(cbind(canc_list_cor[[2]],colnames(canc_data)[canc_list_cor[[1]][,i]]))
  colnames(canc_list_cor[[2]])[i] = strsplit(data_files[i],split = '.',fixed = T)[[1]][1]
}


"""
406 blca.cnv_clinical_blca.new
1080 brca.cnv_clinical_brca.new
612 coadread.cnv_clinical_coadread.new
185 esca.cnv_clinical_esca.new
529 kirc.cnv_clinical_kirc.new
288 kirp.cnv_clinical_kirp.new
513 lgg.cnv_clinical_lgg.new
370 lihc.cnv_clinical_lihc.new
565 ov.cnv_clinical_ov.new
185 paad.cnv_clinical_paad.new
493 prad.cnv_clinical_prad.new
258 sarc.cnv_clinical_sarc.new
359 skcm.cnv_clinical_skcm.new
435 stad.cnv_clinical_stad.new
539 ucec.cnv_clinical_ucec.new
"""

## Fisher Exact Test
PAIRS = list(c(-2,-2),list(c(-2,-1),c(-2,-1)),list(c(-1,-2),c(-1,-2)),list(c(-1,-1),c(-1,-1)))
create_contingency = function(surv_data,pair,lt=FALSE){
  # this function returns a n x 4 matrix of values for a fisher contingency table.
  # n is calculated from the number of pairs of CNV values one wants to generate contingency table info for.
  n = dim(surv_data)[1]
  if (lt) {
    lab = "lt"
    nd1d2 = length(which((surv_data[,1] <= pair[1]) & (surv_data[,2] <= pair[2])))
    nd1 = length(which(surv_data[,1] <= pair[1]))
    nd2 = length(which(surv_data[,2] <= pair[2]))
  } else {
    lab = "eq"
    nd1d2 = length(which((surv_data[,1] == pair[1]) & (surv_data[,2] == pair[2])))
    nd1 = length(which(surv_data[,1] == pair[1]))
    nd2 = length(which(surv_data[,2] == pair[2]))
  }
  return(c(lab,pair,nd1d2,(nd1-nd1d2),(nd2-nd1d2),n-((nd1+nd2)-nd1d2)))
}

contingency_input = as.data.frame(t(matrix(unlist(lapply(PAIRS,function(x) if(class(x)=="list"){ cbind(create_contingency(all_data[,combinations[i,]],x[[1]]),create_contingency(all_data[,combinations[i,]],x[[2]],lt=TRUE)) }else{ create_contingency(all_data[,combinations[i,]],x) })),nrow = 7)))

rb_dict = function(nt){
  if(nt == 'A' | nt == 'C'){
    return('R')
  } else {
    return('B')
  }
}

rb_conv = function(bc){
  return(paste(sapply(strsplit(bc,split = '')[[1]],function(x) rb_dict(x)),collapse = ''))
}

solve_hgd = function(d,e,f,g){
  N = d+e+f+g
  r = d+f
  n = f+g
  max_k = min(c(r,n))
  min_k = max(0,r+n-N)
  cutoff = hgd(f,r,n,N)
  print(paste("cutoff:",cutoff))
  tmp_p = 0.0
  for(k in min_k:(max_k)){
    doub_p = hgd(k,r,n,N)
    #tmp_p = tmp_p + doub_p
    print(paste("doub_p:",doub_p,"tmp_p:",tmp_p,"1-tmp_p:",1-tmp_p,"k:",k))
    if(doub_p<=cutoff){
      tmp_p = tmp_p+doub_p
    } else {
      return(tmp_p)
    }
    #print(paste("doub_p:",doub_p,"tmp_p:",tmp_p,"1-tmp_p:",1-tmp_p,"k:",k))
  }
  return(tmp_p)
}

FUN <- function(x) {
  x <- as.integer(x)
  div <- seq_len(abs(x))
  factors <- div[x %% div == 0L]
  factors <- list(neg = -factors, pos = factors)
  return(factors)
}
library(gmp)
prod(sapply(seq(0,N-1,1),function(x) ))

library(parallel)
library(doParallel)
library(foreach)
the_mat = matrix(c(10,20,30,40),byrow = TRUE,nrow = 2)
foreach(x=list(c(1,1),c(1,2),c(2,1),c(2,2))) %dopar% {
  start_time = Sys.time()
  cat("Start_time:",start_time,"\n",file = paste0("debug_file_",Sys.getpid(),".txt"),append = TRUE)
  for(i in 1:10000){
    cat(i,the_mat[x[1],x[2]],Sys.getpid(),Sys.time(),"\n",file = paste0("debug_file_",Sys.getpid(),".txt"),append = TRUE)
  }
  end_time = Sys.time()
  cat("End_time:",end_time,"\n",file = paste0("debug_file_",Sys.getpid(),".txt"),append = TRUE)
  cat("Time_difference:",end_time-start_time,"\n",file = paste0("debug_file_",Sys.getpid(),".txt"),append = TRUE)
}
