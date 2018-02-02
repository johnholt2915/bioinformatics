options(stringsAsFactors = F)
args <- commandArgs(trailingOnly = TRUE)
cnv_file = args[1]
genes_file = args[2]
workdir = '/sc/orga/work/holtj02/Other_Projects/Reva_SurvAnalysis/TanyaHypothesis_Data'
setwd(workdir)
plots_out = paste('/sc/orga/work/holtj02/Other_Projects/Reva_SurvAnalysis/TanyaHypothesis_Data/RB1_All_plots',cnv_file,sep='/')
dir.create(plots_out)
out_dir = 'r_cnvs'

make_class = function(row) {
  if(row[1] == -2 && row[2] < 0){
    return(1)
  } else {
    return(0)
  }
}

#### Survival Analysis
library(survival)

sink(paste(plots_out,paste(cnv_file,"RB1_All_stats_out.txt",sep='-'),sep='/'))
#print(paste(plots_out,paste(cnv_file,"RB1_All_stats_out.txt",sep='-'),sep='/'))
cat(paste("***********",cnv_file,"************\n",sep=' '))
all_data = read.table(paste(out_dir,cnv_file,sep='/'),sep='\t',header=T,check.names = F)
genes_of_interest = read.csv(genes_file,header=F)
combinations = apply(genes_of_interest,c(1,2),function(x) which(colnames(all_data)==x))
#genes_of_interest = cbind(rep("RB1",24776),colnames(all_data)[4:24779])
#combinations = cbind(rep(which(colnames(all_data)=="RB1"),24776),c(4:24779))
for(i in 1:dim(combinations)[1]){
    cat("===============================================================\n")
    cat(paste("Kaplan-Meier Test between ",genes_of_interest[i,1]," and ",genes_of_interest[i,2],"\n",sep = ""))
    cat("===============================================================\n")
    surv_input = data.frame(cbind(all_data[,1:3],apply(all_data[,combinations[i,]],1,function(x) make_class(x))))
    colnames(surv_input)[4] = "Class"
    if(length(which(surv_input$Class==1))==0 || length(which(surv_input$Class==0))==0){
        print(paste("only one class for",cnv_file,genes_of_interest[i,1],genes_of_interest[i,2]))
        next()
    }
    mfitA=survfit(Surv(survival_days, vital_status)  ~Class, data=surv_input)
    mdiff=survdiff(Surv(survival_days, vital_status)  ~Class, data=surv_input)
    print(mdiff)
    pdf(paste(plots_out,paste(cnv_file,genes_of_interest[i,1],genes_of_interest[i,2],"pdf",sep='.'),sep='/'),height=8.5,width=11)
    plot(mfitA, conf.int=FALSE,col=c("green", "red"),lty=1:1,lwd=2.5,xlab="Days",ylab="Survival %", cex.lab=1.5,cex=1.5,cex.main=1.5)
    title(paste(cnv_file,paste("Chisq:",mdiff$chisq,sep=" "),sep=' ; '),cex=1.5)
    legend('topright', c(paste("RB1==-2 and ",genes_of_interest[i,2],"<0 (",mfitA$n[2],")",sep = ""),paste("RB1>-2 or ",genes_of_interest[i,2],">=0 (",mfitA$n[1],")",sep='')), lty=c(1), lwd=2.5, cex=1.5, col=c("red", "green"))
    dev.off()
}
sink()

