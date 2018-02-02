options(stringsAsFactors = F)
args <- commandArgs(trailingOnly = TRUE)
genes_file = args[1]
canc_data = args[2]
#filter_file = args[2]
workdir = '/sc/orga/work/holtj02/Other_Projects/Reva_SurvAnalysis/TanyaHypothesis_Data'
setwd(workdir)
plots_out = paste('/sc/orga/work/holtj02/Other_Projects/Reva_SurvAnalysis/TanyaHypothesis_Data/Filtered_Plots',canc_data,"survivorship_plots",sep='/')
dir.create(plots_out)

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

#### Survival Analysis
library(survival)

## get list of all cnv data files
out_dir = 'r_cnvs'

## read in the genes file and filter file
print("Loading genes file...")
genes_of_interest = read.csv(genes_file,header=T)

print("Loading cnv data...")
sink(paste(plots_out,"All_stats_out.txt",sep='/'))
cat(paste("***********",canc_data,"************\n",sep=' '))
all_data = read.table(paste(out_dir,canc_data,sep='/'),sep='\t',header=T,check.names = F)
combinations = as.data.frame(cbind(match(genes_of_interest[,2],colnames(all_data)),match(genes_of_interest[,3],colnames(all_data))))
for(j in 1:dim(genes_of_interest)[2]){
    if(j < 11){
        next()
    }
    genes_to_test = which(genes_of_interest[,j]<0.001)
    #print(length(genes_to_test))
    if(length(genes_to_test)==0){
        next()
    }
    for(i in genes_to_test){
        test = colnames(genes_of_interest)[j]
        cat("===============================================================\n")
        cat(paste("Kaplan-Meier Test between ",genes_of_interest[i,2]," and ",genes_of_interest[i,3]," for test ",test,"\n",sep = ""))
        cat("===============================================================\n")
        #print(combinations[i,])
        #print(dim(all_data[,c(combinations[i,1],combinations[i,2])]))
        #print(paste("combinations",combinations[i,]))
        #print(test)
        surv_input = data.frame(cbind(all_data[,1:3],apply(all_data[,c(combinations[i,1],combinations[i,2])],1,function(x) make_class(x,test))))
        colnames(surv_input)[4] = "Class"
        if(length(which(surv_input$Class==1))==0 || length(which(surv_input$Class==0))==0){
            print(paste("only one class for",canc_data,genes_of_interest[i,2],genes_of_interest[i,3]))
            next()
        }
        #print(as.data.frame(cbind(surv_input,all_data[,c(combinations[i,1],combinations[i,2])])))
        mfitA=survfit(Surv(survival_days, vital_status)  ~Class, data=surv_input)
        mdiff=survdiff(Surv(survival_days, vital_status)  ~Class, data=surv_input)
        print(mdiff)
        pdf(paste(plots_out,paste(canc_data,genes_of_interest[i,2],genes_of_interest[i,3],test,"pdf",sep='.'),sep='/'),height=8.5,width=11)
        plot(mfitA, conf.int=FALSE,col=c("green", "red"),lty=1:1,lwd=2.5,xlab="Days",ylab="Survival %", cex.lab=1.5,cex=1.5,cex.main=1.5)
        title(paste(canc_data,test,paste(genes_of_interest[i,2],"and",genes_of_interest[i,3]),paste("Chisq:",mdiff$chisq,sep=" "),sep=' ; '),cex=1.5)
        #if(test=="Pval_D1AD1B"){
        #    legend('topright', c(paste(genes_of_interest[i,2],"=-1 (",mfitA$n[2],")",sep = ""),paste(genes_of_interest[i,3],"=-1 (",mfitA$n[1],")",sep='')), lty=c(1),lwd=2.5,cex=1.5,col=c("red","green"))
        #} else if(test=="Pval_D1AD2B"){
        #    legend('topright', c(paste(genes_of_interest[i,2],"=-1 (",mfitA$n[2],")",sep = ""),paste(genes_of_interest[i,3],"=-2 (",mfitA$n[1],")",sep='')), lty=c(1),lwd=2.5,cex=1.5,col=c("red","green"))
        #} else if(test=="Pval_D1AD3B"){
        #    legend('topright', c(paste(genes_of_interest[i,2],"=-1 (",mfitA$n[2],")",sep = ""),paste(genes_of_interest[i,3],"<0 (",mfitA$n[1],")",sep='')), lty=c(1),lwd=2.5,cex=1.5,col=c("red","green"))
        #} else if(test=="Pval_D2AD1B"){
        #    legend('topright', c(paste(genes_of_interest[i,2],"=-2 (",mfitA$n[2],")",sep = ""),paste(genes_of_interest[i,3],"=-1 (",mfitA$n[1],")",sep='')), lty=c(1),lwd=2.5,cex=1.5,col=c("red","green"))
        #} else if(test=="Pval_D3AD1B"){
        #    legend('topright', c(paste(genes_of_interest[i,2],"<0 (",mfitA$n[2],")",sep = ""),paste(genes_of_interest[i,3],"=-1 (",mfitA$n[1],")",sep='')), lty=c(1),lwd=2.5,cex=1.5,col=c("red","green"))
        #} else if(test=="Pval_D2AD2B"){
        #    legend('topright', c(paste(genes_of_interest[i,2],"=-2 (",mfitA$n[2],")",sep = ""),paste(genes_of_interest[i,3],"=-2 (",mfitA$n[1],")",sep='')), lty=c(1),lwd=2.5,cex=1.5,col=c("red","green"))
        #} else if(test=="Pval_D2AD3B"){
        #    legend('topright', c(paste(genes_of_interest[i,2],"=-2 (",mfitA$n[2],")",sep = ""),paste(genes_of_interest[i,3],"<0 (",mfitA$n[1],")",sep='')), lty=c(1),lwd=2.5,cex=1.5,col=c("red","green"))
        #} else if(test=="Pval_D3AD2B"){
        #    legend('topright', c(paste(genes_of_interest[i,2],"<0 (",mfitA$n[2],")",sep = ""),paste(genes_of_interest[i,3],"=-2 (",mfitA$n[1],")",sep='')), lty=c(1),lwd=2.5,cex=1.5,col=c("red","green"))
        #} else if(test=="Pval_D3AD3B"){
        #    legend('topright', c(paste(genes_of_interest[i,2],"<0 (",mfitA$n[2],")",sep = ""),paste(genes_of_interest[i,3],"<0 (",mfitA$n[1],")",sep='')), lty=c(1),lwd=2.5,cex=1.5,col=c("red","green"))
        #}
        legend('topright', c(paste("Both Depleted (",mfitA$n[2],")",sep = ""),paste("XOR Depleted (",mfitA$n[1],")",sep='')), lty=c(1),lwd=2.5,cex=1.5,col=c("red","green"))
        dev.off()
    }
}
sink()

