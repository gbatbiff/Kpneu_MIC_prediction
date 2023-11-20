### Rscript pop_struct_enet.R pheno_file matrix_file st_file
rm(list = ls())
library(OneR)
library(purrr)
library(scales)
library(ranger)
library(glmnet)
library(caret)
library(doParallel)
#library(doMC)
#registerDoMC(cores = 40)
registerDoParallel(cores = 30)

args=commandArgs(trailingOnly = TRUE)

phen_file <- args[1]
matrix_file<-args[2]
simulation_type<-args[4]
dd<-read.delim(matrix_file,sep='\t')
#dd<-read.delim("../pres_abs_matrix.txt",sep='\t')
#colnames(dd)<-gsub("X","",colnames(dd))
#colnames(dd)<-gsub(".fna","",colnames(dd))
rownames(dd)<-dd$ID
dd$ID<-NULL
###
tdd<-as.data.frame(t(dd))


n_bins<-c(4,6,8,10)
set.seed(1234)

for (n in n_bins){


  
  phen<-read.delim(phen_file,sep='',header=F)
  #phen<-read.delim("simul_1000SNPs_EF_2.5_herit_0.9_traitprev_0.1.phen",sep='',header=F)
  phen$V1<-NULL
  phen$V2<-gsub("X","",phen$V2)
  phen$V2<-gsub(".fna","",phen$V2)
  colnames(phen)<-c("sample","MIC")
  
  phen$sim_ord_mic<-rescale(phen$MIC,to=c(0.25,16))
  phen$sim_ord_mic<-bin(phen$sim_ord_mic, nbins = n,method = "length") 
  phen$MIC<-NULL
  phen$sim_ord_mic<-gsub("[(]","",phen$sim_ord_mic)
  phen$sim_ord_mic<-gsub("[]]","",phen$sim_ord_mic)
  col1<-as.data.frame(unlist(map(strsplit(phen$sim_ord_mic, split = ","), 1)))
  colnames(col1)[1]<-"elem1"
  col1$elem1<-as.numeric(as.character(col1$elem1))
  col2<-as.data.frame(unlist(map(strsplit(phen$sim_ord_mic, split = ","), 2)))
  colnames(col2)[1]<-"elem2"
  col2$elem2<-as.numeric(as.character(col2$elem2))
  midpoint_df<-cbind(col1,col2)
  midpoint_df$midpoint<-(midpoint_df$elem1+midpoint_df$elem2)/2
  midpoint_df$elem1<-NULL
  midpoint_df$elem2<-NULL
  phen_mid<-as.data.frame(cbind(as.character(phen$sample),midpoint_df$midpoint))
  colnames(phen_mid)<-c("sample","midpoint_mic")
  
  #sorted_bins<-sort(as.numeric(as.character(levels(phen_mid$midpoint_mic))))
  sorted_bins<-sort(unique(as.numeric(phen_mid$midpoint_mic)))
  
  phen_mid$midpoint_mic<-as.numeric(as.character(phen_mid$midpoint_mic))
  
  #last<-tail(sort(table(phen_mid$midpoint_mic),decreasing=T),3)
  #last<-as.numeric(rownames(last))
  #censored_class<-sorted_bins[4] ### ULTIMA CLASSE BILANCIATA
  #phen_mid$midpoint_mic<-ifelse(phen_mid$midpoint_mic %in% last,censored_class,phen_mid$midpoint_mic) ### SE CLASSI APPARTENGONO A SBILANCIATE INCORPORA IN ULTIMA (censor)
  
  colnames(phen_mid)<-c("ID","MIC")
  
  tdd$ID<-rownames(tdd)
  tdd$ID<-gsub("_out","",tdd$ID)
  tdd$ID<-gsub(".out","",tdd$ID)
  
  phen_mid$ID<-gsub("_out","",phen_mid$ID)
  phen_mid$ID<-gsub(".out","",phen_mid$ID)
  
 


  ###
  mm<-merge(tdd,phen_mid,by="ID")
  rownames(mm)<-mm$ID
  ###
  st_file<-args[3]
  st<-read.delim(st_file,header = F)
  st$V2<-gsub("-","_",st$V2)
  st$V1<-as.character(st$V1)
  st$V2<-as.character(st$V2)
  colnames(st)<-c("ID","PP")
  ###
  mm<-merge(mm,st,by="ID")
  ids<-mm$ID
  mm1<-mm[,c("group_20046","ubiE","rpmJ","PP","MIC")] ###dummy df to avoid mem limit

  pred<-dummyVars(formula = "~.",data = mm1)
  df<-data.frame(predict(pred,newdata = mm1),row.names=ids)
  colnames(df)<-gsub("PPPP","PP",colnames(df))

  pp_lineages_df<-df[,4:ncol(df)]

  df<-cbind(mm,pp_lineages_df)
  ###

  ###

  
  
  
  ###
  
  ###
  trainIndex <- createDataPartition(df$MIC,  p = .70, list = FALSE, times=1) ###partizione 70%train-30%test
  
  train_y<-df[trainIndex,]$MIC
  test_y<-df[-trainIndex,]$MIC
  
  ###
  
  train_df<-df[trainIndex,]
  test_df<-df[-trainIndex,]
  test_df$MIC<-NULL
  
 



  #1/(size of cluster) for every sample in the cluster
  st_frst_pos<-min(grep("PP",colnames(df)))
  st_last_pos<-max(grep("PP",colnames(df)))
  list_st<-grep("PP",colnames(df),value = T)
  ###
  
  ### dataframe con frequenze pesate
  freqs<-data.frame()
  for (i in list_st[-1]) {
    stdf<-as.data.frame(table(df[,i]))
    stdf$type<-i
    stdf<-stdf[which(stdf$Var1==1),]
    stdf$Var1<-NULL
    freqs[i,1]<-stdf$Freq
  }
  
  freqs$reweight<-1/freqs$V1
  ### mettere nella variabile weight 1 per tutti gli snp e per i st i valori del vettore reweight
  st_frst_pos<-min(grep("PP",colnames(df)))
  freqs$V1<-rownames(freqs$V1)
  freqs$PP<-rownames(freqs)
  ###
  obs<-rownames(train_df) ### ceppi da pesare nel training set
  ###usare il vettore dei ceppi ed associare il reweight calcolato nella tabella mergiata e dare questo weight a cv.glmnet
  st_weight_merged<-merge(freqs,st,by="PP") 
  ###
  rownames(st_weight_merged)<-st_weight_merged$ID
  pesi<-st_weight_merged[obs,]
  reweight<-pesi$reweight
  
  ### FUNZIONE PER TESTARE DI QUANTO SI SCOSTA UNA CLASSE PREDETTA, ERROR=0 OTTIMO, ERROR=1 SI SCOSTA DI UN BIN
  loosely_correct <- function(correct, predicted, bins, error=0) {
    bins <- sort(bins)
    correct_bin <- findInterval(correct, bins)
    predicted_bin <- findInterval(predicted, bins)
    if(abs(correct_bin - predicted_bin) <= error) {
      correct <- TRUE
    } else {
      correct <- FALSE
    }
    return(correct)
  }
  
  r2_df_train<-data.frame()
  r2_df_test<-data.frame()
  
    start_time_runtime <- Sys.time()

  	df_stat<-data.frame()
  	
	rf_perm_regr <- ranger(y=train_df$MIC, x=train_df[colnames(train_df)!="MIC"], num.trees = 500, importance="impurity", write.forest = TRUE,case.weights=reweight,num.threads=50)
	


	r2_train<-rf_perm_regr$r.squared


	start_time_prediction <- Sys.time()


	pred_rf<-predict(rf_perm_regr,test_df,type = "response")
	
end_time_prediction <-	Sys.time()
  
  prediction_time <- end_time_prediction - start_time_prediction

  prediction_time<-as.character(prediction_time)
  
  end_time_runtime <- Sys.time()
  
  runtime_analysis<- end_time_runtime - start_time_runtime
  runtime_analysis<-as.character(runtime_analysis)
  
	sst <- sum((test_y - mean(test_y))^2)
	sse <- sum((pred_rf$predictions - test_y)^2)
	r2_test <- 1 - sse / sst


	
laststr<-gsub("^.*/", "", phen_file)
h2<-unlist(strsplit(laststr,split = "_"))[6]


df_stat[paste(n,"bins",sep="_"),c("R2_train","R2_test","simulation","model","type","h2","prediction_time","train_test_time")]<-c(r2_train,r2_test,simulation_type,"RF","midpoint_regression",h2,prediction_time,runtime_analysis)
    

  
  write.table(df_stat,sprintf("%s_RF_R2_midpoint_regression_%s_bins.txt",phen_file,n),
            sep = '\t',quote = F)



  check_imp<-as.data.frame(importance(rf_perm_regr))


  write.table(check_imp,sprintf("%s_RF_R2_imp_pred_midpoint_regression_%s_bins.txt",phen_file,n),
            sep = '\t',quote = F,col.names=NA)




    
}
  
