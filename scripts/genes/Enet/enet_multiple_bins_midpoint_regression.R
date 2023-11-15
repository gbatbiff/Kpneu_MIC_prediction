### Rscript pop_struct_enet.R pheno_file matrix_file st_file

rm(list = ls())
library(OneR)
library(purrr)
library(scales)
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

rownames(dd)<-dd$ID
dd$ID<-NULL
###
tdd<-as.data.frame(t(dd))


n_bins<-c(4,6,8,10)

set.seed(1234)

for (n in n_bins){
  
  phen<-read.delim(phen_file,sep='',header=F)
  #phen<-read.delim("simul_3betalact_EF_2.5_herit_0.1_traitprev_0.2.phen",sep='',header=F)
  phen$V1<-NULL
  phen$V2<-gsub("X","",phen$V2)
  phen$V2<-gsub(".fna","",phen$V2)
  phen$V2<-gsub("_out","",phen$V2)
  phen$V2<-gsub(".fna.out","",phen$V2)
  phen$V2<-gsub(".out","",phen$V2)
  
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
    
  colnames(phen_mid)<-c("ID","MIC")
  
  tdd$ID<-rownames(tdd)
  ###
  colnames(tdd)<-gsub("~~~","_",colnames(tdd))
  
  mm<-merge(tdd,phen_mid,by="ID")
  rownames(mm)<-mm$ID
  
  ###
  st_file<-args[3]
  st<-read.delim(st_file,header = F)
  #st<-read.delim("renamed_lineages.txt",header = F)
  st$V1<-as.character(st$V1)
  st$V2<-as.character(st$V2)
  colnames(st)<-c("ID","PP")
  ###
  
  mm<-merge(mm,st,by="ID")
  ids<-mm$ID
 
  
  mm1<-mm[,c("group_2462","group_1619","group_1272" ,"PP","MIC")] ###dummy df to avoid mem limit
  
  pred<-dummyVars(formula = "~.",data = mm1)
  df<-data.frame(predict(pred,newdata = mm1),row.names=ids)
  colnames(df)<-gsub("PPPP","PP",colnames(df))
  pp_lineages_df<-df[,4:ncol(df)]
  df<-cbind(mm,pp_lineages_df)
  df$ID<-NULL
  df$PP<-NULL
  ###
  ###
  

trainIndex <- createDataPartition(df$MIC,  p = .70, list = FALSE, times=1) ###partizione 70%train-30%test



train_y<-df[trainIndex,]$MIC
test_y<-df[-trainIndex,]$MIC
df$MIC<-NULL
###
df$ID<-NULL
#df$PP<-NULL

train_df<-data.matrix(df[trainIndex,])
test_df<-data.matrix(df[-trainIndex,])
#1/(size of cluster) for every sample in the cluster
st_frst_pos<-min(grep("PP",colnames(df)))
st_last_pos<-max(grep("PP",colnames(df)))
list_st<-grep("PP",colnames(df),value = T)
###

#train_df<-train_df[,500:750]


### dataframe con frequenze pesateù


freqs<-data.frame()
for (i in list_st) {
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

  
  


  r2_df_train<-data.frame()
  r2_df_test<-data.frame()
  
  
  df_stat<-data.frame()
  
  
  df_pred<-data.frame()
  
  
  
  start_time_runtime <- Sys.time()
  
  
  cv_elnet<- cv.glmnet(x=train_df, y=train_y, alpha=0.01,
                       type.measure="mse", family="gaussian", keep = TRUE,parallel = TRUE,weights = reweight)
  
  start_time_prediction <- Sys.time()

  pred_enet<-predict(cv_elnet,s="lambda.1se",newx = test_df) ### PREDICTION
  
  end_time_prediction <-	Sys.time()
  
  prediction_time <- end_time_prediction - start_time_prediction

  prediction_time<-as.character(prediction_time)
  
  end_time_runtime <- Sys.time()
  
  runtime_analysis<- end_time_runtime - start_time_runtime
  runtime_analysis<-as.character(runtime_analysis)
  
  
  r2_train<-cv_elnet$glmnet.fit$dev.ratio[which(cv_elnet$glmnet.fit$lambda == cv_elnet$lambda.min)] ### R2 train set
  r2_df_train[paste(n,"bins",sep="_"),1]<-r2_train ###DF CREATION WITH ALPHA AND R2
  colnames(r2_df_train)<-"R2_train" 
  
  
  sst <- sum((test_y - mean(test_y))^2)
  sse <- sum((pred_enet - test_y)^2)
  r2_test <- 1 - sse / sst
  r2_df_test[paste(n,"bins",sep="_"),1]<-r2_test
  colnames(r2_df_test)<-"R2_test" 
  
  
  laststr<-gsub("^.*/", "", phen_file)
  h2<-unlist(strsplit(laststr,split = "_"))[6]
  
  df_stat[paste(n,"bins",sep="_"),c("R2_train","R2_test","simulation","model","type","h2","prediction_time","train_test_time")]<-c(r2_train,r2_test,simulation_type,"enet","midpoint_regression",h2,prediction_time,runtime_analysis)
  
  
  
  write.table(df_stat,sprintf("%s_R2_enet_midpoint_regression_%s_bins.txt",phen_file,n),
              sep = '\t',quote = F,col.names=NA)
  
  
  coef_elnet<-coef(cv_elnet,s="lambda.1se")
  coef_elnet<-as.data.frame(data.matrix(coef_elnet))
  
  # coef_df_summary<-as.data.frame(as.matrix(coef_elnet))
    colnames(coef_elnet)<-"coef_importance"
    
  
  
  write.table(coef_elnet,sprintf("%s_summary_coef_R2_enet_midpoint_regression_%s_bins.txt",phen_file,n),
              sep = '\t',quote = F,col.names=NA)
  
  
  
  
}



