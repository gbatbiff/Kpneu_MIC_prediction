### Rscript pop_struct_enet.R pheno_file matrix_file st_file
rm(list = ls())
library(OneR)
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
#dd<-read.delim("pres_abs_matrix.txt",sep='\t')

colnames(dd)<-gsub("X","",colnames(dd))
colnames(dd)<-gsub(".fna","",colnames(dd))
rownames(dd)<-dd$ID
dd$ID<-NULL
###
tdd<-as.data.frame(t(dd))
tdd[tdd>1]<-1
#tdd<-tdd[,c(1:50)]
###

phen<-read.delim(phen_file,sep='',header=F)
#phen<-read.delim("simul_3betalact_EF_2.5_herit_0.9_traitprev_0.2.phen",sep='',header=F)
phen$V1<-NULL
phen$V2<-gsub("X","",phen$V2)
phen$V2<-gsub(".fna","",phen$V2)
phen$V2<-gsub(".out","",phen$V2)
phen$V2<-gsub("_out","",phen$V2)


tdd$ID<-rownames(tdd)
colnames(phen)[1]<-"ID"
colnames(phen)[2]<-"MIC"
###
mm<-merge(tdd,phen,by="ID")
rownames(mm)<-mm$ID
###
st_file<-args[3]
st<-read.delim(st_file,header = F)
#st<-read.delim("renamed_lineages.txt",header=F)

st$V2<-gsub("-","_",st$V2)
st$V1<-as.character(st$V1)
st$V2<-as.character(st$V2)
colnames(st)<-c("ID","PP")

tdd$ID<-rownames(tdd)


###
mm<-merge(mm,st,by="ID")
ids<-mm$ID
mm1<-mm[,c("group_20046","ubiE","rpmJ","PP","MIC")] ###dummy df to avoid mem limit

pred<-dummyVars(formula = "~.",data = mm1)
df<-data.frame(predict(pred,newdata = mm1),row.names=ids)
colnames(df)<-gsub("PPPP","PP",colnames(df))

pp_lineages_df<-df[,4:ncol(df)]

df<-cbind(mm,pp_lineages_df)


df$sim_cont_mic<-rescale(df$MIC,to=c(0.25,16))
#df$sim_cont_mic<-bin(df$sim_cont_mic, nbins = 4,method = "content")
df$MIC<-NULL

###
set.seed(1234)
###
trainIndex <- createDataPartition(df$sim_cont_mic,  p = .70, list = FALSE, times=1) ###partizione 70%train-30%test

train_y<-df[trainIndex,]$sim_cont_mic
test_y<-df[-trainIndex,]$sim_cont_mic
  
  ###
  
train_df<-df[trainIndex,]
test_df<-df[-trainIndex,]
test_df$sim_cont_mic<-NULL
train_df$MIC<-NULL  
 


#1/(size of cluster) for every sample in the cluster
st_frst_pos<-min(grep("PP",colnames(train_df)))
st_last_pos<-max(grep("PP",colnames(train_df)))
list_st<-grep("PP",colnames(train_df),value = T)
###

freqs<-data.frame()
for (i in list_st[-1]) {
  stdf<-table(train_df[,i])	
  
    if (dim(stdf)>1){
  stdf<-as.data.frame(stdf)
  stdf$type<-i
  stdf<-stdf[which(stdf$Var1==1),]
  stdf$Var1<-NULL
  freqs[i,1]<-stdf$Freq
  }
}

freqs$reweight<-1/freqs$V1
### mettere nella variabile weight 1 per tutti gli snp e per i st i valori del vettore reweight
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
	
start_time_runtime <- Sys.time()

rf_perm_regr <- ranger(y=train_df$sim_cont_mic, x=train_df[colnames(train_df)!="sim_cont_mic"], num.trees = 500, importance="impurity", write.forest = TRUE,case.weights=reweight,num.threads=50)



r2_train<-rf_perm_regr$r.squared

start_time_prediction<-Sys.time()

pred_rf<-predict(rf_perm_regr,test_df,type = "response")
sst <- sum((test_y - mean(test_y))^2)
sse <- sum((pred_rf$predictions - test_y)^2)
r2_test <- 1 - sse / sst

end_time_prediction<-Sys.time()

prediction_time <- end_time_prediction - start_time_prediction
prediction_time<-as.character(prediction_time)

end_time_runtime <- Sys.time()

runtime_analysis<- end_time_runtime - start_time_runtime
runtime_analysis<-as.character(runtime_analysis)



laststr<-gsub("^.*/", "", phen_file)
h2<-unlist(strsplit(laststr,split = "_"))[6]
#h2<-unlist(strsplit(phen_file,split = "_"))[7]
  

df_stat[1,c("R2_train","R2_test","simulation","model","type","h2","prediction_time","train_test_time")]<-c(r2_train,r2_test,simulation_type,"RF",                                                "regression",h2,prediction_time,runtime_analysis)


  
  write.table(df_stat,sprintf("%s_RF_R2_continuous.txt",phen_file),
            sep = '\t',quote = F,row.names=F)
  

check_imp<-as.data.frame(importance(rf_perm_regr))




write.table(check_imp,sprintf("%s_RF_imp_pred_R2_continuous.txt",phen_file),
            sep = '\t',quote = F,col.names=NA)




