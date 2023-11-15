### Rscript pop_struct_enet.R pheno_file matrix_file st_file
rm(list = ls())

library(OneR)
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
#dd<-read.delim("../../data/machine_learning_models/data/3run/input/pres_abs_matrix.txt",sep='\t')

dd<-read.delim(matrix_file,sep='\t')
colnames(dd)<-gsub("X","",colnames(dd))
colnames(dd)<-gsub(".fna","",colnames(dd))

dd$ID<-gsub("_out","",dd$ID)
dd$ID<-gsub(".fna.out","",dd$ID)
dd$ID<-gsub(".out","",dd$ID)

rownames(dd)<-dd$ID
dd$ID<-NULL
###
tdd<-as.data.frame(t(dd))


set.seed(1234)

phen<-read.delim(phen_file,sep='',header=F)
#phen<-read.delim("../../data/machine_learning_models/data/3run/1000snps/simul_1000SNPs_EF_10.0_herit_0.9_traitprev_0.1.phen",sep='',header=F)
phen$V1<-NULL
phen$V2<-gsub("X","",phen$V2)
phen$V2<-gsub(".fna","",phen$V2)
tdd$ID<-rownames(tdd)
colnames(phen)[1]<-"ID"
colnames(phen)[2]<-"MIC"
###
mm<-merge(tdd,phen,by="ID")
rownames(mm)<-mm$ID
###

st_file<-args[3]
st<-read.delim(st_file,header = F)
#st<-read.delim("../../data/machine_learning_models/data/3run/1000snps/renamed_lineages.txt",header = F)

#st$V2<-gsub("[.]","_",st$V2)
st$V1<-as.character(st$V1)
st$V2<-as.character(st$V2)
colnames(st)<-c("ID","PP")
###
#mm$ID<-gsub("[.]","_",mm$ID)



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


df$sim_cont_mic<-rescale(df$MIC,to=c(0.25,16))
#df$sim_cont_mic<-bin(df$sim_cont_mic, nbins = 4,method = "content")
df$MIC<-NULL

###

###
trainIndex <- createDataPartition(df$sim_cont_mic,  p = .70, list = FALSE, times=1) ###partizione 70%train-30%test


train_y<-df[trainIndex,]$MIC
test_y<-df[-trainIndex,]$MIC
df$sim_cont_mic<-NULL
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
r2_train<-cv_elnet$glmnet.fit$dev.ratio[which(cv_elnet$glmnet.fit$lambda == cv_elnet$lambda.min)] ### R2 train set
r2_df_train[paste("alpha","0.01",sep="_"),1]<-r2_train ###DF CREATION WITH ALPHA AND R2
colnames(r2_df_train)<-"R2_train" 

start_time_prediction<-Sys.time()

pred_enet<-predict(cv_elnet,s="lambda.1se",newx = test_df) ### PREDICTION

end_time_prediction<-Sys.time()

prediction_time <- end_time_prediction - start_time_prediction
prediction_time<-as.character(prediction_time)

end_time_runtime <- Sys.time()

runtime_analysis<- end_time_runtime - start_time_runtime
runtime_analysis<-as.character(runtime_analysis)

write.table(df_pred,sprintf("%s_prediction_performance_R2_enet_normal_regression.txt",phen_file),
            sep = '\t',quote = F,col.names=NA)


sst <- sum((test_y - mean(test_y))^2)
sse <- sum((pred_enet - test_y)^2)
r2_test <- 1 - sse / sst
r2_df_test[paste("alpha","0.01",sep="_"),1]<-r2_test
colnames(r2_df_test)<-"R2_test" 


laststr<-gsub("^.*/", "", phen_file)
h2<-unlist(strsplit(laststr,split = "_"))[6]



df_stat[1,c("R2_train","R2_test","simulation","model","type","h2","prediction_time","train_test_time")]<-c(r2_train,r2_test,simulation_type,"enet",                                                "regression",h2,prediction_time,runtime_analysis)


write.table(df_stat,sprintf("%s_R2_enet_normal_regression.txt",phen_file),
            sep = '\t',quote = F,row.names=F)


coef_elnet<-coef(cv_elnet,s="lambda.1se")


coef_elnet<-as.data.frame(data.matrix(coef_elnet))

 #coef_df_summary<-as.data.frame(as.matrix(coef_elnet))
colnames(coef_elnet)<-"coef_importance"
    


write.table(coef_elnet,sprintf("%s_summary_coef_R2_enet_normal_regression.txt",phen_file),
            sep = '\t',quote = F,col.names=NA)










