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
registerDoParallel(cores = 20)
#setwd("/home/ghepard/progetto_x_londra/data/machine_learning_models/data/real_data/")
#
args=commandArgs(trailingOnly = TRUE)
# setwd("/home/ghepard/progetto_x_londra/data/machine_learning_models/data/real_data/scripts/")
phen_file <- args[1]
#phen_file<-"amikacin.txt"
laststr<-gsub("^.*/", "", phen_file)
antib<-unlist(strsplit(laststr,split = "[.]"))[1]



matrix_file<-args[2]
# simulation_type<-args[4]
dd<-read.delim(matrix_file,sep='\t')
#dd<-read.delim("../matrix_lineages/pres_abs_matrix.txt",sep='\t')
colnames(dd)<-gsub("X","",colnames(dd))
colnames(dd)<-gsub(".fna","",colnames(dd))
rownames(dd)<-dd$ID
dd$ID<-NULL
# ###
tdd<-as.data.frame(t(dd))
#tdd[tdd>1]<-1


#tdd<-tdd[,c(1:50)]


phen<-read.delim(phen_file,sep='',header=F)
#phen<-read.delim("../../../../MIC/amikacin.txt",sep='')

tdd$ID<-rownames(tdd)
colnames(phen)[1]<-"ID"
colnames(phen)[2]<-"MIC"

###
mm<-merge(tdd,phen,by="ID")
rownames(mm)<-mm$ID
###
st_file<-args[3]
#st<-read.delim("../matrix_lineages/renamed_lineages.txt",header = F)
st<-read.delim(st_file,header = F)
#st$V2<-gsub("[.]","_",st$V2)
st$V1<-as.character(st$V1)
st$V2<-as.character(st$V2)
colnames(st)<-c("ID","PP")

###
#mm$ID<-gsub("[.]","_",mm$ID)
mm<-merge(mm,st,by="ID")
ids<-mm$ID
mm$ID<-NULL
####### ENCODING

mm$MIC<-as.numeric(as.character(mm$MIC))

pred<-dummyVars(formula = "~.",data = mm)



df<-data.frame(predict(pred,newdata = mm),row.names=ids)
colnames(df)<-gsub("PPPP","PP",colnames(df))



###
set.seed(1234)
###


df<-df[!is.na(df$MIC),]


df$MIC<-log2(df$MIC)


trainIndex <- createDataPartition(df$MIC, p = .70, list = FALSE, times=1) ###partizione 70%train-30%test

#table(bin(df$MIC,method = "length",5))

train_y<-df[trainIndex,]$MIC
test_y<-df[-trainIndex,]$MIC
#train_yy<-bin(train_y,nbins=4,method = "content")  


train_df<-df[trainIndex,]
test_df<-df[-trainIndex,]
test_df$MIC<-NULL




#1/(size of cluster) for every sample in the cluster
st_frst_pos<-min(grep("PP",colnames(train_df)))
st_last_pos<-max(grep("PP",colnames(train_df)))
list_st<-grep("PP",colnames(train_df),value = T)
###

freqs<-data.frame()
for (i in list_st) {
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



################################ RANDOM FOREST



r2_df_train<-data.frame()
r2_df_test<-data.frame()

df_stat_rf<-data.frame()
###


start_time_runtime <- Sys.time()

rf_perm_regr <- ranger(MIC~., data = train_df, num.trees = 500, importance="impurity", write.forest = TRUE,case.weights=reweight,num.threads = 20)

r2_train<-rf_perm_regr$r.squared

start_time_prediction<-Sys.time()

pred_rf<-predict(rf_perm_regr,test_df,type = "response")
end_time_prediction<-Sys.time()

prediction_time <- end_time_prediction - start_time_prediction
prediction_time<-as.character(prediction_time)

end_time_runtime <- Sys.time()

runtime_analysis<- end_time_runtime - start_time_runtime
runtime_analysis<-as.character(runtime_analysis)


sst <- sum((test_y - mean(test_y))^2)
sse <- sum((pred_rf$predictions - test_y)^2)
r2_test <- 1 - sse / sst

df_stat_rf[1,c("R2_train","R2_test","antibiotic","model","type","prediction_time","train_test_time")]<-
  c(r2_train,r2_test,antib,"rf","regression",prediction_time,runtime_analysis)



import_nopval<-as.data.frame(importance(rf_perm_regr))
colnames(import_nopval)<-"coef_importance"

write.table(df_stat_rf,sprintf("%s_summary_log2.txt",antib),col.names = NA,quote = F,sep = "\t")

write.table(import_nopval,sprintf("%s_rf_regr_importance.txt",antib),col.names = NA,quote = F,sep = "\t")



################################ ELASTIC NET 



r2_df_train<-data.frame()
r2_df_test<-data.frame()

df_stat_enet<-data.frame()

df$MIC<-NULL
train_df<-data.matrix(df[trainIndex,])
test_df<-data.matrix(df[-trainIndex,])


cv_elnet<- cv.glmnet(x=train_df[,colnames(train_df)!="MIC"], y=train_y, alpha=0.01,
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


sst <- sum((test_y - mean(test_y))^2)
sse <- sum((pred_enet - test_y)^2)
r2_test <- 1 - sse / sst
r2_df_test[paste("alpha","0.01",sep="_"),1]<-r2_test
colnames(r2_df_test)<-"R2_test" 


df_stat_enet[1,c("R2_train","R2_test","antibiotic","model","type","prediction_time","train_test_time")]<-
  c(r2_train,r2_test,antib,"enet","regression",prediction_time,runtime_analysis)


coef_elnet<-coef(cv_elnet,s="lambda.1se")

coef_elnet<-as.data.frame(data.matrix(coef_elnet))
#coef_df_summary<-as.data.frame(as.matrix(coef_elnet))
colnames(coef_elnet)<-"coef_importance"


#write.table(df_stat_enet,sprintf("%s_summary_log10.txt",antib),col.names = NA,quote = F,sep = "\t")

write.table(import_nopval,sprintf("%s_enet_regr_importance.txt",antib),col.names = NA,quote = F,sep = "\t")

######

final_summary<-rbind(df_stat_enet,df_stat_rf)

write.table(final_summary,sprintf("%s_regression_summary_log2.txt",antib),quote = F,sep = "\t",row.names = F)

