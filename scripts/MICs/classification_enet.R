### Rscript pop_struct_enet.R pheno_file matrix_file st_file
rm(list = ls())
library(OneR)
library(reshape2)
library(scales)
library(ranger)
library(glmnet)
library(caret)
library(doParallel)
#library(doMC)
#registerDoMC(cores = 40)
registerDoParallel(cores = 40)
args=commandArgs(trailingOnly = TRUE)
#setwd("/home/ghepard/progetto_x_londra/data/machine_learning_models/data/real_data/scripts/")
phen_file <- args[1]
#phen_file<-"amikacin.txt"
laststr<-gsub("^.*/", "", phen_file)
antib<-unlist(strsplit(laststr,split = "[.]"))[1]

matrix_file<-args[2]
simulation_type<-args[4]
dd<-read.delim(matrix_file,sep='\t')
#dd<-read.delim("../matrix_lineages/pres_abs_matrix.txt",sep='\t')
colnames(dd)<-gsub("X","",colnames(dd))
colnames(dd)<-gsub(".fna","",colnames(dd))
rownames(dd)<-dd$ID
dd$ID<-NULL
###
tdd<-as.data.frame(t(dd))
###
phen<-read.delim(phen_file,sep='')

tdd$ID<-rownames(tdd)
colnames(phen)[1]<-"ID"
colnames(phen)[2]<-"MIC"
###

mm<-merge(tdd,phen,by="ID")
rownames(mm)<-mm$ID
###



sorted_bins<-sort(unique(as.numeric(log2(phen$MIC))))



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



#mm$MIC<-as.numeric(as.character(mm$MIC))


pred<-dummyVars(formula = "~.",data = mm)

df<-data.frame(predict(pred,newdata = mm),row.names=ids)
colnames(df)<-gsub("PPPP","PP",colnames(df))




###
#set.seed(12345)
###

### SELECT ONLY MIC WITH OCCURRENCE > 10


#mic_count<-as.data.frame(melt(table(phen$MIC)))
#mic_count<-mic_count[which(mic_count$value>10),]
#mic_count_selected<-mic_count$Var1
#df$MIC<-ifelse(as.character(df$MIC) %in% as.character(mic_count_selected),df$MIC,NA)





df<-df[!is.na(df$MIC),]

#print (head(df$MIC))



df$MIC<-log2(df$MIC)


trainIndex <- createDataPartition(df$MIC, p = .70, list = FALSE, times=1) ###partizione 70%train-30%test


train_y<-df[trainIndex,]$MIC
test_y<-df[-trainIndex,]$MIC


train_df<-df[trainIndex,]
train_y<-train_df$MIC

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



#### RANDOM FOREST 



start_time_runtime <- Sys.time()



cv_elnet <- cv.glmnet(x=train_df, y=as.factor(train_y), alpha=0.01,
                      type.measure="class", family="multinomial", keep = TRUE,parallel = TRUE,weights = reweight)




start_time_prediction <- Sys.time()

enet_pred_class<-predict(cv_elnet,s="lambda.1se",newx = test_df,type="class")

end_time_prediction <-Sys.time()

prediction_time <- end_time_prediction - start_time_prediction
prediction_time <- as.character(prediction_time)

end_time_runtime <- Sys.time()

runtime_analysis<- end_time_runtime - start_time_runtime
runtime_analysis<-as.character(runtime_analysis)

tab_err0<-as.data.frame(table(mapply(loosely_correct, correct = test_y, predicted = enet_pred_class, MoreArgs = list(bins = sorted_bins, error = 0))))
acc_err0<-tab_err0[which(tab_err0$Var1=="TRUE"),"Freq"]/length(test_y)
tab_err1<-as.data.frame(table(mapply(loosely_correct, correct = test_y, predicted = enet_pred_class, MoreArgs = list(bins = sorted_bins, error = 1))))
acc_err1<-tab_err1[which(tab_err1$Var1=="TRUE"),"Freq"]/length(test_y)


cm<-confusionMatrix(as.factor(enet_pred_class),as.factor(test_y))


print (tab_err1)
print (cm)


cmtab<-as.data.frame(cm$byClass)
cmtab[] <- lapply(cmtab, function(x) format(round(x, 2), nsmall = 2))
colnames(cmtab)<-gsub(" ","_",colnames(cmtab))
bacc<-mean(as.numeric(cmtab$Balanced_Accuracy))

df_stat_enet<-data.frame()

df_stat_enet[paste(n,"bins",sep="_"),1]<-bacc
df_stat_enet[paste(n,"bins",sep="_"),2]<-acc_err0
df_stat_enet[paste(n,"bins",sep="_"),3]<-acc_err1
df_stat_enet[paste(n,"bins",sep="_"),4]<-antib
df_stat_enet[paste(n,"bins",sep="_"),5]<-"enet"
df_stat_enet[paste(n,"bins",sep="_"),6]<-prediction_time
df_stat_enet[paste(n,"bins",sep="_"),7]<-runtime_analysis
df_stat_enet[paste(n,"bins",sep="_"),8]<-n

colnames(df_stat_enet)<-c("bACC","ACC_err0","ACC_err1","antib","model","prediction_time","train_test_time","classes")



write.table(df_stat_enet,sprintf("%s_classification_summary_log2.txt",antib),quote = F,sep = "\t",row.names = F)


