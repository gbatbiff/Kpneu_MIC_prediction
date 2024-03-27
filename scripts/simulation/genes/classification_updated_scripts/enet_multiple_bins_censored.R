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
colnames(dd)<-gsub("X","",colnames(dd))
colnames(dd)<-gsub(".fna","",colnames(dd))
rownames(dd)<-dd$ID
dd$ID<-NULL
###
tdd<-as.data.frame(t(dd))


n_bins<-c(2,4,6)

set.seed(1234)

for (n in n_bins){



  phen<-read.delim(phen_file,sep='',header=F)
  #phen<-read.delim("../simulated_phenotypes/simul_1000SNPs_EF_2.5_herit_0.9_traitprev_0.1.phen",sep='',header=F)
  phen$V1<-NULL
  phen$V2<-gsub("X","",phen$V2)
  phen$V2<-gsub(".fna","",phen$V2)
  colnames(phen)<-c("sample","MIC")
  phen$sim_ord_mic<-rescale(phen$MIC,to=c(0.25,16))
  qt<-quantile(phen$sim_ord_mic, probs = c(0.025,0.975))
  phen$sim_ord_mic_adjusted<-ifelse(phen$sim_ord_mic<qt[[1]],qt[[1]],phen$sim_ord_mic)
  phen$sim_ord_mic_adjusted<-ifelse(phen$sim_ord_mic_adjusted>qt[[2]],qt[[2]],phen$sim_ord_mic_adjusted)
  
  
  phen$sim_ord_mic<-bin(phen$sim_ord_mic_adjusted, nbins = n,method = "length") 
  phen$sim_ord_mic_adjusted<-NULL
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
  
  sorted_bins<-sort(unique(as.numeric(phen_mid$midpoint_mic)))
  
  phen_mid$midpoint_mic<-as.numeric(as.character(phen_mid$midpoint_mic))
  
  #last<-tail(sort(table(phen_mid$midpoint_mic),decreasing=T),3)
  #last<-as.numeric(rownames(last))
  #censored_class<-sorted_bins[4] ### ULTIMA CLASSE BILANCIATA
  #phen_mid$midpoint_mic<-ifelse(phen_mid$midpoint_mic %in% last,censored_class,phen_mid$midpoint_mic) ### SE CLASSI APPARTENGONO A SBILANCIATE INCORPORA IN ULTIMA (censor)


  
  colnames(phen_mid)<-c("ID","MIC")
  colnames(tdd)<-gsub("~~~","_",colnames(tdd))
  
  tdd$ID<-rownames(tdd)
  ###

  mm<-merge(tdd,phen_mid,by="ID")
  rownames(mm)<-mm$ID

    st_file<-args[3]
  st<-read.delim(st_file,header = F)
  #st<-read.delim("renamed_lineages.txt",header = F)
  st$V1<-as.character(st$V1)
  st$V2<-as.character(st$V2)
  colnames(st)<-c("ID","PP")
  ###
  
  mm<-merge(mm,st,by="ID")
  ids<-mm$ID
  
  mm1<-mm[,c(11960:11964)]
  
  
  
  
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

  
  
  
  if (n==2){

    df_stat<-data.frame()

    start_time_runtime <- Sys.time()

    cv_elnet <- cv.glmnet(x=train_df, y=as.factor(train_y), alpha=0.01,
                          type.measure="class", family="binomial", keep = TRUE,parallel = TRUE,weights = reweight)
    
    start_time_prediction <- Sys.time()
    
    pred_class<-predict(cv_elnet,s="lambda.1se",newx = test_df,type="class")
    
    end_time_prediction <-Sys.time()
    
    prediction_time <- end_time_prediction - start_time_prediction
    prediction_time <- as.character(prediction_time)
    
    end_time_runtime <- Sys.time()

    runtime_analysis<- end_time_runtime - start_time_runtime
    runtime_analysis<-as.character(runtime_analysis)
    
    tab_err0<-as.data.frame(table(mapply(loosely_correct, correct = test_y, predicted = pred_class, MoreArgs = list(bins = sorted_bins, error = 0))))
    acc_err0<-tab_err0[which(tab_err0$Var1=="TRUE"),"Freq"]/length(test_y)
    tab_err1<-as.data.frame(table(mapply(loosely_correct, correct = test_y, predicted = pred_class, MoreArgs = list(bins = sorted_bins, error = 1))))
    acc_err1<-tab_err1[which(tab_err1$Var1=="TRUE"),"Freq"]/length(test_y)
    
    
 cm<-confusionMatrix(as.factor(pred_class),as.factor(test_y))
    
### COMPUTE STATISTICS
    cmtab<-as.data.frame(cm$byClass)
    cmtab[] <- lapply(cmtab, function(x) format(round(x, 2), nsmall = 2))
    
    cmtab<-as.data.frame(t(cmtab))
    colnames(cmtab)<-gsub(" ","_",colnames(cmtab))
    bacc<-as.character(cmtab$Balanced_Accuracy)
    Precision<-cmtab$Precision
    Recall<-cmtab$Recall
    F1<-cmtab$F1
    spec<-cmtab$Specificity
    sens<-cmtab$Sensitivity 
    laststr<-gsub("^.*/", "", phen_file)
    h2<-unlist(strsplit(laststr,split = "_"))[6]
   


    df_stat[paste(n,"bins",sep="_"),1]<-bacc
    df_stat[paste(n,"bins",sep="_"),2]<-acc_err0
    df_stat[paste(n,"bins",sep="_"),3]<-acc_err1
    df_stat[paste(n,"bins",sep="_"),4]<-simulation_type
    df_stat[paste(n,"bins",sep="_"),5]<-"enet"
    df_stat[paste(n,"bins",sep="_"),6]<-"censored"
    df_stat[paste(n,"bins",sep="_"),7]<-h2
    df_stat[paste(n,"bins",sep="_"),8]<-prediction_time
    df_stat[paste(n,"bins",sep="_"),9]<-runtime_analysis
    df_stat[paste(n,"bins",sep="_"),10]<-Precision
    df_stat[paste(n,"bins",sep="_"),11]<-Recall
    df_stat[paste(n,"bins",sep="_"),12]<-F1
    df_stat[paste(n,"bins",sep="_"),13]<-spec
    df_stat[paste(n,"bins",sep="_"),14]<-sens


    colnames(df_stat)<-c("bACC","ACC_err0","ACC_err1","simulation","model","type","h2","prediction_time","train_test_time", "precision", "recall", "F1","Specificity", "Sensitivity")

    write.table(df_stat,sprintf("%s_ACC_enet_bins_%s_censored.txt",phen_file,n),
                sep = '\t',quote = F,col.names=NA)
    

       
  }
  
  
  
  
  if (n==4){
 

    df_stat<-data.frame()
	

    start_time_runtime <- Sys.time()


    cv_elnet <- cv.glmnet(x=train_df, y=as.factor(train_y), alpha=0.01,
                          type.measure="class", family="multinomial", keep = TRUE,parallel = TRUE,weights = reweight)
    
    start_time_prediction <- Sys.time()
    
    pred_class<-predict(cv_elnet,s="lambda.1se",newx = test_df,type="class")
    
    end_time_prediction <-Sys.time()
    
    prediction_time <- end_time_prediction - start_time_prediction
    prediction_time <- as.character(prediction_time)
    
    end_time_runtime <- Sys.time()

    runtime_analysis<- end_time_runtime - start_time_runtime
    runtime_analysis<-as.character(runtime_analysis)
    
    tab_err0<-as.data.frame(table(mapply(loosely_correct, correct = test_y, predicted = pred_class, MoreArgs = list(bins = sorted_bins, error = 0))))
    acc_err0<-tab_err0[which(tab_err0$Var1=="TRUE"),"Freq"]/length(test_y)
    tab_err1<-as.data.frame(table(mapply(loosely_correct, correct = test_y, predicted = pred_class, MoreArgs = list(bins = sorted_bins, error = 1))))
    acc_err1<-tab_err1[which(tab_err1$Var1=="TRUE"),"Freq"]/length(test_y)
    
    cm<-confusionMatrix(as.factor(pred_class),as.factor(test_y))
 
     ### COMPUTE STATISTICS
    cmtab<-as.data.frame(cm$byClass)
    cmtab[] <- lapply(cmtab, function(x) format(round(x, 2), nsmall = 2))
    colnames(cmtab)<-gsub(" ","_",colnames(cmtab))


    bacc<-mean(as.numeric(cmtab$Balanced_Accuracy))
    spec<-mean(as.numeric(cmtab$Specificity))
    sens<-mean(as.numeric(cmtab$Sensitivity))


    f1<-gsub("NA", "0", cmtab$F1)
    f1<-mean(as.numeric(f1))
    recall<-gsub("NA", "0", cmtab$Recall)
    recall<-mean(as.numeric(recall))
    precision<-gsub("NA", "0", cmtab$Precision)
    precision<-mean(as.numeric(precision))
    laststr<-gsub("^.*/", "", phen_file)
    h2<-unlist(strsplit(laststr,split = "_"))[6]

    df_stat[paste(n,"bins",sep="_"),1]<-bacc
    df_stat[paste(n,"bins",sep="_"),2]<-acc_err0
    df_stat[paste(n,"bins",sep="_"),3]<-acc_err1
    df_stat[paste(n,"bins",sep="_"),4]<-simulation_type
    df_stat[paste(n,"bins",sep="_"),5]<-"enet"
    df_stat[paste(n,"bins",sep="_"),6]<-"censored"
    df_stat[paste(n,"bins",sep="_"),7]<-h2
    df_stat[paste(n,"bins",sep="_"),8]<-prediction_time
    df_stat[paste(n,"bins",sep="_"),9]<-runtime_analysis
    df_stat[paste(n,"bins",sep="_"),10]<-precision
    df_stat[paste(n,"bins",sep="_"),11]<-recall
    df_stat[paste(n,"bins",sep="_"),12]<-f1
    df_stat[paste(n,"bins",sep="_"),13]<-spec
    df_stat[paste(n,"bins",sep="_"),14]<-sens
   
    colnames(df_stat)<-c("bACC","ACC_err0","ACC_err1","simulation","model","type","h2","prediction_time","train_test_time", "precision", "recall", "F1","Specificity", "Sensitivity")
    

    
    write.table(df_stat,sprintf("%s_ACC_enet_bins_%s_censored.txt",phen_file,n),
                sep = '\t',quote = F,col.names = NA)
    
    

     
  
    
  }
  
  
  if (n==6){
    
    df_stat<-data.frame()

    
    
    cv_elnet <- cv.glmnet(x=train_df, y=as.factor(train_y), alpha=0.01,
                          type.measure="class", family="multinomial", keep = TRUE,parallel = TRUE,weights = reweight)
    
    
    start_time_prediction <- Sys.time()
    
    pred_class<-predict(cv_elnet,s="lambda.1se",newx = test_df,type="class")
    
    end_time_prediction <- Sys.time()
    
    prediction_time <- end_time_prediction - start_time_prediction
     
    
    tab_err0<-as.data.frame(table(mapply(loosely_correct, correct = test_y, predicted = pred_class, MoreArgs = list(bins = sorted_bins, error = 0))))
    acc_err0<-tab_err0[which(tab_err0$Var1=="TRUE"),"Freq"]/length(test_y)
    tab_err1<-as.data.frame(table(mapply(loosely_correct, correct = test_y, predicted = pred_class, MoreArgs = list(bins = sorted_bins, error = 1))))
    acc_err1<-tab_err1[which(tab_err1$Var1=="TRUE"),"Freq"]/length(test_y)
    
    
    cm<-confusionMatrix(as.factor(pred_class),as.factor(test_y))
 
     ### COMPUTE STATISTICS
    cmtab<-as.data.frame(cm$byClass)
    cmtab[] <- lapply(cmtab, function(x) format(round(x, 2), nsmall = 2))
    colnames(cmtab)<-gsub(" ","_",colnames(cmtab))


    bacc<-mean(as.numeric(cmtab$Balanced_Accuracy))
    spec<-mean(as.numeric(cmtab$Specificity))
    sens<-mean(as.numeric(cmtab$Sensitivity))


    f1<-gsub("NA", "0", cmtab$F1)

    f1<-mean(as.numeric(f1))
    recall<-gsub("NA", "0", cmtab$Recall)
    recall<-mean(as.numeric(recall))
    precision<-gsub("NA", "0", cmtab$Precision)
    precision<-mean(as.numeric(precision))
    laststr<-gsub("^.*/", "", phen_file)
    h2<-unlist(strsplit(laststr,split = "_"))[6]

    df_stat[paste(n,"bins",sep="_"),1]<-bacc
    df_stat[paste(n,"bins",sep="_"),2]<-acc_err0
    df_stat[paste(n,"bins",sep="_"),3]<-acc_err1
    df_stat[paste(n,"bins",sep="_"),4]<-simulation_type
    df_stat[paste(n,"bins",sep="_"),5]<-"enet"
    df_stat[paste(n,"bins",sep="_"),6]<-"censored"
    df_stat[paste(n,"bins",sep="_"),7]<-h2
    df_stat[paste(n,"bins",sep="_"),8]<-prediction_time
    df_stat[paste(n,"bins",sep="_"),9]<-runtime_analysis
    df_stat[paste(n,"bins",sep="_"),10]<-precision
    df_stat[paste(n,"bins",sep="_"),11]<-recall
    df_stat[paste(n,"bins",sep="_"),12]<-f1
    df_stat[paste(n,"bins",sep="_"),13]<-spec
    df_stat[paste(n,"bins",sep="_"),14]<-sens
   
    colnames(df_stat)<-c("bACC","ACC_err0","ACC_err1","simulation","model","type","h2","prediction_time","train_test_time", "precision", "recall", "F1","Specificity", "Sensitivity")
    

    
    write.table(df_stat,sprintf("%s_ACC_enet_bins_%s_censored.txt",phen_file,n),
                sep = '\t',quote = F,col.names = NA)
    
      
     
  }
  
  
}
