#Getting r-sq of insample and out of sample prediction
rm(list=ls())

#-----------------------------------For pooled data------------------------------------------------------------
#Load data file's preds insamp and outsamp
filename = paste('AllBTData_Post_IF', 'mUGM','3','zaUTerusign_var')

load(paste('SCMuscimol_Study_datafiles/datafiles/prediction stats files/',filename,'_preds',sep=""))
load(paste('SCMuscimol_Study_datafiles/datafiles/insamp stats files/',filename,'_insamp',sep=""))
load(paste('SCMuscimol_Study_datafiles/datafiles/outsamp stats files/',filename,'_outsamp',sep=""))


#6.run this once youve done all the datafiles and the pooled data has all the injections data for the subject
#preds vs insamp
SE_insamppred_rtmean = list()
SE_insamppred_rtmedian = list()
SE_insamppred_rtperc25 = list()
SE_insamppred_rtperc75 = list()
SE_insamppred_acc = list()
SE_insamp_rtmean = list()
SE_insamp_rtmedian = list()
SE_insamp_rtperc25 = list()
SE_insamp_rtperc75=list()
SE_insamp_acc = list()

SE_outsamppred_rtmean = list()
SE_outsamppred_rtmedian = list()
SE_outsamppred_rtperc25 = list()
SE_outsamppred_rtperc75 = list()
SE_outsamppred_acc = list()
SE_outsamp_rtmean = list()
SE_outsamp_rtmedian = list()
SE_outsamp_rtperc25 = list()
SE_outsamp_rtperc75=list()
SE_outsamp_acc = list()

for (c in 1:dim(preds['ncohs'])[1]) { #for each coh
  #cumulative list of (observations - predictions)^2 -- insamp
  SE_insamppred_rtmean = c(SE_insamppred_rtmean,(df_observ_insamp[['rt_mean']][[c]] - preds[['rt_mean']][[c]])^2)
  SE_insamppred_rtmedian = c(SE_insamppred_rtmedian, (df_observ_insamp[['rt_median']][[c]] - preds[['rt_median']][[c]])^2)
  SE_insamppred_rtperc25 = c(SE_insamppred_rtperc25,(df_observ_insamp[['rt_perc_25']][[c]] - preds[['rt_perc_25']][[c]])^2)
  SE_insamppred_rtperc75 = c(SE_insamppred_rtperc75, (df_observ_insamp[['rt_perc_75']][[c]] - preds[['rt_perc_75']][[c]])^2)
  SE_insamppred_acc = c(SE_insamppred_acc, (df_observ_insamp[['acc']][[c]] - preds[['acc']][[c]])^2)
  #cumulative list of (observations - observation_mean)^2 -- insamp
  SE_insamp_rtmean = c(SE_insamp_rtmean,(df_observ_insamp[['rt_mean']][[c]]  - mean(df_observ_insamp[['rt_mean']]))^2)
  SE_insamp_rtmedian = c(SE_insamp_rtmedian,(df_observ_insamp[['rt_median']][[c]] - mean(df_observ_insamp[['rt_median']]))^2)
  SE_insamp_rtperc25 = c(SE_insamp_rtperc25,(df_observ_insamp[['rt_perc_25']][[c]]  - mean(df_observ_insamp[['rt_perc_25']]))^2)
  SE_insamp_rtperc75 = c(SE_insamp_rtperc75,(df_observ_insamp[['rt_perc_75']][[c]] - mean(df_observ_insamp[['rt_perc_75']]))^2)
  SE_insamp_acc = c(SE_insamp_acc,(df_observ_insamp[['acc']][[c]] - mean(df_observ_insamp[['acc']]))^2)
  
  #cumulative list of (observations - predictions)^2 -- outsamp
  SE_outsamppred_rtmean = c(SE_outsamppred_rtmean,(df_observ_outsamp[['rt_mean']][[c]] - preds[['rt_mean']][[c]])^2)
  SE_outsamppred_rtmedian = c(SE_outsamppred_rtmedian, (df_observ_outsamp[['rt_median']][[c]] - preds[['rt_median']][[c]])^2)
  SE_outsamppred_rtperc25 = c(SE_outsamppred_rtperc25,(df_observ_outsamp[['rt_perc_25']][[c]] - preds[['rt_perc_25']][[c]])^2)
  SE_outsamppred_rtperc75 = c(SE_outsamppred_rtperc75, (df_observ_outsamp[['rt_perc_75']][[c]] - preds[['rt_perc_75']][[c]])^2)
  SE_outsamppred_acc = c(SE_outsamppred_acc, (df_observ_outsamp[['acc']][[c]] - preds[['acc']][[c]])^2)
  #cumulative list of (observations - observation_mean)^2 -- outsamp
  SE_outsamp_rtmean = c(SE_outsamp_rtmean,(df_observ_outsamp[['rt_mean']][[c]]  - mean(df_observ_outsamp[['rt_mean']]))^2)
  SE_outsamp_rtmedian = c(SE_outsamp_rtmedian,(df_observ_outsamp[['rt_median']][[c]] - mean(df_observ_outsamp[['rt_median']]))^2)
  SE_outsamp_rtperc25 = c(SE_outsamp_rtperc25,(df_observ_outsamp[['rt_perc_25']][[c]]  - mean(df_observ_outsamp[['rt_perc_25']]))^2)
  SE_outsamp_rtperc75 = c(SE_outsamp_rtperc75,(df_observ_outsamp[['rt_perc_75']][[c]] - mean(df_observ_outsamp[['rt_perc_75']]))^2)
  SE_outsamp_acc = c(SE_outsamp_acc,(df_observ_outsamp[['acc']][[c]] - mean(df_observ_outsamp[['acc']]))^2)
}

#get Rsq for insample prediction 
Rsq_insamppred_rtmean = 1-(sum(unlist(SE_insamppred_rtmean))/sum(unlist(SE_insamp_rtmean)))
Rsq_insamppred_rtmedian = 1-(sum(unlist(SE_insamppred_rtmedian))/sum(unlist(SE_insamp_rtmedian)))
Rsq_insamppred_rtperc25 = 1-(sum(unlist(SE_insamppred_rtperc25 ))/sum(unlist(SE_insamp_rtperc25)))
Rsq_insamppred_rtperc75 = 1-(sum(unlist(SE_insamppred_rtperc75 ))/sum(unlist(SE_insamp_rtperc75)))
Rsq_insamppred_acc = 1-(sum(unlist(SE_insamppred_acc ))/sum(unlist(SE_insamp_acc)))

#get Rsq for outsample prediction 
Rsq_outsamppred_rtmean = 1-sum(unlist(SE_outsamppred_rtmean))/sum(unlist(SE_outsamp_rtmean))
Rsq_outsamppred_rtmedian = 1-sum(unlist(SE_outsamppred_rtmedian))/sum(unlist(SE_outsamp_rtmedian))
Rsq_outsamppred_rtperc25 = 1-sum(unlist(SE_outsamppred_rtperc25 ))/sum(unlist(SE_outsamp_rtperc25))
Rsq_outsamppred_rtperc75 = 1-sum(unlist(SE_outsamppred_rtperc75 ))/sum(unlist(SE_outsamp_rtperc75))
Rsq_outsamppred_acc = 1-sum(unlist(SE_outsamppred_acc ))/sum(unlist(SE_outsamp_acc))


# For plotting the preds vs insamp
coh = c(-24, -17, -10, -5, -3, 0, 3, 5, 10, 17, 24)
plot(coh, df_observ_insamp$rt_mean, ylim=c(0.7, 1.2), ylab = 'RT Mean (s)', xlab = 'Coherence (%)') 
points(coh, preds$rt_mean, col='red')


plot(coh, df_observ_insamp$rt_median, ylim=c(0.7, 1.2), ylab = 'RT Median (s)', xlab = 'Coherence (%)') 
points(coh, preds$rt_median, col='red')


#plot(coh, df_observ_insamp$rt_perc_25, ylim=c(0.7, 1.2), ylab = 'RT 25th Percentile (s)', xlab = 'Coherence (%)') 
plot(coh, df_observ_insamp$rt_perc_25, ylim=c(0.7, 1.2), ylab = 'RT 25th Percentile (s)', xlab = 'Coherence (%)') 
points(coh, preds$rt_perc_25, col='red')


#plot(coh, df_observ_insamp$rt_perc_75, ylim=c(0.9, 1.5), ylab = 'RT 75th Percentile (s)', xlab = 'Coherence (%)') 
plot(coh, df_observ_insamp$rt_perc_75, ylim=c(0.9, 1.4), ylab = 'RT 75th Percentile (s)', xlab = 'Coherence (%)') 
points(coh, preds$rt_perc_75, col='red')


plot(coh, df_observ_insamp$acc, ylim=c(0, 1), ylab = 'Proportion to IF Choice', xlab = 'Coherence (%)') 
points(coh, preds$acc, col='red')



##--------------------------------------------------------------------------------------------------
#For individual data that have their own predicitions
rm(list=ls())
individual_dataname = c(paste('SP010720GPPre_IF_1','mUGM','1','aU'),paste('SP011620GPPre_IF_1','mUGM','2', 'aU'), paste('SP012120GPPre_IF_1', 'mUGM','1','aU'),paste('SP013020GPPre_IF_1','mUGM','1','aU'),paste('SP020620GPPre_IF_1','mUGM','2','aU'),paste('SP021320GPPre_IF_1','mUGM','1','aU'), paste('SP121719GPPre_IF_1','mUGM','5','aU'))


rtmeanlist_insamp = list()
rtmedianlist_insamp = list()
rtperc25meanlist_insamp = list()
rtperc75meanlist_insamp = list()
accmeanlist_insamp = list()
rtmeanlist_outsamp = list()
rtmedianlist_outsamp = list()
rtperc25meanlist_outsamp = list()
rtperc75meanlist_outsamp = list()
accmeanlist_outsamp = list()
#get total mean of RT measures and acc across all injections first 
for (i in 1:length(individual_dataname)) {
  load(paste('SCMuscimol_Study_datafiles/datafiles/insamp stats files/', individual_dataname[i],'_insamp',sep=""))
  rtmeanlist_insamp = c(rtmeanlist_insamp,df_observ_insamp$rt_mean)
  rtmedianlist_insamp = c(rtmedianlist_insamp, df_observ_insamp$rt_median)
  rtperc25meanlist_insamp = c(rtperc25meanlist_insamp,df_observ_insamp$rt_perc_25)
  rtperc75meanlist_insamp = c(rtperc75meanlist_insamp,df_observ_insamp$rt_perc_75)
  accmeanlist_insamp = c(accmeanlist_insamp,df_observ_insamp$acc)
  load(paste('SCMuscimol_Study_datafiles/datafiles/outsamp stats files/',individual_dataname[i],'_outsamp',sep=""))
  rtmeanlist_outsamp = c(rtmeanlist_outsamp,df_observ_outsamp$rt_mean)
  rtmedianlist_outsamp = c(rtmedianlist_outsamp, df_observ_outsamp$rt_median)
  rtperc25meanlist_outsamp = c(rtperc25meanlist_outsamp,df_observ_outsamp$rt_perc_25)
  rtperc75meanlist_outsamp = c(rtperc75meanlist_outsamp,df_observ_outsamp$rt_perc_75)
  accmeanlist_outsamp = c(accmeanlist_outsamp,df_observ_outsamp$acc)
}
rtmean_mean_insamp = mean(unlist(rtmeanlist_insamp))
rtmedian_mean_insamp  = mean(unlist(rtmedianlist_insamp))
rtperc25_mean_insamp = mean(unlist(rtperc25meanlist_insamp))
rtperc75_mean_insamp = mean(unlist(rtperc75meanlist_insamp))
acc_mean_insamp  = mean(unlist(accmeanlist_insamp))
rtmean_mean_outsamp = mean(unlist(rtmeanlist_outsamp))
rtmedian_mean_outsamp  = mean(unlist(rtmedianlist_outsamp))
rtperc25_mean_outsamp = mean(unlist(rtperc25meanlist_outsamp))
rtperc75_mean_outsamp = mean(unlist(rtperc75meanlist_outsamp))
acc_mean_outsamp  = mean(unlist(accmeanlist_outsamp))


SE_pred_rtmean_insamp = list()
SE_pred_rtmedian_insamp= list()
SE_pred_rtperc25_insamp = list()
SE_pred_rtperc75_insamp = list()
SE_pred_acc_insamp = list()
SE_pred_rtmean_outsamp = list()
SE_pred_rtmedian_outsamp= list()
SE_pred_rtperc25_outsamp = list()
SE_pred_rtperc75_outsamp = list()
SE_pred_acc_outsamp = list()

SE_mean_rtmean_insamp = list()
SE_mean_rtmedian_insamp= list()
SE_mean_rtperc25_insamp = list()
SE_mean_rtperc75_insamp = list()
SE_mean_acc_insamp = list()
SE_mean_rtmean_outsamp = list()
SE_mean_rtmedian_outsamp= list()
SE_mean_rtperc25_outsamp = list()
SE_mean_rtperc75_outsamp = list()
SE_mean_acc_outsamp = list()
for (i in 1:length(individual_dataname)) {
  load(paste('SCMuscimol_Study_datafiles/datafiles/prediction stats files/',individual_dataname[i],'_preds',sep=""))
  load(paste('SCMuscimol_Study_datafiles/datafiles/insamp stats files/',individual_dataname[i],'_insamp',sep=""))
  load(paste('SCMuscimol_Study_datafiles/datafiles/outsamp stats files/',individual_dataname[i],'_outsamp',sep=""))
  
  for (c in 1:dim(df_observ_insamp)[1]){
    SE_pred_rtmean_insamp = c(SE_pred_rtmean_insamp, (df_observ_insamp$rt_mean[c] - preds$rt_mean[c])^2)
    SE_pred_rtmedian_insamp = c(SE_pred_rtmedian_insamp, (df_observ_insamp$rt_median[c] - preds$rt_median[c])^2)
    SE_pred_rtperc25_insamp = c(SE_pred_rtperc25_insamp, (df_observ_insamp$rt_perc_25[c] - preds$rt_perc_25[c])^2)
    SE_pred_rtperc75_insamp = c(SE_pred_rtperc75_insamp, (df_observ_insamp$rt_perc_75[c] - preds$rt_perc_75[c])^2)
    SE_pred_acc_insamp = c(SE_pred_acc_insamp, (df_observ_insamp$acc[c] - preds$acc[c])^2)
    SE_pred_rtmean_outsamp = c(SE_pred_rtmean_outsamp, (df_observ_outsamp$rt_mean[c] - preds$rt_mean[c])^2)
    SE_pred_rtmedian_outsamp = c(SE_pred_rtmedian_outsamp, (df_observ_outsamp$rt_median[c] - preds$rt_median[c])^2)
    SE_pred_rtperc25_outsamp = c(SE_pred_rtperc25_outsamp, (df_observ_outsamp$rt_perc_25[c] - preds$rt_perc_25[c])^2)
    SE_pred_rtperc75_outsamp = c(SE_pred_rtperc75_outsamp, (df_observ_outsamp$rt_perc_75[c] - preds$rt_perc_75[c])^2)
    SE_pred_acc_outsamp = c(SE_pred_acc_outsamp, (df_observ_outsamp$acc[c] - preds$acc[c])^2)
    
    SE_mean_rtmean_insamp = c(SE_mean_rtmean_insamp, (df_observ_insamp$rt_mean[c] - rtmean_mean_insamp)^2)
    SE_mean_rtmedian_insamp = c(SE_mean_rtmedian_insamp, (df_observ_insamp$rt_median[c] - rtmedian_mean_insamp)^2)
    SE_mean_rtperc25_insamp = c(SE_mean_rtperc25_insamp, (df_observ_insamp$rt_perc_25[c] - rtperc25_mean_insamp)^2)
    SE_mean_rtperc75_insamp = c(SE_mean_rtperc75_insamp, (df_observ_insamp$rt_perc_75[c] - rtperc75_mean_insamp)^2)
    SE_mean_acc_insamp = c(SE_mean_acc_insamp, (df_observ_insamp$acc[c] - acc_mean_insamp)^2)
    SE_mean_rtmean_outsamp = c(SE_mean_rtmean_outsamp, (df_observ_outsamp$rt_mean[c] - rtmean_mean_outsamp)^2)
    SE_mean_rtmedian_outsamp = c(SE_mean_rtmedian_outsamp, (df_observ_outsamp$rt_median[c] - rtmedian_mean_outsamp)^2)
    SE_mean_rtperc25_outsamp = c(SE_mean_rtperc25_outsamp, (df_observ_outsamp$rt_perc_25[c] - rtperc25_mean_outsamp)^2)
    SE_mean_rtperc75_outsamp = c(SE_mean_rtperc75_outsamp, (df_observ_outsamp$rt_perc_75[c] - rtperc75_mean_outsamp)^2)
    SE_mean_acc_outsamp = c(SE_mean_acc_outsamp, (df_observ_outsamp$acc[c] - acc_mean_outsamp)^2)
  }
  
}

#get Rsq for insample prediction 
Rsq_insamppred_rtmean = 1-(sum(unlist(SE_pred_rtmean_insamp))/sum(unlist(SE_mean_rtmean_insamp)))
Rsq_insamppred_rtmedian = 1-(sum(unlist(SE_pred_rtmedian_insamp))/sum(unlist(SE_mean_rtmedian_insamp)))
Rsq_insamppred_rtperc25 = 1-(sum(unlist(SE_pred_rtperc25_insamp ))/sum(unlist(SE_mean_rtperc25_insamp)))
Rsq_insamppred_rtperc75 = 1-(sum(unlist(SE_pred_rtperc75_insamp ))/sum(unlist(SE_mean_rtperc75_insamp)))
Rsq_insamppred_acc = 1-(sum(unlist(SE_pred_acc_insamp ))/sum(unlist(SE_mean_acc_insamp)))

#get Rsq for outsample prediction 
Rsq_outsamppred_rtmean = 1-sum(unlist(SE_pred_rtmean_outsamp))/sum(unlist(SE_mean_rtmean_outsamp))
Rsq_outsamppred_rtmedian = 1-sum(unlist(SE_pred_rtmedian_outsamp))/sum(unlist(SE_mean_rtmedian_outsamp))
Rsq_outsamppred_rtperc25 = 1-sum(unlist(SE_pred_rtperc25_outsamp ))/sum(unlist(SE_mean_rtperc25_outsamp))
Rsq_outsamppred_rtperc75 = 1-sum(unlist(SE_pred_rtperc75_outsamp ))/sum(unlist(SE_mean_rtperc75_outsamp))
Rsq_outsamppred_acc = 1-sum(unlist(SE_pred_acc_outsamp ))/sum(unlist(SE_mean_acc_outsamp))


