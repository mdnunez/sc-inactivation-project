# Compute a set of variables from the the choice and RT data for each coherence.
# Note that in the original chartr scripts, they fit correct/incorrect choices, where correct = 1, and incorrect = 0.
# This makes an accumulate to bound model where one of the bound is correct choice, and the other bound is incorrect choice. 
# But if you want to fit a model with one bound being right choice whereas the other bound is the left choice, 
# then you can model that, but keep in mind that "accuracy" no longer means percent correct. It means the proportion 
# of choices to whatever the choice (i.e. right or to inactivated field choice) is = 1. -EJ

#loading some packages and setting working directory
library("DEoptim", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("ggthemes", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("ggplot2", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("gridExtra", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("reshape2", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("tictoc", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("R.matlab", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
# Remove variables
rm(list=ls())


# -----------MODIFY MAIN SETTINGS HERE--------------------------------------------------------------------------

datafile = 'BT120819GPPost'   # Name of the csv datafile without the ".csv"
path_of_datafile = paste('SCMuscimol_Study_datafiles/datafiles/csv files/awayIF_vs_toIF/Musc_rt/', sep="")    #specify the path of where the csv datafile is
remove_easy_coh =1             # 1 = remove the easy coherences; 0 = don't easy coherences; we want to remove the easy coherences when theres not much error trials for them
easy_coh = c(-24, 24, -36, 36)         #if you are removing easy coh, specify the coherences (both negative and positive coherences)
prop_leaveout = 0.2              # proportion of trials to leave out for out of sample prediction: input any proportion 0 to 1

# --------------------------------------------------------------------------------------------------------------

filename = paste(datafile,'.csv',sep="")    # input name of the datafile you want to process
Trials_all = read.csv(paste(path_of_datafile, '/',filename,sep=""));   # open csv

# to remove the easy conditions for modeling bc they don't have enough error trials
if (remove_easy_coh == 1) {
  index_easycoh = which(Trials_all$coh %in% easy_coh)
  Trials_all= Trials_all[-index_easycoh,]
}

# to leave out __% of trials for out of sample prediction
num_leaveout = prop_leaveout*lengths(Trials_all)[1]
leaveout_trials = round(runif(2*num_leaveout)*lengths(Trials_all)[1])
leaveout_trials = unique(leaveout_trials)[1:num_leaveout]
if (length(leaveout_trials) > 1) {
  Trials= Trials_all[-leaveout_trials,]
  Trials_leaveout = Trials_all[leaveout_trials,]
} else {
  Trials = Trials_all
}

# defining coherence, RT quantiles, and choice data to make quantile plots 
cohValues = unique(Trials$coh);
cohValues = sort(cohValues,decreasing=FALSE)
dNameRows = c("cor", "err");
dNameCols = cohValues;
quantileNames = c(seq(from=0.1,to=0.9,by=0.1),0.975)
quantileNames = c(seq(from=0.1,to=0.9,by=0.1))
counts = matrix(data=NA, nrow=2,ncol=length(cohValues), dimnames=list(dNameRows, dNameCols));
accuracy = c(cohValues);
quantiles = array(data=NA, dim = c(length(quantileNames),2,length(cohValues)), dimnames=list(quantileNames, dNameRows, dNameCols))
pb = array(data=NA, dim = c(10,2,length(cohValues)), dimnames=list(1:10, dNameRows, dNameCols))
plot(y=NULL,x=NULL,ylim=c(0, 1),xlim=c(0,1),ylab="RT",xlab="Proportion Correct")
whichQ = seq(1,10,2)


print("Writing down different coherence values")
for(cid in c(cohValues)){
  print(cid)
  # logical array of true/false for each trial  depending on current coh
  iX = Trials$coh==cid;
  temp = sum(Trials$response[iX]==1);   #how many are right/toIF/correct, depending on what 0 means for the response column in your csv file
  counts[1,cnt] = temp;
  temp = sum(Trials$response[iX]==0);  #how many are left/awayIF/incorrect, depending on what 0 means for response column in your csv file
  counts [2,cnt] = temp;
  accuracy[cnt] = counts[1,cnt]/(counts[1,cnt] + counts[2,cnt]);
  ixCorrect = Trials$response==1 & Trials$coh==cid;
  ixWrong = Trials$response==0 & Trials$coh==cid;
  quantiles[1:length(quantileNames),1,cnt] = quantile(Trials$RT[ixCorrect],quantileNames);
  quantiles[1:length(quantileNames),2,cnt] = quantile(Trials$RT[ixWrong],quantileNames);
  pb[1:10,1,cnt] = counts[1,cnt]/10; # Divide by 10 because you have quantiles :)!
  pb[1:10,2,cnt] = counts[2,cnt]/10;
  points(rep(accuracy[cnt],1,length(whichQ)), quantiles[whichQ,1,cnt],col="green3", pch=4)
  points(1-rep(accuracy[cnt],1,length(whichQ)), quantiles[whichQ,2,cnt],col="red3",pch=4)
  cnt = cnt + 1;
  # Get percentiles for the right/toIF/correct and left/awayIF/incorrect values for each value of coherence
}
names(accuracy)<-cohValues
p = accuracy;
q = quantiles;
n = counts;
x = subset(Trials,select=c("coh","response","RT"))
dat = list(p=p,n=n,q=q,pb=pb,x=x)   #the dat file you'll need for running the modeling script


# Saving file 
saveFileName = paste(datafile,'_dat',sep="")
if (prop_leaveout ==0) {     #if no trials were left out for out-of-sample prediction
  save(file=saveFileName,dat)   #just save the dat file
} else {                     #if trials were left out for out-of-sample prediction
  save(file=saveFileName,dat)   #save the dat file
  save(file=paste(saveFileName, '_leaveout',sep=""),Trials_leaveout)  #also saving the indices of the trials left out for out of sample prediction
}
