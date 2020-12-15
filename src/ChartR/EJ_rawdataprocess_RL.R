# Compute a set of variables from the raw data.
# This is only for awayIF vs toIF organized csv files. The response variable in the csv files
# represents the side of the choice. 1 = chose toIF and 0 = chose awayIF

#loading some packages and setting working directory
library("DEoptim", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("ggthemes", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("ggplot2", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("gridExtra", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("reshape2", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("tictoc", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
#setwd("/CHaRTr-master_1")
# Remove variables
rm(list=ls())


# -----------MODIFY MAIN SETTINGS HERE--------------------------------------------------------------------------
filename = 'BT120819GPPre.csv'    # input name of the datafile you want to process
Trials = read.csv("../data/BT120819GPPre.csv");   # directory of datafile
toIF = 1                          # 1 = only analyze trials that were set toIF; 0 = only analyze trials that were set awayIF (not choice, but the tg direction)
remove_easy_coh = 1               # 1 = remove the 24 and 36 coherences; 0 = don't remove 24 and 36 coherences; we want to remove the easy coherences bc theres not much error trials for them
prop_leaveout = .2               # proportion of trials to leave out for out of sample prediction: input any proportion 0 to 1
# --------------------------------------------------------------------------------------------------------------


# to remove the 24 and 36 condition for correct vs incorrect modeling bc they don't have enough error trials
if (remove_easy_coh == 1) {
  index_24_36 = which(Trials$coh %in% c(24, 36))
  Trials= Trials[-index_24_36,]
}
# to only choose trials that are set awayIF or toIF and keep all 0 coh trials
if (toIF == 1) {
  index_toIF = which(Trials$tgdir %in% 1)
  index_0coh = which(Trials$coh %in% 0)
  index_tokeep = unique(c(index_toIF, index_0coh))
  Trials = Trials[index_tokeep,]
  saveFileName=paste(filename,'toIF')
} else if (toIF == 0) {
  index_awayIF = which(Trials$tgdir %in% 0)
  index_0coh = which(Trial$coh %in% 0)
  index_tokeep = unique(c(index_awayIF, index_0coh))
  Trials = Trials[index_awayIF,]
  saveFileName=paste(filename,'awayIF')
}

# to leave out __% of trials for out of sample prediction
num_leaveout = prop_leaveout*lengths(Trials)[1]
leaveout_trials = round(runif(num_leaveout)*lengths(Trials)[1])
if (length(leaveout_trials) > 0) {
  Trials= Trials[-leaveout_trials,]
}
  



cohValues = unique(Trials$coh);
cohValues = sort(cohValues,decreasing=FALSE)
# Now count number of trials for each condition
# Should be much easier than this but I am struggling :)
dNameRows = c("cor", "err");
dNameCols = cohValues;
quantileNames = c(seq(from=0.1,to=0.9,by=0.1),0.975)
quantileNames = c(seq(from=0.1,to=0.9,by=0.1))
counts = matrix(data=NA, nrow=2,ncol=length(cohValues), dimnames=list(dNameRows, dNameCols));
accuracy = c(cohValues);
quantiles = array(data=NA, dim = c(length(quantileNames),2,length(cohValues)), dimnames=list(quantileNames, dNameRows, dNameCols))
pb = array(data=NA, dim = c(10,2,length(cohValues)), dimnames=list(1:10, dNameRows, dNameCols))
plot(y=NULL,x=NULL,ylim=c(0.3, 1.5),xlim=c(0,1),ylab="RT",xlab="Proportion Correct")
whichQ = seq(1,10,2)

cnt = 1;


print("Writing down different coherence values")
# Here we had 7 levels of coherence
for(cid in c(cohValues)){
  print(cid)
  # logical array of true/false for each trial depending on current coh
  iX = Trials$coh==cid;
  #print(sum(Trials$response[iX]==1))
  temp = sum(Trials$response[iX]==1);   #how many are toIF
  counts[1,cnt] = temp;
  temp = sum(Trials$response[iX]==0);  #how many are awayIF
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
  # Get percentiles for the correct and incorrect values for each value of coherence
}
names(accuracy)<-cohValues
p = accuracy;
q = quantiles;
n = counts;
x = subset(Trials,select=c("coh","response","RT"));
dat = list(p=p,n=n,q=q,pb=pb,x=x);
save(file=saveFileName,dat);
# Now create the pb