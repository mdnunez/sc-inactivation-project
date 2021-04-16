# Get predictions for a data session which was already fitted for parameter estimates. And also get the 
# in-sample and out-sample choice and RT stats too for Rsq calculation --> getRsq.R for insample and 
# out-sample prediction
# To get Rsq calculations once you have the dataset saved at the end of this script, run getRsq.R to
# obtain Rsq calculations and visual plot comparing prediction and in-sample data


rm(list=ls())
# load packages needed
library("DEoptim", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("ggthemes", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("ggplot2", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("gridExtra", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("reshape2", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("tictoc", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
# Add some helper functions
source("fitroutines_preds.R")
source("helperfunctions.R")

# -----------MODIFY MAIN SETTINGS HERE--------------------------------------------------------------------------
filename = paste('AllSPData_Post_IF', 'mUGM', '4' , 'aU') 
dataname = 'AllSPData_Post_IF'    # input name of the datafile 
model= "mUGM"                      #DDM vs cDDM vs UGM
jobnum=4 #choose the best fit from each of the models 
fullrangefit = 1;

#Get predictions for single sessions? 1 - yes, get predictions for single sessions, 0 - no, get predictions for pooled data
single_session = 0

#Are there fixed params? 0 for no, 1 for yes
#1-11 = drift rates, 12 = start point, 13 = bounds, 14 = Ter (non decision time), 15 = lambda parameter for collapsing bound,  16 = k parameter for collapsing bound, 17 = aprime parameter for collapsing bound, 18 = intercept parameter for urgency signal, 19 = slope parameter for urgency signal
fixed_params = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 1, 0, 0, 0, 0, 0, 0) # 0 = not a fixed param, 1 = fixed param
names = c('v1','v2','v3', 'v4','v5','v6','v7','v8','v9','v10','v11', 'z','aU','Ter','lambda','k','aprime','intercept','usign_var')
#Specify any fixed values if you want any values fixed
fixed_values = c(-6.120330e+00,-5.357206e+00,-3.806973e+00,-3.444725e+00,-3.548213e+00,-2.744311e+00, -3.116225e+00,-1.779184e+00 , -1.894910e-01,1.662550e+00,3.693622e+00 , 4.894953e-01, 1.087356e+04 , 5.111463e-02 , 0, 0, 0, 0, 0.98482311)  #from monkey S (SP) post-injection parameter estimates for mUGM aU

# --------------------------------------------------------------------------------------------------------------
saveFileName=paste(filename)
fixed_params_df = data.frame(names,fixed_params,fixed_values)

fixed_params_names = character()
if  (sum(fixed_params_df['fixed_params']) != 0) {    #If there are some fixed parameters
  for (n in 1:sum(fixed_params_df['fixed_params'])) {
    indices_of_1 = which(fixed_params_df['fixed_params'] ==1)
    fixed_params_names = paste(fixed_params_names,fixed_params_df[['names']][[indices_of_1[n]]],sep="")
  }
  saveFileName = paste(saveFileName,fixed_params_names,sep=" ")
}

# Which job number to run, this tells you the run number
set.seed(jobnum)
fnam = letters[jobnum];

# Setup some baseline parameters
contp = list(p=0)  # if estcontp==FALSE, give proportion of contaminant responses
maxits = 750  # number of iterations of DEoptim to run
nparticles = 400  # number of particles/chains for DEoptim
nmc =10000  # number of MC samples for simulating the models at each iteration
estcontp=FALSE  # estimate contaminant mixture probability from data. usually set to false
bailouttime=4  # time after which to bail out of diffusion sims, in seconds
maxTimeStep=as.double(length(seq(0,bailouttime,.001)))   # max iterations for fitting routine, then bail out
pred=T  # generate model predictions from fitting routine (only use this once you have estimtated parameters from data using DEoptim)
nreps = 5;
# cutoff for very slow RTs
# gub = global upper bound in seconds - no RTs should be slower than this (and none were in data)
gub=4


#load dat file
if (single_session == 1) {
  load(paste(paste('SCMuscimol_Study_datafiles/datafiles/dat files/for single session fitting/',dataname,'.csv',sep=""), 'IF_1'))
} else {
  load(paste('SCMuscimol_Study_datafiles/datafiles/dat files/for pooled session fitting/',dataname,'.csv',sep=""))
}

#load model fit data
load(paste("SCMuscimol_Study_datafiles/datafiles/out files/",filename, sep=""))

# for simple switching in fitting routine
qps=as.numeric(dimnames(dat$q)[[1]])
ncohs=1:dim(dat$q)[3]
nds=length(ncohs)
actualParams = paramsandlims(model, nds,fakePars=TRUE,fullrangefit=fullrangefit,fixed_params_df=fixed_params_df)
fitUGM = unname(actualParams$fitUGM)
lowers = actualParams$lowers
uppers = actualParams$uppers
parnames = actualParams$parnames
stepsize=ifelse(fitUGM < 0,.001,1)  # time step for diffusion process: .001s (Stone, DDM), 1ms (UGM)
stoch.s=ifelse(fitUGM < 0,.1,100)   # diffusion constant
timecons=ifelse(fitUGM < 0,0,100)   # time constant for low pass filter
usign=ifelse(fitUGM < 0,0,1)        # scaling value for linear urgency function (usign=1 is linear with time). EAM needs usign=0, timecons=0


params= out$pars
p=makeparamlist(params,fitUGM,ncohs, zeroCoh=FALSE,fixed_params_df) 

#get the predictions
preds=getpreds(x=out$pars,nmc,contp,ncohs,fitUGM,pred=F,gub,
               qps,stepsize,stoch.s,timecons,usign,parnames,maxTimeStep,fixed_params_df)


alltmprtlist = list()
alltmpresplist = list()
#Get the stats of the preds 
coherences = sort(unique(dat[['x']][['coh']]))
for (c in 1:dim(preds)[1]) {
  preds[c,'rt_mean'] = mean(unlist(preds[c,'rt']))
  preds[c,'rt_median'] = median(unlist(preds[c,'rt']))
  preds[c,'rt_perc_25'] = quantile(unlist(preds[c,'rt']),.25)
  preds[c,'rt_perc_75'] = quantile(unlist(preds[c,'rt']),.75)
  preds[c,'acc'] = sum(unlist(preds[c,'resp']) == 1)/length(unlist(preds[c,'resp']))
  
  alltmprtlist = c(alltmprtlist, list(dat$x['RT'][dat$x['coh'] == coherences[c]]))
  alltmpresplist = c(alltmpresplist, list(dat$x['response'][dat$x['coh'] == coherences[c]]))
}


#Get the observed stats for in-sample data if you haven't already
df_observ_insamp = data.frame(ncohs)
df_observ_insamp$rt = alltmprtlist
df_observ_insamp$resp = alltmpresplist
df_observ_insamp$ntrial = lengths(alltmpresplist)
for (c in 1:dim(df_observ_insamp)[1]) {
  df_observ_insamp[c,'rt_mean'] = mean(unlist(df_observ_insamp[c,'rt']))
  df_observ_insamp[c,'rt_median'] = median(unlist(df_observ_insamp[c,'rt']))
  df_observ_insamp[c,'rt_perc_25'] = quantile(unlist(df_observ_insamp[c,'rt']),.25)
  df_observ_insamp[c,'rt_perc_75'] = quantile(unlist(df_observ_insamp[c,'rt']),.75)
  df_observ_insamp[c,'acc'] = sum(unlist(df_observ_insamp[c,'resp']) == 1)/length(unlist(df_observ_insamp[c,'resp']))
}

#Get the observed stats for out-of-sample data if you haven't already
alltmprtlist = list()
alltmpresplist = list()
allntriallist = list()
load(paste(paste('SCMuscimol_Study_datafiles/datafiles/Trials files/Excluded for outsample/',dataname,'.csv',sep=""),'IF','leaveout'))   # directory of datafile for everything else

if (single_session == 1) {
  Trials_leaveout = Trials_leaveout[-(which(abs(Trials_leaveout['coh']) == 24)),]
} 

coherences = sort(unlist(unique(Trials_leaveout['coh'])))
for (c in 1:length(coherences)) {
  alltmprtlist = c(alltmprtlist, list(Trials_leaveout['RT'][Trials_leaveout['coh'] == coherences[c]]))
  alltmpresplist = c(alltmpresplist, list(Trials_leaveout['response'][Trials_leaveout['coh'] == coherences[c]]))
}    
df_observ_outsamp = data.frame(ncohs)
df_observ_outsamp$rt = alltmprtlist
df_observ_outsamp$resp = alltmpresplist
df_observ_outsamp$ntrial = lengths(alltmpresplist)
for (c in 1:length(coherences)) {
  df_observ_outsamp[c,'rt_mean'] = mean(unlist(df_observ_outsamp[c,'rt']))
  df_observ_outsamp[c,'rt_median'] = median(unlist(df_observ_outsamp[c,'rt']))
  df_observ_outsamp[c,'rt_perc_25'] = quantile(unlist(df_observ_outsamp[c,'rt']),.25)
  df_observ_outsamp[c,'rt_perc_75'] = quantile(unlist(df_observ_outsamp[c,'rt']),.75)
  df_observ_outsamp[c,'acc'] = sum(unlist(df_observ_outsamp[c,'resp']) == 1)/length(unlist(df_observ_outsamp[c,'resp']))
}


save(preds, file=paste(saveFileName,'_preds',sep=""))
save(df_observ_insamp, file=paste(saveFileName,'_insamp',sep=""))
save(df_observ_outsamp, file=paste(saveFileName,'_outsamp',sep=""))