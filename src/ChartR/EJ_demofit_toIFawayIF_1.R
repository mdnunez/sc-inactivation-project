# Remove variables
rm(list=ls())
# load packages needed
library("DEoptim", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("ggthemes", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("ggplot2", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("gridExtra", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("reshape2", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("tictoc", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library('stringr')
#library("R.matlab", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
matfilepath = '../data/'
# set working directory
#setwd("/CHaRTr-master_1")
# Add some helper functions
source("EJ_fitroutines_toIFawayIF_1.R")
source("EJ_helperfunctions_toIFawayIF_1.R")

jobnumbers = c(1,2,3,4,5)
for (currjobnum in jobnumbers)
{
  # -----------MODIFY MAIN SETTINGS HERE--------------------------------------------------------------------------
datafile = 'SP012120GPPre'       #AllBTData_Pre_IF for pooled BT data or BT_ParamRecSimData for param recovery
model= "mUGM"                      #Model to fit: DDM vs cDDM vs UGM
jobnum=currjobnum
fullrangefit = 1

#Is this a parameter recovery file?
ParamRecFile = 0 
fakeparsnum = 2
model_simulated = 'DDM'

#Is it single session?
single_session = 1

#Are there fixed params? 0 for no, 1 for yes
#1 - drift rates, 2 - start point, 3 - bounds, 4 - Ter, 5 - 
fixed_params = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 1, 0, 0, 0, 0, 0, 0)
names = c('v1','v2','v3', 'v4','v5','v6','v7','v8','v9','v10','v11', 'z','aU','Ter','lambda','k','aprime','intercept','usign_var')
#fixed_values = c(-5.82936467, -5.31630420, -3.17205091, -1.67667785, -0.36760022, 0.49816691, 1.40051163 , 1.74577969, 5.04039534, 6.15790066, 7.57566033, 0.49285957, 9730.84304149, 0.08761421, 0, 0, 0, 0, 0.98333457) #fron butters's pre-injection parameter estimates for mUGM aU
fixed_values = c(-6.120330e+00,-5.357206e+00,-3.806973e+00,-3.444725e+00,-3.548213e+00,-2.744311e+00, -3.116225e+00,-1.779184e+00 , -1.894910e-01,1.662550e+00,3.693622e+00 , 4.894953e-01, 1.087356e+04 , 5.111463e-02 , 0, 0, 0, 0, 0.98482311)  #fron sparky's post-injection parameter estimates for mUGM aU
fixed_values = c(-7.924324e+00, -8.280229e+00, -6.143492e+00, -3.634831e+00, -3.688837e+00, -3.304528e+00 , -1.543795e+00 , -9.410696e-01, 3.268737e-01, 3.726943e+00 , 4.603061e+00, 5.059867e-01, 1.288811e+04, 5.328652e-02, 0, 0, 0, 0, 0.98333457) #fron butters's pre-injection parameter estimates for mUGM aU



#fixed_values = c(-5.82936467, -5.31630420, -3.17205091, -1.67667785, -0.36760022, 0.49816691, 1.40051163 , 1.74577969, 5.04039534, 6.15790066, 7.57566033, 0.49285957, 9730.84304149, 0.08761421, 0, 0, 0, 0, 0.98333457) #fron butters's pre-injection parameter estimates for mUGM aU
#fixed_values = c(-5.33924206,-4.09048342,-2.68155036,-1.88956904,-1.75060943,-0.51704293, 0.97633100,0.92245143, 2.24383979,3.24078443,5.09557487, 0.51042833, 9786.12900470, 0.05054732, 0, 0, 0, 0, 0.98482311)  #fron sparky's pre-injection parameter estimates for mUGM aU
# --------------------------------------------------------------------------------------------------------------

fixed_params_df = data.frame(names,fixed_params,fixed_values)


if (ParamRecFile){
  filename = paste(paste(datafile,'.csv',sep=""),model_simulated, fakeparsnum, sep=" ")
  saveFileName=paste(filename,model,jobnum)
} else if (single_session) {
  filename = paste((paste(datafile,'.csv',sep = "")), 'IF_1') 
  saveFileName = paste(paste(datafile,'_IF_1',sep=""),model,jobnum)
} else {
  filename = paste(datafile,'.csv',sep="")    # input name of the datafile you want to process
  saveFileName=paste(datafile,model,jobnum)
}



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
estcontp=TRUE  # estimate contaminant mixture probability from data. usually set to false
bailouttime=4  # time after which to bail out of diffusion sims, in seconds
maxTimeStep=as.double(length(seq(0,bailouttime,.001)))   # max iterations for fitting routine, then bail out
pred=F  # generate model predictions from fitting routine (only use this once you have estimtated parameters from data using DEoptim)
nreps = 5;
# cutoff for very slow RTs
# gub = global upper bound in seconds - no RTs should be slower than this (and none were in data)
gub=4

#load the data
load(filename)

# for simple switching in fitting routine
qps=as.numeric(dimnames(dat$q)[[1]])
ncohs=1:dim(dat$q)[3]
# make parameter vector - different scaling for Stone and Stone+UGM
#    order: drifts (5), eta, upper boundary, Ter
nds=length(ncohs)
actualParams = paramsandlims(model, nds,fakePars=TRUE, fullrangefit=fullrangefit,fixed_params_df=fixed_params_df)
fitUGM = unname(actualParams$fitUGM)
lowers = actualParams$lowers
uppers = actualParams$uppers
parnames = actualParams$parnames
stepsize=ifelse(fitUGM < 0,.001,1)  # time step for diffusion process: .001s (Stone, DDM), 1ms (UGM)
stoch.s=ifelse(fitUGM < 0,.1,100)   # diffusion constant
timecons=ifelse(fitUGM < 0,0,100)   # time constant for low pass filter
usign=ifelse(fitUGM < 0,0,1)        # scaling value for linear urgency function (usign=1 is linear with time). EAM needs usign=0, timecons=0

params= actualParams$fakeParams
p=makeparamlist(params,fitUGM,ncohs, zeroCoh=FALSE,fixed_params_df) 

  
# DEoptim parameter estimation 
tic
print(paste("Starting optimization ...", saveFileName))
library(DEoptim)
system.time({
  tmp=DEoptim(
    fn=obj,
    lower=lowers,
    upper=uppers,
    dat=dat,
    nmc=nmc,
    contp=contp,
    ncohs=ncohs,
    fitUGM=fitUGM,
    gub=gub,
    pred=FALSE,
    qps=qps,
    stepsize=stepsize,
    stoch.s=stoch.s,
    timecons=timecons,
    usign=usign,
    parnames=parnames,
    maxTimeStep=maxTimeStep,
    fixed_params_df=fixed_params_df,
    control=DEoptim.control(itermax=maxits,NP=nparticles,trace=TRUE,
                            parallelType=1,reltol=1e-6,steptol=200,
                            # load objects used for fitting, for parallelType==1
                            parVar=list("dat","lowers","uppers","nmc","contp","ncohs", "fitUGM","pred",
                                        "qps", "stepsize","stoch.s","timecons","usign","parnames","maxTimeStep","maxits","nparticles","gub",
                                        "diffusionC","makeparamlist","contaminantmixresps","qmpouts","getpreds","obj","returnListOfModels")
                            # same again, but for functions
    ))})

cat(paste("dataset: ", filename, ", model:", model, fnam, "\n\n",sep=" "))
out=tmp$optim$bestmem
names(out)=parnames
print(round(out,4))
print(round(tmp$optim$bestval,4))

# re-calculate obj for best fitting parameters, to determine amount of noise
# in the obj value for the best fit
mcsforpreds=50000
reobj=obj(x=tmp$optim$bestmem,dat=dat,nmc=mcsforpreds,
          contp=contp,ncohs=ncohs,fitUGM=fitUGM,gub=gub,pred=FALSE,
          qps=qps,stepsize=stepsize,stoch.s=stoch.s,timecons=timecons,usign=usign,
          parnames=parnames,maxTimeStep=maxTimeStep,fixed_params_df=fixed_params_df)
print(round(reobj,4))
# Now compute it for each level of quantile. Suspicion is that you go awry for the hardest coherences and you really need to think
# about what goes on there. Life is not easy there :)
mcsforpreds=50000
reobjperpoint=objPerQ(x=tmp$optim$bestmem,dat=dat,nmc=mcsforpreds,
                      contp=contp,ncohs=ncohs,fitUGM=fitUGM,gub=gub,pred=FALSE,
                      qps=qps,stepsize=stepsize,stoch.s=stoch.s,timecons=timecons,usign=usign,
                      parnames=parnames,maxTimeStep=maxTimeStep,fixed_params_df=fixed_params_df)
print(round(reobjperpoint,4))
out=list(dataset=filename,model=model,ndataset=fnam,pars=out,
         obj=-tmp$optim$bestval,reobj=-reobj, reobjperpoint=reobjperpoint)

out$fixed_params = fixed_params_df
toc()
save(out,file=saveFileName)
}


#Changes I made for making toIF vs awayIF
# -made correct/incorrect to toIF/awayIF
# - made starting point,z, an explicit free parameter to vary in DDM, cDDM, UGM, and bUGM. the z parameter is the 
# % of the distance between the lower and upper bound.
# -made estcontp = TRUE to NOT estimate contaminant RT reponses, because we'll be modeling the muscimol 
# -unfixed the 0 coherence drift rate value so it's not fixed to 0 anymore
#     - nstart=1 in paramsandlims for toIF/awayIF, whereas in correct/incorrect, nstart=2
#     - zeroCoh=FALSE in makeparamslist
# -for SP and BT data, I excluded the 36 condition bc rhe doesn't have enough trials for that 
# -for modeling the UGM with free slope and intercept for the urgency signal, I used modified bUGM to have the 
# usign_var, so now it has both free intercept and slope
# -made a fullrangefit feature, where if fullrangefit =1, the lowers and uppers of the parameters are the full 
# range of the drift rates pre-set by chartr. If =0, then the lowers and uppers of each drift rates are customized 
# manually in EJ_helperfunctions_toIFawayIF.R, to a tighter range for faster fitting
# - 4/1/20 Made a free/fixed controls for parameters. The demofit, helperfunctions, and fitroutine scripts
#   that go with this update has a "_1" at the end of the name. I also made mbUGM as a separate model where it
#   varies the intercept and slope, while bUGM is just the slope.
