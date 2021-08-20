# clear all variables 
rm(list=ls())
# load packages needed
library("DEoptim", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("ggthemes", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("ggplot2", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("gridExtra", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("reshape2", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("tictoc", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library('stringr')
# Add some dependent functions in other scripts
source("fitroutines.R")
source("helperfunctions.R")

jobnumbers = c(1,2,3,4,5)  #this is how many iterations the same model fitting will do with different random seeds to make sure parameter estimation is consistent -- default: the model will be fit to the data 5 times

for (currjobnum in jobnumbers) #for each round of model fitting 
{
  jobnum=currjobnum
  
  # -----------MODIFY MAIN SETTINGS HERE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  datafile = 'BT120819GPPost'       #Must be a "dat" file processed by csv_to_dat.R 
  model= "DDM"                    #Specify which model to fit: i.e. DDM, UGM, etc
  fullrangefit = 1                 #For specifying the constraint intervals/range of each parameter, use the full range recommended by CHaRTr? 1 = yes; 0 = no, I want to specify my own ranges in fitroutines.R
  
  #Are you modeling a parameter recovery data (generated data)? 
  ParamRecFile = 0      #0 for no, 1 for yes 
  model_simulated = 'DDM'  #if yes, then specify which model you simulated the parameter recovery datafile with  
  
  #Are you modeling the data from a single experimental session? 1 = yes, only single session; 0 = no, I'm fitting pooled data across experimental sessions
  single_session = 1
  
  #Names of the paramters and the order of parameters in the vector 
  #in names... 1 - drift rates (v1-vN for N number of drift rates or coherences), 2 - start point (z), 3 - upper bound (aU), 4 - non-decision time (Ter), 5 - parameter for collapsing bound (lambda) in cDDM
  #6 - parameter for collapsing bound (k) in cDDM, 7 - parameter for collapsing bound (aprime) in cDDM, 8 - parameter for urgency signal (intercept) in UGM variants, 9 - parameter for urgency signal (usign_var) for UGM variants
  names = c('v1','v2','v3', 'v4','v5','v6','v7','v8','v9','v10','v11', 'z','aU','Ter','lambda','k','aprime','intercept','usign_var') 
  
  #Are there fixed params? 0 for no, 1 for yes, in the order corresponding to "names" vector above 
  fixed_params = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0)
  
  ## Specifying the value of any fixed values if you're fixing any parameters (defined by fixed_params) -- the order of the parameters in this vector follows "names" vector. if fixed_params = 0 on a particular index, the value of the index in this vector won't matter since you're not fixing them
  fixed_values = c(-6.120330e+00,-5.357206e+00,-3.806973e+00,-3.444725e+00,-3.548213e+00,-2.744311e+00, -3.116225e+00,-1.779184e+00 , -1.894910e-01,1.662550e+00,3.693622e+00 , 4.894953e-01, 1.087356e+04 , 5.111463e-02 , 0, 0, 0, 0, 0.98482311)  #for monkey S post-injection parameter estimates for mUGM aU
  #fixed_values = c(-7.924324e+00, -8.280229e+00, -6.143492e+00, -3.634831e+00, -3.688837e+00, -3.304528e+00 , -1.543795e+00 , -9.410696e-01, 3.268737e-01, 3.726943e+00 , 4.603061e+00, 5.059867e-01, 1.288811e+04, 5.328652e-02, 0, 0, 0, 0, 0.98333457) #for monkey B pre-injection parameter estimates for mUGM aU

  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  fixed_params_df = data.frame(names,fixed_params,fixed_values)
  
  
  # getting the filename of the data to load (filename) and creating the name of the output file (saveFileName)
  if (ParamRecFile){
    filename = paste(paste(datafile,'_dat',sep=""),model_simulated, fakeparsnum, sep=" ")
    saveFileName=paste(datafile,model,jobnum)
  } else if (single_session) {
    filename = paste((paste(datafile,'_dat',sep = ""))) 
    saveFileName = paste(paste(datafile,'_IF_1',sep=""),model,jobnum)
  } else {
    filename = paste(datafile,'_dat',sep="")    
    saveFileName=paste(datafile,model,jobnum)
  }
  
  
  #If there are some fixed parameters, get the names of the parameters you want to fix 
  fixed_params_names = character()
  if  (sum(fixed_params_df['fixed_params']) != 0) {    
    for (n in 1:sum(fixed_params_df['fixed_params'])) {
      indices_of_1 = which(fixed_params_df['fixed_params'] ==1)
      fixed_params_names = paste(fixed_params_names,fixed_params_df[['names']][[indices_of_1[n]]],sep="")
    }
    saveFileName = paste(saveFileName,fixed_params_names,sep=" ")   #the saving file name will have the name of the parameters you have fixed
  }
  
  # Which job number to run, this tells you the run number
  set.seed(jobnum)    #set the random seed to the model fitting iteration number 
  fnam = letters[jobnum];    # letter index for each fit; a = 1st fit iteration, b = 2nd fit iteration, c = 3rd fit iteration... etc; for later model fit comparison measures
  
  # Setup some baseline parameters for fitting 
  contp = list(p=0)  # if estcontp==FALSE, give proportion of contaminant responses
  maxits = 750  # number of iterations of DEoptim to run
  nparticles = 400  # number of particles/chains for DEoptim
  nmc =10000  # number of MC samples for simulating the models at each iteration
  estcontp=FALSE  # estimate contaminant mixture probability from data. usually set to false
  bailouttime=4  # time after which to bail out of diffusion simulation of a trial, in seconds
  maxTimeStep=as.double(length(seq(0,bailouttime,.001)))   # max iterations for fitting routine, then bail out
  pred=F  # True if you want to generate model predictions from fitting routine (only use this once you have estimtated parameters from data using DEoptim)
  #nreps = 5;
  
  # cutoff for very slow RTs
  # gub = global upper bound in seconds - no RTs should be slower than this (and none were in data)
  gub=4
  
  #load the data
  load(filename)
  
  # for simple switching in fitting routine
  qps=as.numeric(dimnames(dat$q)[[1]]) 
  ncohs=1:dim(dat$q)[3]    #coherence index  
  # make parameter vector - different scaling for Stone and Stone+UGM
  nds=length(ncohs)        #number of signed coherences 
  actualParams = paramsandlims(model, nds,fakePars=TRUE, fullrangefit=fullrangefit,fixed_params_df=fixed_params_df)  #parameter dataframe for upper and lower constraints, fake parameters, and fitUGM flag
  fitUGM = unname(actualParams$fitUGM)   #fitUGM flag; if fitUGM is negative, then we're not fitting UGM or any of its variants
  lowers = actualParams$lowers           #lower constraints of each parameter
  uppers = actualParams$uppers           #upper constraits of each parameter 
  parnames = actualParams$parnames       #list of parameter names 
  stepsize=ifelse(fitUGM < 0,.001,1)  # time step for diffusion process: .001s (Stone, DDM), 1ms (UGM)
  stoch.s=ifelse(fitUGM < 0,.1,100)   # diffusion noise constant
  timecons=ifelse(fitUGM < 0,0,100)   # time constant for low pass filter if fitting UGM/UGM variants
  usign=ifelse(fitUGM < 0,0,1)        # scaling value for linear urgency function (usign=1 is linear with time) if fitting UGM/UGM variants. Evidence accumulation models needs usign=0, timecons=0
  
  
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


#Changes I made from the original chartr-demofit.R, chartr-HelperFcuntions.R, and chartr-FitRoutines.R
# -made correct/incorrect to toIF/awayIF
# - made starting point,z, an explicit free parameter to vary in DDM, cDDM, UGM, and bUGM. the z parameter is the 
# % of the distance between the lower and upper bound.
# -made estcontp = TRUE to NOT estimate contaminant RT reponses, because we'll be modeling the muscimol data
# -unfixed the 0 coherence drift rate value so it's not fixed to 0 anymore
#     - nstart=1 in paramsandlims for toIF/awayIF, whereas in correct/incorrect, nstart=2
#     - zeroCoh=FALSE in makeparamslist
# -for monkey S and monkey B data, I excluded the 36 condition bc they doesn't have enough error trials
# -for modeling the UGM with free slope and intercept for the urgency signal, I used modified bUGM to have the 
# usign_var, so now it has both free intercept and slope
# -made a fullrangefit feature, where if fullrangefit =1, the lowers and uppers of the parameters are the full 
# range of the drift rates pre-set by chartr. If =0, then the lowers and uppers of each drift rates are customized 
# manually in helperfunctions.R, to a tighter range for faster fitting
# - Made a free/fixed controls for parameters. The fitting_masterscript.R, helperfunctions.R, and fitroutine.R scripts
#   go with this update. I also made mbUGM as a separate model where it
#   varies the intercept and slope, while bUGM is just the slope.
