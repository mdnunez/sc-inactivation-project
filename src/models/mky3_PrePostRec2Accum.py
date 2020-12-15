# mky3_PrePostRec2Accum.py - Loads ALL behavioral data from Pre- Post- and Recovery-from Muscimol recordings 
#                    Model with all variables free across sessions, only modeling rightwardRT, with 2 explicit accumulators
#
# Copyright (C) 2020 Michael D. Nunez, <mdnunez@mednet.ucla.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Record of Revisions
#
# Date            Programmers                         Descriptions of Change
# ====         ================                       ======================
# 08/10/20      Michael Nunez          Converted from mky3_PrePostRec2Accum.py
# 08/11/20      Michael Nunez           Replace NAN with masked array

#Notes:
#Simple injection conditions: 3*2 = 6 (Pre, Post, Rec)*(Muscimol, Saline)
#Key:
# 1 - Pre, Muscimol
# 2 - Post, Muscimol
# 3 - Rec, Muscimol
# 4 - Pre, Saline
# 5 - Post, Saline
# 6 - Rec, Saline

import numpy as np
import scipy.io as sio
import scipy as sp
import os
import sys
from scipy.optimize import curve_fit
from time import strftime
import pyjags
import pandas as pd
import mky3_modelresults

#Flag for simplified injection conditions
simple = True

#Load data
if len(sys.argv) > 1:
    dataname = sys.argv[1]
else:
    dataname = 'fixedexpSCInjectGPRTdata'
data = sio.loadmat(f'../data/{dataname}.mat')

if simple:
    allrightwardRT = ((data['RF_choice']*2 - 3)*data['responsetime']).squeeze()
    rightwardRT = allrightwardRT[data['training'][0,:]==1]
    allinjectlarge = np.asarray(data['simpleinject'],dtype='int16').squeeze()
    injectlarge = allinjectlarge[data['training'][0,:]==1]
else:
    allrightwardRT = ((data['choice']*2 - 3)*data['responsetime']).squeeze()
    rightwardRT = allrightwardRT[data['training'][0,:]==1]
    allinjectlarge = np.asarray(data['injectioncond'],dtype='int16').squeeze()
    injectlarge = allinjectlarge[data['training'][0,:]==1]
allcondition = np.asarray(data['condition'],dtype='int16').squeeze()
condition = allcondition[data['training'][0,:]==1]
allsession = np.asarray(data['sessions'],dtype='int16').squeeze()
session = allsession[data['training'][0,:]==1]
allmonkeylarge = np.asarray(data['monkey'],dtype='int16').squeeze()
monkeylarge = allmonkeylarge[data['training'][0,:]==1]

N = np.size(rightwardRT)
nconditions = np.max(condition)
nsessions = np.max(session)

monkey = np.empty((nsessions))
injectfull = np.empty((nsessions))
maxrt = np.empty((nsessions))
minrt = np.empty((nsessions))
for s in np.unique(session):
    monkey[s-1] = monkeylarge[session == s][0]
    injectfull[s-1] = injectlarge[session == s][0]
    maxrt[s-1] = np.max(np.abs(rightwardRT[session == s]))
    minrt[s-1] = np.min(np.abs(rightwardRT[session == s]))


monkey = np.asarray(monkey, dtype='int16')
injectfull = np.asarray(injectfull, dtype='int16')
uniqueinject, inject = np.unique(injectfull,return_inverse=True) #Remove injection conditions without any samples
inject += 1
nmonkeys = np.max(monkey)
ninject = np.max(inject)

if simple:
    possible_labels = ['Pre Muscimol',
    'Post Muscimol',
    'Rec Muscimol',
    'Pre Saline',
    'Post Saline',
    'Rec Saline']
else:
    possible_labels = ['Pre, Right, Muscimol',
    'Post, Right, Muscimol',
    'Rec, Right, Muscimol',
    'Pre, Left, Muscimol',
    'Post, Left, Muscimol',
    'Rec, Left, Muscimol',
    'Pre, Right, Saline',
    'Post, Right, Saline',
    'Rec, Right, Saline',
    'Pre, Left, Saline',
    'Post, Left, Saline',
    'Rec, Left, Saline']
inject_labels = [possible_labels[i] for i in (uniqueinject-1)]



#Basic hierarchical drift-diffusion (DDM) and Psychometric models with lapse trials
thismodel = '''
model {

    #Between-session variability in non-decision time
    tersd ~ dgamma(.3,1)

    #Between-session variability in Speed-accuracy trade-off
    alphasd ~ dgamma(1,1)

    #Between-session variability in Rightward start point bias
    betasd ~ dgamma(.3,1)

    #Between-session variability in lapse trial probability
    problapsesd ~ dgamma(.3,1)

    #Gain effect of an injection experiment on probability of out IF response
    gainoutIF <- 1

    #Fix the gain effect of Pre-muscimol on probability of in IF response to 1
    #The constraint is for model identifiability
    gaininIF[1] <- 1

    for (e in 2:ninject) {
        #Gain effect of an injection experiment on probability of in IF response
        gaininIF[e] ~ dnorm(1, pow(3,-2))
    }

    for (e in 1:ninject) {
        for (m in 1:nmonkeys) {

            ##########
            #Hierarchical DDM parameters
            ##########

            #Hierarchical Non-decision time
            terhier[e,m] ~ dnorm(.5, pow(.25,-2))

            #Hierarchical Speed-accuracy trade-off
            alphahier[e,m] ~ dnorm(1, pow(.5,-2))

            #Hierarchical Rightward start point bias
            betahier[e,m] ~ dnorm(.5, pow(.25,-2))

            #Hierarchical lapse trial probability
            problapsehier[e,m] ~ dnorm(.3, pow(.15,-2))
        }
    }


    for (m in 1:nmonkeys) {
        #Number of seconds to accumulate an away from IF response when the probability of an IF response is 0%
        rateparam[m] ~ dnorm(.5, pow(.25, -2))T(.1, 4)
    } 

    for (c in 1:nconditions) {
        #Probability of in IF response as encoded by the IF SC
        SCinIF[c] ~ dnorm(.5, pow(.25, -2))T(0, 1)
    } 


    for (s in 1:nsessions) {

        ##########
        #DDM parameters
        ##########

        #Non-decision time
        ter[s] ~ dnorm(terhier[inject[s],monkey[s]], pow(tersd,-2))T(0, 1)

        #Speed-accuracy trade-off
        alpha[s] ~ dnorm(alphahier[inject[s],monkey[s]], pow(alphasd,-2))T(0, 3)

        #Rightward start point bias
        beta[s] ~ dnorm(betahier[inject[s],monkey[s]], pow(betasd,-2))T(0, 1)

        #Probability of a lapse trial
        problapse[s] ~ dnorm(problapsehier[inject[s],monkey[s]], pow(problapsesd,-2))T(0, 1)
        probDDM[s] <- 1 - problapse[s]

        for (c in 1:nconditions) {

            #Rightward/Leftward drift rate
            deltaright[s,c] <- ( SCinIF[c]*gaininIF[inject[s]] - (1 - SCinIF[c])*gainoutIF ) / ( rateparam[monkey[s]] )

        }

    }

    ##########
    # Wiener likelihoods
    for (i in 1:N) {

        # Log density for DDM process of rightward/leftward RT
        ld_comp2[i, 1] <- dlogwiener(rightwardRT[i], alpha[session[i]], ter[session[i]], beta[session[i]], deltaright[session[i],condition[i]])

        # Log density for lapse trials (negative max RT to positive max RT)
        ld_comp2[i, 2] <- logdensity.unif(rightwardRT[i], -3, 3)

        # Select one of these two densities (Mixture of nonlapse and lapse trials)
        selected_density[i] <- exp(ld_comp2[i, DDMorLapse[i]] - Constant)
        
        # Generate a likelihood for the MCMC sampler using a trick to maximize density value
        Ones2[i] ~ dbern(selected_density[i])

        # Probability of mind wandering trials (lapse trials)
        DDMorLapse[i] ~ dcat( c(probDDM[session[i]], problapse[session[i]]) )
    }
}
'''


# Input for mixture modeling
Ones = np.ones(N)
Constant = 20

# pyjags code

# Make sure $LD_LIBRARY_PATH sees /usr/local/lib
# Make sure that the correct JAGS/modules-4/ folder contains wiener.so and wiener.la
pyjags.modules.load_module('wiener')
pyjags.modules.load_module('dic')
pyjags.modules.list_modules()

nchains = 6
burnin = 2000  # Note that scientific notation breaks pyjags
nsamps = 50000

trackvars = ['rateparam', 'SCinIF', 'gaininIF',  
            'tersd', 'alphasd', 'betasd', 'problapsesd',
            'terhier', 'alphahier', 'betahier', 'problapsehier',
            'ter', 'alpha', 'beta', 'problapse',
             'deltaright']

# Initialize non-decision time with mininum RT across all conditions
# Use maximum RT for the bounds on the lapse process, modeled by a uniform distribution
initials = []
for c in range(0, nchains):
    chaininit = {
        'rateparam': np.random.uniform(.1 ,1., size=(nmonkeys)),
        'SCinIF': np.random.uniform(0., 1., size=(nconditions)),
        'gaininIF': np.random.uniform(0., 2., size=(ninject)),
        'tersd': np.random.uniform(.01, .5),
        'alphasd': np.random.uniform(.1, 1.),
        'betasd': np.random.uniform(.01, .5),
        'problapsesd': np.random.uniform(.01, .5),
        'terhier': np.random.uniform(0., .5, size=(ninject,nmonkeys)),
        'alphahier': np.random.uniform(.5, 2., size=(ninject,nmonkeys)),
        'betahier': np.random.uniform(.2, .8, size=(ninject,nmonkeys)),
        'problapsehier': np.random.uniform(.01, .1, size=(ninject,nmonkeys)),
        'ter': np.random.uniform(0., .5, size=(nsessions)),
        'alpha': np.random.uniform(.5, 2., size=(nsessions)),
        'beta': np.random.uniform(.2, .8, size=(nsessions)),
        'problapse': np.random.uniform(.01, .1,  size=(nsessions)),
    }
    for s in range(0, nsessions):
            chaininit['ter'][s] = np.random.uniform(0., minrt[s]/2)
    tempmask = np.zeros((ninject),dtype=bool)
    tempmask[0] = True
    chaininit['gaininIF'] = np.ma.array(chaininit['gaininIF'], mask=tempmask)
    initials.append(chaininit)

# Run JAGS model

savedir = '../modelfits/'
# Choose JAGS model type
modelname = 'PrePostRec2Accum'

# Save model
timestart = strftime('%b') + '_' + strftime('%d') + '_' + \
    strftime('%y') + '_' + strftime('%H') + '_' + strftime('%M')
filename = dataname + '_' + modelname + timestart
modelfile = filename + '.jags'
f = open(savedir+modelfile, 'w')
f.write(thismodel)
f.close()
print('Fitting model %s ...' % (modelfile))


indata = dict(N=N, rightwardRT=rightwardRT, condition=condition, session=session, monkey=monkey, inject=inject,
                                nconditions=nconditions, nsessions=nsessions, nmonkeys=nmonkeys, ninject=ninject,
                                Ones2=Ones, Constant=Constant)

threaded = pyjags.Model(file=savedir+modelfile, init=initials,
                        data=indata,
                        chains=nchains, adapt=burnin, threads=6,
                        progress_bar=True)


samples = threaded.sample(nsamps, vars=trackvars, thin=10)
samples['inject_labels'] = inject_labels

savestring = filename + ".mat"

print('Saving results to: \n %s' % (savestring))

sio.savemat(savedir + savestring, samples)

print('Running diagnostics and saving figures...\n')
(samples, diags, filename) = mky3_modelresults.load_data(filename)
mky3_modelresults.save_figures(samples, filename)

