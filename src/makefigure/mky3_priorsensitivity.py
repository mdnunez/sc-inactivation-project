# mky3_priorsensitivity.py - Generates a table of full HDDM results with different priors
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
# 11/18/20     Michael Nunez                              Original code

# Imports
import numpy as np
import scipy.io as sio
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats.stats import pearsonr
import os
import sys
# from IPython import get_ipython  # Run magic functions from script
# get_ipython().magic('pylab')  # Initialize ipython matplotlib plotting graphics

def diagnostic(insamples):
    """
    Returns Rhat (measure of convergence, less is better with an approximate
    1.10 cutoff) and Neff, number of effective samples).

    Reference: Gelman, A., Carlin, J., Stern, H., & Rubin D., (2004).
              Bayesian Data Analysis (Second Edition). Chapman & Hall/CRC:
              Boca Raton, FL.


    Parameters
    ----------
    insamples: dic
        Sampled values of monitored variables as a dictionary where keys
        are variable names and values are numpy arrays with shape:
        (dim_1, dim_n, iterations, chains). dim_1, ..., dim_n describe the
        shape of variable in JAGS model.

    Returns
    -------
    dict:
        Rhat for each variable. Prints Maximum Rhat
    """

    result = {}  # Initialize dictionary
    maxrhats = np.zeros((len(insamples.keys())), dtype=float)
    maxrhatsnew = np.zeros((len(insamples.keys())), dtype=float)
    minneff = np.ones((len(insamples.keys())), dtype=float)*np.inf
    allkeys ={} # Initialize dictionary
    keyindx = 0
    for key in insamples.keys():
        if key[0] != '_':
            result[key] = {}
            
            possamps = insamples[key]
            
            # Number of chains
            nchains = possamps.shape[-1]
            
            # Number of samples per chain
            nsamps = possamps.shape[-2]
            
            # Number of variables per key
            nvars = np.prod(possamps.shape[0:-2])
            
            # Reshape data
            allsamps = np.reshape(possamps, possamps.shape[:-2] + (nchains * nsamps,))

            # Reshape data to preduce R_hatnew
            possampsnew = np.empty(possamps.shape[:-2] + (int(nsamps/2), nchains * 2,))
            newc=0
            for c in range(nchains):
                possampsnew[...,newc] = np.take(np.take(possamps,np.arange(0,int(nsamps/2)),axis=-2),c,axis=-1)
                possampsnew[...,newc+1] = np.take(np.take(possamps,np.arange(int(nsamps/2),nsamps),axis=-2),c,axis=-1)
                newc += 2

            # Index of variables
            varindx = np.arange(nvars).reshape(possamps.shape[0:-2])
            
            # Reshape data
            alldata = np.reshape(possamps, (nvars, nsamps, nchains))
                    
            # Mean of each chain for rhat
            chainmeans = np.mean(possamps, axis=-2)
            # Mean of each chain for rhatnew
            chainmeansnew = np.mean(possampsnew, axis=-2)
            # Global mean of each parameter for rhat
            globalmean = np.mean(chainmeans, axis=-1)
            globalmeannew = np.mean(chainmeansnew, axis=-1)
            result[key]['mean'] = globalmean
            result[key]['std'] = np.std(allsamps, axis=-1)
            globalmeanext = np.expand_dims(
                globalmean, axis=-1)  # Expand the last dimension
            globalmeanext = np.repeat(
                globalmeanext, nchains, axis=-1)  # For differencing
            globalmeanextnew = np.expand_dims(
                globalmeannew, axis=-1)  # Expand the last dimension
            globalmeanextnew = np.repeat(
                globalmeanextnew, nchains*2, axis=-1)  # For differencing
            # Between-chain variance for rhat
            between = np.sum(np.square(chainmeans - globalmeanext),
                             axis=-1) * nsamps / (nchains - 1.)
            # Mean of the variances of each chain for rhat
            within = np.mean(np.var(possamps, axis=-2), axis=-1)
            # Total estimated variance for rhat
            totalestvar = (1. - (1. / nsamps)) * \
                within + (1. / nsamps) * between
            # Rhat (original Gelman-Rubin statistic)
            temprhat = np.sqrt(totalestvar / within)
            maxrhats[keyindx] = np.nanmax(temprhat) # Ignore NANs
            allkeys[keyindx] = key
            result[key]['rhat'] = temprhat
            # Between-chain variance for rhatnew
            betweennew = np.sum(np.square(chainmeansnew - globalmeanextnew),
                             axis=-1) * (nsamps/2) / ((nchains*2) - 1.)
            # Mean of the variances of each chain for rhatnew
            withinnew = np.mean(np.var(possampsnew, axis=-2), axis=-1)
            # Total estimated variance
            totalestvarnew = (1. - (1. / (nsamps/2))) * \
                withinnew + (1. / (nsamps/2)) * betweennew
            # Rhatnew (Gelman-Rubin statistic from Gelman et al., 2013)
            temprhatnew = np.sqrt(totalestvarnew / withinnew)
            maxrhatsnew[keyindx] = np.nanmax(temprhatnew) # Ignore NANs
            result[key]['rhatnew'] = temprhatnew
            # Number of effective samples from Gelman et al. (2013) 286-288
            neff = np.empty(possamps.shape[0:-2])
            for v in range(0, nvars):
                whereis = np.where(varindx == v)
                rho_hat = []
                rho_hat_even = 0
                rho_hat_odd = 0
                t = 2
                while ((t < nsamps - 2) & (float(rho_hat_even) + float(rho_hat_odd) >= 0)):
                    variogram_odd = np.mean(np.mean(np.power(alldata[v,(t-1):nsamps,:] - alldata[v,0:(nsamps-t+1),:],2),axis=0)) # above equation (11.7) in Gelman et al., 2013
                    rho_hat_odd = 1 - np.divide(variogram_odd, 2*totalestvar[whereis]) # Equation (11.7) in Gelman et al., 2013
                    rho_hat.append(rho_hat_odd)
                    variogram_even = np.mean(np.mean(np.power(alldata[v,t:nsamps,:] - alldata[v,0:(nsamps-t),:],2),axis=0)) # above equation (11.7) in Gelman et al., 2013
                    rho_hat_even = 1 - np.divide(variogram_even, 2*totalestvar[whereis]) # Equation (11.7) in Gelman et al., 2013
                    rho_hat.append(rho_hat_even)
                    t += 2
                rho_hat = np.asarray(rho_hat)
                neff[whereis] = np.divide(nchains*nsamps, 1 + 2*np.sum(rho_hat)) # Equation (11.8) in Gelman et al., 2013
            result[key]['neff'] = np.round(neff) 
            minneff[keyindx] = np.nanmin(np.round(neff))
            keyindx += 1

            # Geweke statistic?
    # print("Maximum Rhat was %3.2f for variable %s" % (np.max(maxrhats),allkeys[np.argmax(maxrhats)]))
    maxmaxrhat = np.max(maxrhatsnew)
    print("Maximum Rhatnew was %3.2f for variable %s" % (maxmaxrhat,allkeys[np.argmax(maxrhatsnew)]))
    minminneff = np.min(minneff)
    print("Minimum number of effective samples was %d for variable %s" % (minminneff,allkeys[np.argmin(minneff)]))
    return result, maxrhatsnew, minminneff



def load_data(filename):
    # Load model results
    modelloc = (f'../modelfits/{filename}.mat')
    samples = sio.loadmat(modelloc)
    samples_diagrelevant = samples.copy()
    samples_diagrelevant.pop('DDMorLapse', None) #Remove variable DDMorLapse to obtain Rhat diagnostics
    samples_diagrelevant.pop('inject_labels',None)
    diags = diagnostic(samples_diagrelevant)

    nchains = samples['ter'].shape[-1]
    nsamplesperchain = samples['ter'].shape[-2]
    nconditions = samples['ter'].shape[0]
    return (samples, diags, filename)


saveloc = '../../figures/fixedexpSCInjectGPRTdata/PriorSensitivityTable.xlsx'
# Load full HDDM fits
original = 'fixedexpSCInjectGPRTdata_PrePostRecHalfApr_10_20_15_50'
shifted = 'fixedexpSCInjectGPRTdata_PrePostRecHalfShiftedNov_02_20_13_48'
narrow = 'fixedexpSCInjectGPRTdata_PrePostRecHalfNarrowNov_02_20_11_37'
wide = 'fixedexpSCInjectGPRTdata_PrePostRecHalfWideNov_04_20_13_41'

tablemodels = [original, shifted, narrow, wide]

modelnames = ['Original', 'Shifted', 'Narrow', 'Wide']
paramnames = ['Prob. of mean drift rate decrease in monkey B', 'Prob. of mean drift rate decrease in monkey S', 'BF for a mean drift rate not equal to zero in monkey B',
     'BF for a mean drift rate not equal to zero in monkey S', 'Prob. of start point decrease in monkey B', 'Prob. of start point decrease in monkey S',
     'BF for a start point equal to .5 in monkey B', 'BF for a start point equal to .5 in monkey S', 'Prob. of an increase in non-decision time in monkey B', 'Prob. of an increase in non-decision time in monkey S',
     'Prob. of an increase in symmetric boundary in monkey B', 'Prob. of an increase in symmetric boundary in monkey S', 'Prob. of an increase in lapse proportion in monkey B', 'Prob. of an increase in lapse proportion in monkey S']
tabledata = np.empty((len(paramnames),len(modelnames)))


driftpriormeans = np.array([0, 1, 0, 0])
driftpriorstd = np.array([2, 2, 1, 4])

betapriormeans = np.array([.5, .25, .5, .5])
betapriorstd = np.array([.25, .25, .125, .5])

stepsize = .001 #Step size to calculate probabilities from differentiating approximate posterior distributions

modeltrack = 0
for model in tablemodels:
    (samples, diags, filename) = load_data(filename=model)
    deltarightsamps = samples['deltarighthier']
    betasamps = samples['betahier']
    tersamps = samples['terhier']
    alphasamps = samples['alphahier']
    problapsesamps = samples['problapsehier']*100

    nchains = deltarightsamps.shape[-1]
    nsamps = deltarightsamps.shape[-2]
    ninject = deltarightsamps.shape[0]
    nmonkey = deltarightsamps.shape[1]

    
    #Computate probabilities and Bayes factors from posterior distributions of drift rate
    reshapedsamps = np.reshape(deltarightsamps, (ninject, nmonkey, nsamps*nchains))
    diff_pre_post_BT = reshapedsamps[1,1,:] - reshapedsamps[0,1,:]
    kde_diff_BT = stats.gaussian_kde(diff_pre_post_BT)
    wherecalc_neg_BT = np.arange(np.min(diff_pre_post_BT[diff_pre_post_BT <= 0]),0,step=stepsize)
    prob_muscimol_negchange_BT = np.sum(kde_diff_BT(wherecalc_neg_BT)*stepsize)
    print("The probability of a negative change in deltaright for monkey B was %3.2f%%" % (prob_muscimol_negchange_BT*100))
    tabledata[0,modeltrack] = (prob_muscimol_negchange_BT*100)
    diff_pre_post_SP = reshapedsamps[1,0,:] - reshapedsamps[0,0,:]
    kde_diff_SP = stats.gaussian_kde(diff_pre_post_SP)
    wherecalc_neg_SP = np.arange(np.min(diff_pre_post_SP[diff_pre_post_SP <= 0]),0,step=stepsize)
    prob_muscimol_negchange_SP = np.sum(kde_diff_SP(wherecalc_neg_SP)*stepsize)
    print("The probability of a negative change in deltaright for monkey S was %3.2f%%" % (prob_muscimol_negchange_SP*100))
    tabledata[1,modeltrack] = (prob_muscimol_negchange_SP*100)
    # Calculate Bayes Factors using the Savage-Dickey density ratio
    denom = stats.norm.pdf(0, loc=driftpriormeans[modeltrack], scale=driftpriorstd[modeltrack]) #Matches the Hierarchical Rightward/Leftward drift rate prior in each models' JAGS code
    kde_post_BT = stats.gaussian_kde(reshapedsamps[1,1,:])
    num_post_BT = kde_post_BT(0)
    bf_post_BT = num_post_BT / denom
    print("The Bayes Factor for a mean drift of NOT 0 for monkey B in Post-data was %3.4f" % (1/bf_post_BT))
    tabledata[2,modeltrack] = (1/bf_post_BT)
    kde_post_SP = stats.gaussian_kde(reshapedsamps[1,0,:])
    num_post_SP = kde_post_SP(0)
    bf_post_SP = num_post_SP / denom
    print("The Bayes Factor for a mean drift of NOT 0 for monkey S in Post-data was %3.4f" % (1/bf_post_SP))
    tabledata[3,modeltrack] = (1/bf_post_SP)

    #Computate probabilities and Bayes factors from posterior distribution of start point
    reshapedsamps = np.reshape(betasamps, (ninject, nmonkey, nsamps*nchains))
    diff_pre_post_BT = reshapedsamps[1,1,:] - reshapedsamps[0,1,:]
    kde_diff_BT = stats.gaussian_kde(diff_pre_post_BT)
    wherecalc_neg_BT = np.arange(np.min(diff_pre_post_BT[diff_pre_post_BT <= 0]),0,step=stepsize)
    prob_muscimol_negchange_BT = np.sum(kde_diff_BT(wherecalc_neg_BT)*stepsize)
    print("The probability of a negative change in start point for monkey B was %3.2f%%" % (prob_muscimol_negchange_BT*100))
    tabledata[4,modeltrack] = (prob_muscimol_negchange_BT*100)
    diff_pre_post_SP = reshapedsamps[1,0,:] - reshapedsamps[0,0,:]
    kde_diff_SP = stats.gaussian_kde(diff_pre_post_SP)
    wherecalc_neg_SP = np.arange(np.min(diff_pre_post_SP[diff_pre_post_SP <= 0]),0,step=stepsize)
    prob_muscimol_negchange_SP = np.sum(kde_diff_SP(wherecalc_neg_SP)*stepsize)
    print("The probability of a negative change in start point for monkey S was %3.2f%%" % (prob_muscimol_negchange_SP*100))
    tabledata[5,modeltrack] = (prob_muscimol_negchange_SP*100)
    # Calculate Bayes Factors using the Savage-Dickey density ratio
    denom = stats.norm.pdf(.5, loc=betapriormeans[modeltrack], scale=betapriorstd[modeltrack]) #Matches the Hierarchical Rightward start point prior in each models' JAGS code
    kde_post_BT = stats.gaussian_kde(reshapedsamps[1,1,:])
    num_post_BT = kde_post_BT(.5)
    bf_post_BT = num_post_BT / denom
    print("The Bayes Factor for a start point of .5 for monkey B in Post-data was %3.4f" % (bf_post_BT))
    tabledata[6,modeltrack] = (bf_post_BT)
    kde_post_SP = stats.gaussian_kde(reshapedsamps[1,0,:])
    num_post_SP = kde_post_SP(.5)
    bf_post_SP = num_post_SP / denom
    print("The Bayes Factor for a start point of .5 for monkey S in Post-data was %3.4f" % (bf_post_SP))
    tabledata[7,modeltrack] = (bf_post_SP)

    # Calculate the probability of an increase in non-decision time
    reshapedsamps = np.reshape(tersamps, (ninject, nmonkey, nsamps*nchains))
    diff_pre_post_BT = reshapedsamps[1,1,:] - reshapedsamps[0,1,:]
    kde_diff_BT = stats.gaussian_kde(diff_pre_post_BT)
    wherecalc_pos_BT = np.arange(0, np.max(diff_pre_post_BT[diff_pre_post_BT >= 0]),step=stepsize)
    prob_muscimol_poschange_BT = np.sum(kde_diff_BT(wherecalc_pos_BT)*stepsize)
    print("The probability of a positive change in NDT for monkey B was %3.2f%%" % (prob_muscimol_poschange_BT*100))
    tabledata[8,modeltrack] = prob_muscimol_poschange_BT*100
    diff_pre_post_SP = reshapedsamps[1,0,:] - reshapedsamps[0,0,:]
    kde_diff_SP = stats.gaussian_kde(diff_pre_post_SP)
    alphasamps = samples['alphahier']
    wherecalc_pos_SP = np.arange(0,np.max(diff_pre_post_SP[diff_pre_post_SP >= 0]),step=stepsize)
    prob_muscimol_poschange_SP = np.sum(kde_diff_SP(wherecalc_pos_SP)*stepsize)
    print("The probability of a positive change in NDT for monkey S was %3.2f%%" % (prob_muscimol_poschange_SP*100))
    tabledata[9,modeltrack] = prob_muscimol_poschange_SP*100

    #Calculate the probability of positive change in symmetric boundary
    reshapedsamps = np.reshape(alphasamps, (ninject, nmonkey, nsamps*nchains))
    diff_pre_post_BT = reshapedsamps[1,1,:] - reshapedsamps[0,1,:]
    kde_diff_BT = stats.gaussian_kde(diff_pre_post_BT)
    wherecalc_pos_BT = np.arange(0, np.max(diff_pre_post_BT[diff_pre_post_BT >= 0]),step=stepsize)
    prob_muscimol_poschange_BT = np.sum(kde_diff_BT(wherecalc_pos_BT)*stepsize)
    print("The probability of a positive change in symmetric boundary for monkey B was %3.2f%%" % (prob_muscimol_poschange_BT*100))
    tabledata[10,modeltrack] = prob_muscimol_poschange_BT*100
    diff_pre_post_SP = reshapedsamps[1,0,:] - reshapedsamps[0,0,:]
    kde_diff_SP = stats.gaussian_kde(diff_pre_post_SP)
    wherecalc_pos_SP = np.arange(0,np.max(diff_pre_post_SP[diff_pre_post_SP >= 0]),step=stepsize)
    prob_muscimol_poschange_SP = np.sum(kde_diff_SP(wherecalc_pos_SP)*stepsize)
    print("The probability of a positive change in symmetric boundary for monkey S was %3.2f%%" % (prob_muscimol_poschange_SP*100))
    tabledata[11,modeltrack] = prob_muscimol_poschange_SP*100

    #Calculate the probability of positive change in lapse proportion
    reshapedsamps = np.reshape(problapsesamps, (ninject, nmonkey, nsamps*nchains))
    diff_pre_post_BT = reshapedsamps[1,1,:] - reshapedsamps[0,1,:]
    kde_diff_BT = stats.gaussian_kde(diff_pre_post_BT)
    wherecalc_pos_BT = np.arange(0, np.max(diff_pre_post_BT[diff_pre_post_BT >= 0]),step=stepsize)
    prob_muscimol_poschange_BT = np.sum(kde_diff_BT(wherecalc_pos_BT)*stepsize)
    print("The probability of a positive change in lapse proportion for monkey B was %3.2f%%" % (prob_muscimol_poschange_BT*100))
    tabledata[12,modeltrack] = prob_muscimol_poschange_BT*100
    diff_pre_post_SP = reshapedsamps[1,0,:] - reshapedsamps[0,0,:]
    kde_diff_SP = stats.gaussian_kde(diff_pre_post_SP)
    wherecalc_pos_SP = np.arange(0,np.max(diff_pre_post_SP[diff_pre_post_SP >= 0]),step=stepsize)
    prob_muscimol_poschange_SP = np.sum(kde_diff_SP(wherecalc_pos_SP)*stepsize)
    print("The probability of a positive change in lapse proportion for monkey S was %3.2f%%" % (prob_muscimol_poschange_SP*100))
    tabledata[13,modeltrack] = prob_muscimol_poschange_SP*100

    modeltrack += 1


# Convert data to table
tabledf = pd.DataFrame(tabledata, columns=modelnames, index=paramnames)
tabledf.to_excel(saveloc,float_format="%.2f")