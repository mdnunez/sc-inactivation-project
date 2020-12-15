# mky3_modelresults.py - Evaluates results of models fit to RT & choice data RT data from Glass Pattern stimuli (e.g. Glass_Pattern_Master.d)
#
# Copyright (C) 2020 Michael D. Nunez, <mdnunez1@uci.edu>
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
# 01/08/20      Michael Nunez                       Converted from mky2_modelresults.py
# 01/14/20      Michael Nunez                        Ability to handle different models
# 01/30/20      Michael Nunez               Plot hierarchical parameters from model PrePostRecFull
# 02/04/20      Michael Nunez                     Add figures for paper
# 02/07/20      Michael Nunez              Change default plots for PrePostRecFull and PrePostRecHalf models
# 02/11/20      Michael Nunez                Return samples and diags with save figures, split into multiple definitions
# 02/19/20      Michael Nunez             Plots of parameters for individual sessions
# 02/20/20      Michael Nunez                Includes results from more models
# 03/02/20      Michael Nunez                Includes results from 'DayFullDriftFree'
# 03/03/20      Michael Nunez                Includes results from 'DayFullSPFree'
# 03/05/20      Michael Nunez                Fixes for models 'DayFullDriftFree' and 'DayFullSPFree'
# 03/12/20      Michael Nunez                Adaptive home directory
# 03/23/20      Michael Nunez               Plot single session parameters on top
# 04/17/20      Michael Nunez      Calculate Rhatnew with new method given by Gelman et al., 2013 (split chains)
# 04/27/20      Michael Nunez             Print out probability of change of deltasamps from pre to post for each monkey
# 04/28/20      Michael Nunez           Separate single-session estimates for each monkey
# 05/18/20      Michael Nunez              Add model 'DayHalfNDTFree'
# 06/02/20      Michael Nunez                More details in subplots of parameters
# 06/03/20      Michael Nunez               Changes to subplots
# 06/05/20      Michael Nunez          Include model mky3_DayHalfSPBound.py
# 06/16/20      Michael Nunez         Include positive and negative changes in all parameter Pre- to Post-
# 07/02/20      Michael Nunez                  Change color scheme for Extended Data Fig. 5
# 07/23/20      Michael Nunez    Change color scheme, split into more subplots, change subplot limits
# 08/06/20      Michael Nunez               Include results from 'PrePostRecDC'
# 08/10/20      Michael Nunez               Include results from 'PrePostRec2Accum'
# 08/12/20      Michael Nunez      Only plot samples of gaininIF from non-fixed values, report BFs
# 08/13/20      Michael Nunez          Invert BF for the first gaininIF parameter for precision
# 11/03/20      Michael Nunez        Adding models prior sensitivity models: 'PrePostRecHalfNarrow', 'PrePostRecHalfShifted', 'PrePostRecHalfWide'
# 11/24/20      Michael Nunez           Invert BF results for start point

#References:
#https://stackoverflow.com/questions/419163/what-does-if-name-main-do
#https://stackoverflow.com/questions/1186789/what-is-the-best-way-to-call-a-script-from-another-script
#https://matplotlib.org/2.0.2/api/_as_gen/matplotlib.axes.Axes.fill_between.html
# http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=4
#https://stackoverflow.com/questions/46183967/how-to-reshape-only-last-dimensions-in-numpy
#https://github.com/stan-dev/pystan/blob/2ce07365061004c1efec8230e8adc353a2092608/pystan/_chains.pyx
#https://stackoverflow.com/questions/4028904/how-to-get-the-home-directory-in-python
#https://docs.scipy.org/doc/numpy/reference/generated/numpy.take.html

# Imports
import numpy as np
import scipy.io as sio
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
    print("Maximum Rhat was %3.2f for variable %s" % (np.max(maxrhats),allkeys[np.argmax(maxrhats)]))
    print("Maximum Rhatnew was %3.2f for variable %s" % (np.max(maxrhatsnew),allkeys[np.argmax(maxrhatsnew)]))
    print("Minimum number of effective samples was %d for variable %s" % (np.min(minneff),allkeys[np.argmin(minneff)]))
    return result


def jellyfish(possamps):  # jellyfish plots
    """Plots posterior distributions of given posterior samples in a jellyfish
    plot. Jellyfish plots are posterior distributions (mirrored over their
    horizontal axes) with 99% and 95% credible intervals (currently plotted
    from the .5% and 99.5% & 2.5% and 97.5% percentiles respectively.
    Also plotted are the median and mean of the posterior distributions"

    Parameters
    ----------
    possamps : ndarray of posterior chains where the last dimension is
    the number of chains, the second to last dimension is the number of samples
    in each chain, all other dimensions describe the shape of the parameter
    """

    # Number of chains
    nchains = possamps.shape[-1]

    # Number of samples per chain
    nsamps = possamps.shape[-2]

    # Number of dimensions
    ndims = possamps.ndim - 2

    # Number of variables to plot
    nvars = np.prod(possamps.shape[0:-2])

    # Index of variables
    varindx = np.arange(nvars).reshape(possamps.shape[0:-2])

    # Reshape data
    alldata = np.reshape(possamps, (nvars, nchains, nsamps))
    alldata = np.reshape(alldata, (nvars, nchains * nsamps))

    # Plot properties
    LineWidths = np.array([2, 5])
    teal = np.array([0, .7, .7])
    blue = np.array([0, 0, 1])
    orange = np.array([1, .3, 0])
    Colors = [teal, blue]

    # Initialize ylabels list
    ylabels = ['']

    for v in range(0, nvars):
        # Create ylabel
        whereis = np.where(varindx == v)
        newlabel = ''
        for l in range(0, ndims):
            newlabel = newlabel + ('_%i' % whereis[l][0])

        ylabels.append(newlabel)

        # Compute posterior density curves
        kde = stats.gaussian_kde(alldata[v, :])
        bounds = stats.scoreatpercentile(alldata[v, :], (.5, 2.5, 97.5, 99.5))
        for b in range(0, 2):
            # Bound by .5th percentile and 99.5th percentile
            x = np.linspace(bounds[b], bounds[-1 - b], 100)
            p = kde(x)

            # Scale distributions down
            maxp = np.max(p)

            # Plot jellyfish
            upper = .25 * p / maxp + v + 1
            lower = -.25 * p / maxp + v + 1
            lines = plt.plot(x, upper, x, lower)
            plt.setp(lines, color=Colors[b], linewidth=LineWidths[b])
            if b == 1:
                # Mark mode
                wheremaxp = np.argmax(p)
                mmode = plt.plot(np.array([1., 1.]) * x[wheremaxp],
                                 np.array([lower[wheremaxp], upper[wheremaxp]]))
                plt.setp(mmode, linewidth=3, color=orange)
                # Mark median
                mmedian = plt.plot(np.median(alldata[v, :]), v + 1, 'ko')
                plt.setp(mmedian, markersize=10, color=[0., 0., 0.])
                # Mark mean
                mmean = plt.plot(np.mean(alldata[v, :]), v + 1, '^')
                plt.setp(mmean, markersize=10, color=teal)

    # Display plot
    plt.setp(plt.gca(), yticklabels=ylabels, yticks=np.arange(0, nvars + 1))


def mm2inch(*tupl):
    mmperinch = 25.4
    if isinstance(tupl[0], tuple):
        return tuple(i/mmperinch for i in tupl[0])
    else:
        return tuple(i/mmperinch for i in tupl)



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


def save_figures(samples, filename):
    # Save figures

    # # =============================================================================
    # # Plots
    # # =============================================================================
    dataname = filename[0:filename.find('_')]
    model = 0
    modelname = 'unknown'
    if filename.find('PrePostRecmodel1') != -1:
        model = 1
        modelname = 'PrePostRecmodel1'
    elif filename.find('PrePostRecmodel2') != -1:
        model = 2
        modelname = 'PrePostRecmodel2'
    elif filename.find('PrePostRecmodel3') != -1:
        model = 3
        modelname = 'PrePostRecmodel3'
    elif filename.find('PrePostRecFull') != -1:
        model = 4
        modelname = 'PrePostRecFull'
    elif filename.find('PrePostRecHalfNarrow') != -1:
        model = 5 #Produce the same results as for 'PrePostRecHalf' but with different file names
        modelname = 'PrePostRecHalfNarrow'
    elif filename.find('PrePostRecHalfShifted') != -1:
        model = 5 #Produce the same results as for 'PrePostRecHalf' but with different file names
        modelname = 'PrePostRecHalfShifted'
    elif filename.find('PrePostRecHalfWide') != -1:
        model = 5 #Produce the same results as for 'PrePostRecHalf' but with different file names
        modelname = 'PrePostRecHalfWide'
    elif filename.find('PrePostRecHalf') != -1:
        model = 5
        modelname = 'PrePostRecHalf'
    elif filename.find('DayHalfDriftFree') != -1:
        model = 6
        modelname = 'DayHalfDriftFree'
    elif filename.find('DayHalfSPFree') != -1:
        model = 7
        modelname = 'DayHalfSPFree'
    elif filename.find('DayFullDriftFree') != -1:
        model = 8
        modelname = 'DayFullDriftFree'
    elif filename.find('DayFullSPFree') != -1:
        model = 9
        modelname = 'DayFullSPFree'
    elif filename.find('DayHalfNDTFree') != -1:
        model = 10
        modelname = 'DayHalfNDTFree'
    elif filename.find('DayHalfSPBound') != -1:
        model = 11
        modelname = 'DayHalfSPBound'
    elif filename.find('PrePostRecDC') != -1:
        model = 12
        modelname = 'PrePostRecDC'
    elif filename.find('PrePostRec2Accum') != -1:
        model = 13
        modelname = 'PrePostRec2Accum'


    if not os.path.exists((f'../../figures/{dataname}')):
       os.makedirs((f'../../figures/{dataname}'));


    fontsize = 7
    panelsize = 8
    markersize = 2
    legendsize = 5
    precolor = 'k'
    postcolor ='#e86800' 
    reccolor = '#008080'
    BTcolor = 'k'
    SPcolor = 'k'

    np.random.seed(seed=2021)

    dataloaded = False

    plt.figure()
    if model > 3:
        jellyfish(samples['terhier']*1000)
    else:
        jellyfish(samples['ter']*1000)
    plt.xlabel('Milliseconds', fontsize=fontsize)
    plt.title('Non-decision time posteriors', fontsize=fontsize)
    plt.savefig((f'../../figures/{dataname}/{modelname}_NDTPosteriors.png'), dpi=300, format='png',bbox_inches="tight")

    plt.figure()
    if (model != 12):
        if (model > 3):
            jellyfish(samples['alphahier'])
        else:
            jellyfish(samples['alpha'])
        plt.xlabel('Evidence units', fontsize=fontsize)
        plt.title('Boundary (Speed-accuracy tradeoff) posteriors', fontsize=fontsize)
        plt.savefig((f'../../figures/{dataname}/{modelname}_AlphaPosteriors.png'), dpi=300, format='png',bbox_inches="tight")
    else:
        jellyfish(samples['varsigmahier'])
        plt.xlabel('Evidence units per second', fontsize=fontsize)
        plt.title('Diffusion coefficient posteriors', fontsize=fontsize)
        plt.savefig((f'../../figures/{dataname}/{modelname}_DCPosteriors.png'), dpi=300, format='png',bbox_inches="tight")

    plt.figure()
    if model > 3:
        jellyfish(samples['betahier'])
    else:
        jellyfish(samples['beta'])
    plt.xlabel('Evidence units', fontsize=fontsize)
    plt.title('Initial evidence bias towards right posteriors', fontsize=fontsize)
    plt.savefig((f'../../figures/{dataname}/{modelname}_BetaPosteriors.png'), dpi=300, format='png',bbox_inches="tight")

    if (model < 5) | ( (model > 7) & (model <10)):
        plt.figure()
        if model > 3:
            jellyfish(samples['deltaacchier'])
        else:
            jellyfish(samples['deltaacc'])
        plt.xlabel('Evidence units per second', fontsize=fontsize)
        plt.title('Accuracy drift rate posteriors', fontsize=fontsize)
        plt.savefig((f'../../figures/{dataname}/{modelname}_DeltaAccPosteriors.png'), dpi=300, format='png',bbox_inches="tight")

    if (model != 13):
        plt.figure()
        if model > 3:
            jellyfish(samples['deltarighthier'])
        else:
            jellyfish(samples['deltaright'])
        plt.xlabel('Evidence units per second', fontsize=fontsize)
        plt.title('Rightward drift rate posteriors', fontsize=fontsize)
        plt.savefig((f'../../figures/{dataname}/{modelname}_DeltaRightPosteriors.png'), dpi=300, format='png',bbox_inches="tight")
    else:
        plt.figure()
        jellyfish(samples['rateparam'])
        plt.xlabel('Seconds to accumulate away IF response when IF probability 0%', fontsize=fontsize)
        plt.title('Accumulation time', fontsize=fontsize)
        plt.savefig((f'../../figures/{dataname}/{modelname}_AccumulationTime.png'), dpi=300, format='png',bbox_inches="tight")
        plt.figure()
        jellyfish(samples['SCinIF']*100)
        plt.xlabel('Pre-muscimol probability of IF response', fontsize=fontsize)
        plt.title('Pre-injection SC encoded probability', fontsize=fontsize)
        plt.savefig((f'../../figures/{dataname}/{modelname}_SCinIF.png'), dpi=300, format='png',bbox_inches="tight")
        plt.figure()
        jellyfish(samples['gaininIF'][1:,:,:])
        plt.xlabel('Gain scalar', fontsize=fontsize)
        plt.title('Effect of each experimental condition on SC', fontsize=fontsize)
        plt.savefig((f'../../figures/{dataname}/{modelname}_InjectionGain.png'), dpi=300, format='png',bbox_inches="tight")
        # Calculate the probability of change of parameter gaininIF
        stepsize = .001
        ninject = samples['gaininIF'][1:,:,:].shape[0] + 1
        nchains = samples['gaininIF'][1:,:,:].shape[-1]
        nsamps = samples['gaininIF'][1:,:,:].shape[-2]
        reshapedsamps = np.reshape(samples['gaininIF'][1:,:,:], (ninject-1, nsamps*nchains))
        median_gaininIF = np.quantile(reshapedsamps,0.5,axis=1)
        lower_gaininIF = np.quantile(reshapedsamps,0.025,axis=1)
        upper_gaininIF = np.quantile(reshapedsamps,0.975,axis=1)
        # Calculate Bayes Factors using the Savage-Dickey density ratio
        denom = stats.norm.pdf(1, loc=1, scale=3) #Matches the gain effect of an injection prior in mky3_PrePostRec2Accum.py
        kde_post = stats.gaussian_kde(reshapedsamps[0,:])
        num_post = kde_post(1)
        print(np.format_float_scientific(num_post, unique=False, precision=15))
        bf_post = num_post / denom
        print("The Bayes Factor of a gaintoIF of NOT 1 in Post-data was %3.4f" % (1/bf_post))
        print("The Bayes Factor of a gaintoIF of 1 in Post-data was %3.4f" % (bf_post))
        print(np.format_float_scientific(bf_post, unique=False, precision=15))
        print("The posterior median of gaintoIF in Post-data was %3.4f" % (median_gaininIF[0]))
        kde_rec = stats.gaussian_kde(reshapedsamps[1,:])
        num_rec = kde_rec(1)
        bf_rec = num_rec / denom
        print("The Bayes Factor of a gaintoIF of NOT 1 in Post-data was %3.4f" % (1/bf_rec))
        print("The posterior median of gaintoIF in Post-data was %3.4f" % (median_gaininIF[1]))
        kde_pre_saline = stats.gaussian_kde(reshapedsamps[2,:])
        num_pre_saline = kde_pre_saline(1)
        bf_pre_saline = num_pre_saline / denom
        print("The Bayes Factor of a gaintoIF of NOT 1 in Pre-saline-data was %3.4f" % (1/bf_pre_saline))
        print("The posterior median of gaintoIF in Pre-saline-data was %3.4f" % (median_gaininIF[2]))
        kde_post_saline = stats.gaussian_kde(reshapedsamps[3,:])
        num_post_saline = kde_post_saline(1)
        bf_post_saline = num_post_saline / denom
        print("The Bayes Factor of a gaintoIF of NOT 1 in Post-saline-data was %3.4f" % (1/bf_post_saline))
        print("The posterior median of gaintoIF in Post-saline-data was %3.4f" % (median_gaininIF[3]))
        kde_rec_saline = stats.gaussian_kde(reshapedsamps[4,:])
        num_rec_saline = kde_rec_saline(1)
        bf_rec_saline = num_rec_saline / denom
        print("The Bayes Factor of a gaintoIF of NOT 1 in Rec-saline-data was %3.4f" % (1/bf_rec_saline))
        print("The posterior median of gaintoIF in Rec-saline-data was %3.4f" % (median_gaininIF[4]))

    plt.figure()
    if model > 3:
        jellyfish(samples['problapsehier']*100)
    else:
        jellyfish(samples['problapse']*100)
    plt.xlabel('Probability', fontsize=fontsize)
    plt.title('Probability of lapse trial posteriors', fontsize=fontsize)
    plt.savefig((f'../../figures/{dataname}/{modelname}_LapsePosteriors.png'), dpi=300, format='png',bbox_inches="tight")


    fig, axs = plt.subplots(4, 3,figsize=mm2inch(244,110), sharex='col')
    
    if ((model < 6) & (model != 0)) | (model == 12):
        if model > 3:
            deltarightsamps = samples['deltarighthier']
        else:
            deltarightsamps = samples['deltarightses']
        nchains = deltarightsamps.shape[-1]
        nsamps = deltarightsamps.shape[-2]
        ninject = deltarightsamps.shape[0]
        if len(deltarightsamps.shape) > 3:
            nmonkey = deltarightsamps.shape[1]
            reshapedsamps = np.reshape(deltarightsamps, (ninject, nmonkey, nsamps*nchains))
            median_deltaright = np.quantile(reshapedsamps,0.5,axis=2)
            lower_deltaright = np.quantile(reshapedsamps,0.025,axis=2)
            upper_deltaright = np.quantile(reshapedsamps,0.975,axis=2)

            axs[0,0].plot(np.arange(3), median_deltaright[3:6,0], color=SPcolor, lw=1, linestyle='dashdot')
            axs[0,0].plot(np.arange(3), lower_deltaright[3:6,0], color=SPcolor, lw=.5, linestyle='dashdot')
            axs[0,0].plot(np.arange(3), upper_deltaright[3:6,0], color=SPcolor, lw=.5, linestyle='dashdot')
            axs[2,0].plot(np.arange(3), median_deltaright[0:3,1], color=BTcolor, lw=1)
            # axs[0,0].plot(np.arange(3), lower_deltaright[0:3,1], color=BTcolor, .5)
            # axs[0,0].plot(np.arange(3), upper_deltaright[0:3,1], color=BTcolor, .5)
            axs[2,0].fill_between(np.arange(3), lower_deltaright[0:3,1], upper_deltaright[0:3,1], color=BTcolor, alpha=0.125)
            axs[0,0].plot(np.arange(3), median_deltaright[0:3,0], color=SPcolor, lw=1)
            # axs[0,0].plot(np.arange(3), lower_deltaright[0:3,0], color=SPcolor, .5)
            # axs[0,0].plot(np.arange(3), upper_deltaright[0:3,0], color=SPcolor, .5)
            axs[0,0].fill_between(np.arange(3), lower_deltaright[0:3,0], upper_deltaright[0:3,0], color=SPcolor, alpha=0.125)
            axs[0,0].set_ylabel('Mean drift rate across coherences', fontsize=fontsize)
            axs[2,0].set_ylabel('Mean drift rate across coherences', fontsize=fontsize)
            # Calculate the probability of change of parameter deltaright
            stepsize = .001
            diff_pre_post_BT = reshapedsamps[1,1,:] - reshapedsamps[0,1,:]
            kde_diff_BT = stats.gaussian_kde(diff_pre_post_BT)
            wherecalc_neg_BT = np.arange(np.min(diff_pre_post_BT[diff_pre_post_BT <= 0]),0,step=stepsize)
            prob_muscimol_negchange_BT = np.sum(kde_diff_BT(wherecalc_neg_BT)*stepsize)
            print("The probability of a negative change in deltaright for monkey B was %3.2f%%" % (prob_muscimol_negchange_BT*100))
            wherecalc_pos_BT = np.arange(0, np.max(diff_pre_post_BT[diff_pre_post_BT >= 0]),step=stepsize)
            prob_muscimol_poschange_BT = np.sum(kde_diff_BT(wherecalc_pos_BT)*stepsize)
            print("The probability of a positive change in deltaright for monkey B was %3.2f%%" % (prob_muscimol_poschange_BT*100))
            diff_pre_post_SP = reshapedsamps[1,0,:] - reshapedsamps[0,0,:]
            kde_diff_SP = stats.gaussian_kde(diff_pre_post_SP)
            wherecalc_neg_SP = np.arange(np.min(diff_pre_post_SP[diff_pre_post_SP <= 0]),0,step=stepsize)
            prob_muscimol_negchange_SP = np.sum(kde_diff_SP(wherecalc_neg_SP)*stepsize)
            print("The probability of a negative change in deltaright for monkey S was %3.2f%%" % (prob_muscimol_negchange_SP*100))
            wherecalc_pos_SP = np.arange(0,np.max(diff_pre_post_SP[diff_pre_post_SP >= 0]),step=stepsize)
            prob_muscimol_poschange_SP = np.sum(kde_diff_SP(wherecalc_pos_SP)*stepsize)
            print("The probability of a positive change in deltaright for monkey S was %3.2f%%" % (prob_muscimol_poschange_SP*100))
            # Calculate Bayes Factors using the Savage-Dickey density ratio
            denom = stats.norm.pdf(0, loc=0, scale=2) #Matches the Hierarchical Rightward/Leftward drift rate prior in mky3_PrePostRecHalf.py
            kde_pre_BT = stats.gaussian_kde(reshapedsamps[0,1,:])
            num_pre_BT = kde_pre_BT(0)
            bf_pre_BT = num_pre_BT / denom
            print("The Bayes Factor for a mean drift of NOT 0 for monkey B in Pre-data was %3.4f" % (1/bf_pre_BT))
            print("The posterior median of  deltaright for monkey B in Pre-data was %3.4f" % (median_deltaright[0,1]))
            kde_post_BT = stats.gaussian_kde(reshapedsamps[1,1,:])
            num_post_BT = kde_post_BT(0)
            bf_post_BT = num_post_BT / denom
            print("The Bayes Factor for a mean drift of NOT 0 for monkey B in Post-data was %3.4f" % (1/bf_post_BT))
            print("The posterior median of  deltaright for monkey B in Post-data was %3.4f" % (median_deltaright[1,1]))
            kde_pre_SP = stats.gaussian_kde(reshapedsamps[0,0,:])
            num_pre_SP = kde_pre_SP(0)
            bf_pre_SP = num_pre_SP / denom
            print("The Bayes Factor for a mean drift of NOT 0 for monkey S in Pre-data was %3.4f" % (1/bf_pre_SP))
            print("The posterior median of  deltaright for monkey S in Pre-data was %3.4f" % (median_deltaright[0,0]))
            kde_post_SP = stats.gaussian_kde(reshapedsamps[1,0,:])
            num_post_SP = kde_post_SP(0)
            bf_post_SP = num_post_SP / denom
            print("The Bayes Factor for a mean drift of NOT 0 for monkey S in Post-data was %3.4f" % (1/bf_post_SP))
            print("The posterior median of  deltaright for monkey S in Post-data was %3.4f" % (median_deltaright[1,0]))
        else:
            median_deltaright = np.quantile(np.reshape(deltarightsamps, (ninject, nsamps*nchains)),0.5, axis=1)
            lower_deltaright = np.quantile(np.reshape(deltarightsamps, (ninject, nsamps*nchains)),0.025,axis=1)
            upper_deltaright = np.quantile(np.reshape(deltarightsamps, (ninject, nsamps*nchains)),0.975,axis=1)
            axs[0,0].plot(np.arange(3), median_deltaright, color='g', lw=1)
            axs[0,0].plot(np.arange(3), lower_deltaright, color='g', lw=.5)
            axs[0,0].plot(np.arange(3), upper_deltaright, color='g', lw=.5)
            axs[0,0].fill_between(np.arange(3), lower_deltaright, upper_deltaright, color='g', alpha=0.125)
            axs[0,0].set_ylabel('Rightward drift rate (evidence units / sec)', fontsize=fontsize)

        if not ((model > 3) & (len(deltarightsamps.shape) > 3)):
            axs[0,0].xticks(np.arange(3), ('Pre-Injection', 'Post-Injection', 'Recovery'),fontsize=fontsize)
            axs[0,0].tick_params(axis='both', which='major', labelsize=fontsize)
            # plt.savefig((f'../../figures/{dataname}/{modelname}_DriftBiasChange.png'), dpi=300, format='png',bbox_inches="tight")

        if (model > 3) & (len(deltarightsamps.shape) > 3):
            deltarightsamps = samples['deltarightses']
            nchains = deltarightsamps.shape[-1]
            nsamps = deltarightsamps.shape[-2]
            nses = deltarightsamps.shape[0]

            # jitter = np.random.uniform(-.1,.1,size=nses)
            jitterSP = np.hstack((np.linspace(-.1,.1,num=7),np.zeros((nses-7))))
            jitterBT = np.hstack((np.linspace(-.05,.05,num=2),np.zeros((nses-2))))
            jitterSPsaline = np.hstack((np.linspace(-.1,.1,num=4),np.zeros((nses-4))))
            
            dataloaded = True
            data = sio.loadmat(f'../data/{dataname}.mat')
            injectlarge = np.asarray(data['simpleinject'],dtype='int16').squeeze()
            session = np.asarray(data['sessions'],dtype='int16').squeeze()
            monkeylarge = np.asarray(data['monkey'],dtype='int16').squeeze()
            monkey = np.empty((nses))
            injectfull = np.empty((nses))
            for s in np.unique(session):
                monkey[s-1] = monkeylarge[session == s][0]
                injectfull[s-1] = injectlarge[session == s][0]
            monkey = np.asarray(monkey, dtype='int16') - 1
            injectfull = np.asarray(injectfull, dtype='int16')
            uniqueinject, inject = np.unique(injectfull,return_inverse=True) #Remove injection conditions without any samples
            
            median_deltaright = np.ones((ninject, nmonkey, nses))*np.nan
            tracker = np.zeros((ninject,nmonkey),dtype='int16')
            for s in range(nses):
                median_deltaright[inject[s],monkey[s],tracker[inject[s],monkey[s]]] = np.quantile(np.reshape(deltarightsamps[s,:,:], (nsamps*nchains)),0.5)
                tracker[inject[s],monkey[s]] += 1
            
            axs[0,0].plot(np.ones(nses)*0 + jitterSPsaline, median_deltaright[3,0,:], 'v', markersize=markersize, color=precolor)
            axs[0,0].plot(np.ones(nses)*1 + jitterSPsaline, median_deltaright[4,0,:], 'v', markersize=markersize, color=postcolor)
            axs[0,0].plot(np.ones(nses)*2 + jitterSPsaline, median_deltaright[5,0,:], 'v', markersize=markersize, color=reccolor)
            axs[2,0].plot(np.ones(nses)*0 + jitterBT, median_deltaright[0,1,:], 'o', markersize=markersize, color=precolor)
            axs[2,0].plot(np.ones(nses)*1 + jitterBT, median_deltaright[1,1,:], 'o', markersize=markersize, color=postcolor)
            axs[2,0].plot(np.ones(nses)*2 + jitterBT, median_deltaright[2,1,:], 'o', markersize=markersize, color=reccolor)
            axs[0,0].plot(np.ones(nses)*0 + jitterSP, median_deltaright[0,0,:], '^', markersize=markersize, color=precolor)
            axs[0,0].plot(np.ones(nses)*1 + jitterSP, median_deltaright[1,0,:], '^', markersize=markersize, color=postcolor)
            axs[0,0].plot(np.ones(nses)*2 + jitterSP, median_deltaright[2,0,:], '^', markersize=markersize, color=reccolor)
            axs[0,0].set_ylabel('Mean drift rate across coherences', fontsize=fontsize)
            axs[0,0].set_ylim(-1.5,0.75)
            axs[0,0].set_yticks(np.array([-1.5, -1.0, -0.5, 0.0, 0.5]))
            axs[0,0].set_xticks(np.arange(3))
            axs[0,0].set_xticklabels(('Pre-Injection', 'Post-Injection', 'Recovery'),fontsize=fontsize)
            axs[0,0].tick_params(axis='both', which='major', labelsize=fontsize)
            custom_lines = [Line2D([0], [0], marker='^', markersize=markersize, color=SPcolor, lw=.5), 
                Line2D([0], [0], color=SPcolor, marker='v', markersize=markersize, linestyle='dashdot', lw=.5)]
            axs[0,0].legend(custom_lines, ['Monkey S Muscimol', 'Monkey S Saline'],  loc=3, fontsize=legendsize)

            axs[2,0].set_ylabel('Mean drift rate across coherences', fontsize=fontsize)
            axs[2,0].set_ylim(-1.5,0.75)
            axs[2,0].set_yticks(np.array([-1.5, -1.0, -0.5, 0.0, 0.5]))
            axs[2,0].set_xticks(np.arange(3))
            axs[2,0].set_xticklabels(('Pre-Injection', 'Post-Injection', 'Recovery'),fontsize=fontsize)
            axs[2,0].tick_params(axis='both', which='major', labelsize=fontsize)
            custom_lines = [Line2D([0], [0], marker='o', markersize=markersize, color=BTcolor, lw=.5)]
            axs[2,0].legend(custom_lines, ['Monkey B Muscimol'],  loc=3, fontsize=legendsize)
            # plt.savefig((f'../../figures/{dataname}/{modelname}_DriftBiasChange_Session.png'), dpi=300, format='png',bbox_inches="tight")
        
            if model==4:
                deltaaccsamps = samples['deltaaccses']
                median_deltaacc = np.ones((ninject, nmonkey, nses))*np.nan
                tracker = np.zeros((ninject,nmonkey),dtype='int16')
                for s in range(nses):
                    median_deltaacc[inject[s],monkey[s],tracker[inject[s],monkey[s]]] = np.quantile(np.reshape(deltaaccsamps[s,:,:], (nsamps*nchains)),0.5)
                    tracker[inject[s],monkey[s]] += 1
                
                plt.figure(dpi=300)
                plt.plot(np.arange(3), median_deltaacc[3:6,0,:], color=SPcolor, lw=1, linestyle='dashdot')
                plt.plot(np.arange(3), median_deltaacc[0:3,1,:], color=BTcolor, lw=1)
                plt.plot(np.arange(3), median_deltaacc[0:3,0,:], color=SPcolor, lw=1)
                plt.ylabel('Drift rate towards correct response (evidence units / sec)', fontsize=fontsize)
                plt.xticks(np.arange(3), ('Pre-Injection', 'Post-Injection', 'Recovery'),fontsize=fontsize)
                plt.tick_params(axis='both', which='major', labelsize=fontsize)
                plt.savefig((f'../../figures/{dataname}/{modelname}_DriftCorrectChange_Session.png'), dpi=300, format='png',bbox_inches="tight")

        if model > 3:
            betasamps = samples['betahier']
        else:
            betasamps = samples['beta']
        nchains = betasamps.shape[-1]
        nsamps = betasamps.shape[-2]
        ninject = betasamps.shape[0]
        plt.axes(axs[0,1])
        if len(betasamps.shape) > 3:
            nmonkey = betasamps.shape[1]
            reshapedsamps = np.reshape(betasamps, (ninject, nmonkey, nsamps*nchains))
            median_beta = np.quantile(reshapedsamps,0.5,axis=2)
            lower_beta = np.quantile(reshapedsamps,0.025,axis=2)
            upper_beta = np.quantile(reshapedsamps,0.975,axis=2)
            axs[0,1].plot(np.arange(3), median_beta[3:6,0], color=SPcolor, lw=1, linestyle='dashdot')
            axs[0,1].plot(np.arange(3), lower_beta[3:6,0], color=SPcolor, lw=.5, linestyle='dashdot')
            axs[0,1].plot(np.arange(3), upper_beta[3:6,0], color=SPcolor, lw=.5, linestyle='dashdot')
            axs[2,1].plot(np.arange(3), median_beta[0:3,1], color=BTcolor, lw=1)
            # axs[0,1].plot(np.arange(3), lower_beta[0:3,1], color=BTcolor, .5)
            # axs[0,1].plot(np.arange(3), upper_beta[0:3,1], color=BTcolor, .5)
            axs[2,1].fill_between(np.arange(3), lower_beta[0:3,1], upper_beta[0:3,1], color=BTcolor, alpha=0.125)
            axs[0,1].plot(np.arange(3), median_beta[0:3,0], color=SPcolor, lw=1)
            # axs[0,1].plot(np.arange(3), lower_beta[0:3,0], color=SPcolor, .5)
            # axs[0,1].plot(np.arange(3), upper_beta[0:3,0], color=SPcolor, .5)
            axs[0,1].fill_between(np.arange(3), lower_beta[0:3,0], upper_beta[0:3,0], color=SPcolor, alpha=0.125)
            axs[0,1].set_ylabel('Start point of evidence (proportion)', fontsize=fontsize)
            axs[2,1].set_ylabel('Start point of evidence (proportion)', fontsize=fontsize)
            # Calculate the probability of change of parameter beta
            stepsize = .001
            diff_pre_post_BT = reshapedsamps[1,1,:] - reshapedsamps[0,1,:]
            kde_diff_BT = stats.gaussian_kde(diff_pre_post_BT)
            wherecalc_neg_BT = np.arange(np.min(diff_pre_post_BT[diff_pre_post_BT <= 0]),0,step=stepsize)
            prob_muscimol_negchange_BT = np.sum(kde_diff_BT(wherecalc_neg_BT)*stepsize)
            print("The probability of a negative change in start point for monkey B was %3.2f%%" % (prob_muscimol_negchange_BT*100))
            wherecalc_pos_BT = np.arange(0, np.max(diff_pre_post_BT[diff_pre_post_BT >= 0]),step=stepsize)
            prob_muscimol_poschange_BT = np.sum(kde_diff_BT(wherecalc_pos_BT)*stepsize)
            print("The probability of a positive change in start point for monkey B was %3.2f%%" % (prob_muscimol_poschange_BT*100))
            diff_pre_post_SP = reshapedsamps[1,0,:] - reshapedsamps[0,0,:]
            kde_diff_SP = stats.gaussian_kde(diff_pre_post_SP)
            wherecalc_neg_SP = np.arange(np.min(diff_pre_post_SP[diff_pre_post_SP <= 0]),0,step=stepsize)
            prob_muscimol_negchange_SP = np.sum(kde_diff_SP(wherecalc_neg_SP)*stepsize)
            print("The probability of a negative change in start point for monkey S was %3.2f%%" % (prob_muscimol_negchange_SP*100))
            wherecalc_pos_SP = np.arange(0,np.max(diff_pre_post_SP[diff_pre_post_SP >= 0]),step=stepsize)
            prob_muscimol_poschange_SP = np.sum(kde_diff_SP(wherecalc_pos_SP)*stepsize)
            print("The probability of a positive change in start point for monkey S was %3.2f%%" % (prob_muscimol_poschange_SP*100))
            # Calculate Bayes Factors using the Savage-Dickey density ratio
            denom = stats.norm.pdf(.5, loc=.5, scale=.25) #Matches the Hierarchical Rightward start point bias prior in mky3_PrePostRecHalf.py
            kde_pre_BT = stats.gaussian_kde(reshapedsamps[0,1,:])
            num_pre_BT = kde_pre_BT(.5)
            bf_pre_BT = num_pre_BT / denom
            print("The Bayes Factor for a start point of .5 for monkey B in Pre-data was %3.4f" % (bf_pre_BT))
            print("The posterior median of  start point for monkey B in Pre-data was %3.4f" % (median_beta[0,1]))
            kde_post_BT = stats.gaussian_kde(reshapedsamps[1,1,:])
            num_post_BT = kde_post_BT(.5)
            bf_post_BT = num_post_BT / denom
            print("The Bayes Factor for a start point of .5 for monkey B in Post-data was %3.4f" % (bf_post_BT))
            print("The posterior median of  start point for monkey B in Post-data was %3.4f" % (median_beta[1,1]))
            kde_pre_SP = stats.gaussian_kde(reshapedsamps[0,0,:])
            num_pre_SP = kde_pre_SP(.5)
            bf_pre_SP = num_pre_SP / denom
            print("The Bayes Factor for a start point of .5 for monkey S in Pre-data was %3.4f" % (bf_pre_SP))
            print("The posterior median of  start point for monkey S in Pre-data was %3.4f" % (median_beta[0,0]))
            kde_post_SP = stats.gaussian_kde(reshapedsamps[1,0,:])
            num_post_SP = kde_post_SP(.5)
            bf_post_SP = num_post_SP / denom
            print("The Bayes Factor for a start point of .5 for monkey S in Post-data was %3.4f" % (bf_post_SP))
            print("The posterior median of  start point for monkey S in Post-data was %3.4f" % (median_beta[1,0]))
        else:
            median_beta = np.quantile(np.reshape(betasamps, (ninject, nsamps*nchains)),0.5, axis=1)
            lower_beta = np.quantile(np.reshape(betasamps, (ninject, nsamps*nchains)),0.025,axis=1)
            upper_beta = np.quantile(np.reshape(betasamps, (ninject, nsamps*nchains)),0.975,axis=1)
            axs[0,1].plot(np.arange(3), median_beta, color='g', lw=1)
            axs[0,1].plot(np.arange(3), lower_beta, color='g', lw=.5)
            axs[0,1].plot(np.arange(3), upper_beta, color='g', lw=.5)
            axs[0,1].fill_between(np.arange(3), lower_beta, upper_beta, color='g', alpha=0.125)
            axs[0,1].set_ylabel('Start point of evidence', fontsize=fontsize)

        if not ((model > 3) & (len(betasamps.shape) > 3)):
            axs[0,1].set_xticks(np.arange(3))
            axs[0,1].set_xticklabels(('Pre-Injection', 'Post-Injection', 'Recovery'),fontsize=fontsize)
            axs[0,1].tick_params(axis='both', which='major', labelsize=fontsize)
            # plt.savefig((f'../../figures/{dataname}/{modelname}_InitialBiasChange.png'), dpi=300, format='png',bbox_inches="tight")
        
        if (model > 3) & (len(betasamps.shape) > 3):
            betasamps = samples['beta']
            nchains = betasamps.shape[-1]
            nsamps = betasamps.shape[-2]
            nses = betasamps.shape[0]
            
            if dataloaded == False:
                dataloaded = True
                data = sio.loadmat(f'../data/{dataname}.mat')
                injectlarge = np.asarray(data['simpleinject'],dtype='int16').squeeze()
                session = np.asarray(data['sessions'],dtype='int16').squeeze()
                monkeylarge = np.asarray(data['monkey'],dtype='int16').squeeze()
                monkey = np.empty((nses))
                injectfull = np.empty((nses))
                for s in np.unique(session):
                    monkey[s-1] = monkeylarge[session == s][0]
                    injectfull[s-1] = injectlarge[session == s][0]
                monkey = np.asarray(monkey, dtype='int16') - 1
                injectfull = np.asarray(injectfull, dtype='int16')
                uniqueinject, inject = np.unique(injectfull,return_inverse=True) #Remove injection conditions without any samples
                
            median_beta = np.ones((ninject, nmonkey, nses))*np.nan
            tracker = np.zeros((ninject,nmonkey),dtype='int16')
            for s in range(nses):
                median_beta[inject[s],monkey[s],tracker[inject[s],monkey[s]]] = np.quantile(np.reshape(betasamps[s,:,:], (nsamps*nchains)),0.5)
                tracker[inject[s],monkey[s]] += 1
            
            # axs[0,1].plot(np.arange(3), median_beta[3:6,0,:], color=SPcolor, linestyle=None)
            # axs[0,1].plot(np.arange(3), median_beta[0:3,1,:], color=BTcolor, lw=1)
            # axs[0,1].plot(np.arange(3), median_beta[0:3,0,:], color=SPcolor, lw=1)
            axs[0,1].plot(np.ones(nses)*0 + jitterSPsaline, median_beta[3,0,:], 'v', markersize=markersize, color=precolor)
            axs[0,1].plot(np.ones(nses)*1 + jitterSPsaline, median_beta[4,0,:], 'v', markersize=markersize, color=postcolor)
            axs[0,1].plot(np.ones(nses)*2 + jitterSPsaline, median_beta[5,0,:], 'v', markersize=markersize, color=reccolor)
            axs[2,1].plot(np.ones(nses)*0 + jitterBT, median_beta[0,1,:], 'o', markersize=markersize, color=precolor)
            axs[2,1].plot(np.ones(nses)*1 + jitterBT, median_beta[1,1,:], 'o', markersize=markersize, color=postcolor)
            axs[2,1].plot(np.ones(nses)*2 + jitterBT, median_beta[2,1,:], 'o', markersize=markersize, color=reccolor)
            axs[0,1].plot(np.ones(nses)*0 + jitterSP, median_beta[0,0,:], '^', markersize=markersize, color=precolor)
            axs[0,1].plot(np.ones(nses)*1 + jitterSP, median_beta[1,0,:], '^', markersize=markersize, color=postcolor)
            axs[0,1].plot(np.ones(nses)*2 + jitterSP, median_beta[2,0,:], '^', markersize=markersize, color=reccolor)
            axs[0,1].set_ylabel('Start point of evidence', fontsize=fontsize)
            axs[0,1].set_ylim(0.35,.70)
            axs[0,1].set_xticks(np.arange(3))
            axs[0,1].set_xticklabels(('Pre-Injection', 'Post-Injection', 'Recovery'),fontsize=fontsize)
            axs[0,1].tick_params(axis='both', which='major', labelsize=fontsize)
            axs[2,1].set_ylabel('Start point of evidence', fontsize=fontsize)
            axs[2,1].set_ylim(0.35,0.70)
            axs[2,1].set_xticks(np.arange(3))
            axs[2,1].set_xticklabels(('Pre-Injection', 'Post-Injection', 'Recovery'),fontsize=fontsize)
            axs[2,1].tick_params(axis='both', which='major', labelsize=fontsize)
            # axs[0,1].legend(custom_lines, ['Monkey S', 'Monkey B', 'Saline (Monkey S)'],  loc=4, fontsize=legendsize)
            # plt.savefig((f'../../figures/{dataname}/{modelname}_InitialBiasChange_Session.png'), dpi=300, format='png',bbox_inches="tight")

        if model > 3:
            tersamps = samples['terhier']*1000
        else:
            tersamps = samples['ter']*1000
        nchains = tersamps.shape[-1]
        nsamps = tersamps.shape[-2]
        ninject = tersamps.shape[0]
        if len(tersamps.shape) > 3:
            nmonkey = tersamps.shape[1]
            reshapedsamps = np.reshape(tersamps, (ninject, nmonkey, nsamps*nchains))
            median_ter = np.quantile(reshapedsamps,0.5,axis=2)
            lower_ter = np.quantile(reshapedsamps,0.025,axis=2)
            upper_ter = np.quantile(reshapedsamps,0.975,axis=2)
            axs[0,2].plot(np.arange(3), median_ter[3:6,0], color=SPcolor, lw=1, linestyle='dashdot')
            axs[0,2].plot(np.arange(3), lower_ter[3:6,0], color=SPcolor, lw=.5, linestyle='dashdot')
            axs[0,2].plot(np.arange(3), upper_ter[3:6,0], color=SPcolor, lw=.5, linestyle='dashdot')
            axs[2,2].plot(np.arange(3), median_ter[0:3,1], color=BTcolor, lw=1)
            # axs[0,2].plot(np.arange(3), lower_ter[0:3,1], color=BTcolor, .5)
            # axs[0,2].plot(np.arange(3), upper_ter[0:3,1], color=BTcolor, .5)
            axs[2,2].fill_between(np.arange(3), lower_ter[0:3,1], upper_ter[0:3,1], color=BTcolor, alpha=0.125)
            axs[0,2].plot(np.arange(3), median_ter[0:3,0], color=SPcolor, lw=1)
            # axs[0,2].plot(np.arange(3), lower_ter[0:3,0], color=SPcolor, .5)
            # axs[0,2].plot(np.arange(3), upper_ter[0:3,0], color=SPcolor, .5)
            axs[0,2].fill_between(np.arange(3), lower_ter[0:3,0], upper_ter[0:3,0], color=SPcolor, alpha=0.125)
            # Calculate the probability of change of parameter ter
            stepsize = .001
            diff_pre_post_BT = reshapedsamps[1,1,:] - reshapedsamps[0,1,:]
            kde_diff_BT = stats.gaussian_kde(diff_pre_post_BT)
            wherecalc_neg_BT = np.arange(np.min(diff_pre_post_BT[diff_pre_post_BT <= 0]),0,step=stepsize)
            prob_muscimol_negchange_BT = np.sum(kde_diff_BT(wherecalc_neg_BT)*stepsize)
            print("The probability of a negative change in NDT for monkey B was %3.2f%%" % (prob_muscimol_negchange_BT*100))
            wherecalc_pos_BT = np.arange(0, np.max(diff_pre_post_BT[diff_pre_post_BT >= 0]),step=stepsize)
            prob_muscimol_poschange_BT = np.sum(kde_diff_BT(wherecalc_pos_BT)*stepsize)
            print("The probability of a positive change in NDT for monkey B was %3.2f%%" % (prob_muscimol_poschange_BT*100))
            diff_pre_post_SP = reshapedsamps[1,0,:] - reshapedsamps[0,0,:]
            kde_diff_SP = stats.gaussian_kde(diff_pre_post_SP)
            wherecalc_neg_SP = np.arange(np.min(diff_pre_post_SP[diff_pre_post_SP <= 0]),0,step=stepsize)
            prob_muscimol_negchange_SP = np.sum(kde_diff_SP(wherecalc_neg_SP)*stepsize)
            print("The probability of a negative change in NDT for monkey S was %3.2f%%" % (prob_muscimol_negchange_SP*100))
            wherecalc_pos_SP = np.arange(0,np.max(diff_pre_post_SP[diff_pre_post_SP >= 0]),step=stepsize)
            prob_muscimol_poschange_SP = np.sum(kde_diff_SP(wherecalc_pos_SP)*stepsize)
            print("The probability of a positive change in NDT for monkey S was %3.2f%%" % (prob_muscimol_poschange_SP*100))
            print("The posterior median of NDT for monkey B in Pre-data was %3.4f" % (median_ter[0,1]))
            print("The posterior median of NDT for monkey B in Post-data was %3.4f" % (median_ter[1,1]))
            print("The posterior median of NDT for monkey S in Pre-data was %3.4f" % (median_ter[0,0]))
            print("The posterior median of NDT for monkey S in Post-data was %3.4f" % (median_ter[1,0]))
        else:
            median_ter = np.quantile(np.reshape(tersamps, (ninject, nsamps*nchains)),0.5, axis=1)
            lower_ter = np.quantile(np.reshape(tersamps, (ninject, nsamps*nchains)),0.025,axis=1)
            upper_ter = np.quantile(np.reshape(tersamps, (ninject, nsamps*nchains)),0.975,axis=1)
            axs[0,2].plot(np.arange(3), median_ter, color='g', lw=1)
            axs[0,2].plot(np.arange(3), lower_ter, color='g', lw=.5)
            axs[0,2].plot(np.arange(3), upper_ter, color='g', lw=.5)
            axs[0,2].fill_between(np.arange(3), lower_ter, upper_ter, color='g', alpha=0.125)

        axs[0,2].set_xticks(np.arange(3))
        axs[0,2].set_xticklabels(('Pre-Injection', 'Post-Injection', 'Recovery'),fontsize=fontsize)
        axs[0,2].set_ylabel('Non-decision time (ms)', fontsize=fontsize)
        axs[0,2].tick_params(axis='both', which='major', labelsize=fontsize)
        axs[0,2].set_ylim(325,650)

        axs[2,2].set_xticks(np.arange(3))
        axs[2,2].set_xticklabels(('Pre-Injection', 'Post-Injection', 'Recovery'),fontsize=fontsize)
        axs[2,2].set_ylabel('Non-decision time (ms)', fontsize=fontsize)
        axs[2,2].tick_params(axis='both', which='major', labelsize=fontsize)
        axs[2,2].set_ylim(325,650)
        # plt.savefig((f'../../figures/{dataname}/{modelname}_NDTChange.png'), dpi=300, format='png',bbox_inches="tight")

        if (model > 3) & (len(tersamps.shape) > 3):
            tersamps = samples['ter']*1000
            nchains = tersamps.shape[-1]
            nsamps = tersamps.shape[-2]
            nses = tersamps.shape[0]
            
            if dataloaded == False:
                dataloaded = True
                data = sio.loadmat(f'../data/{dataname}.mat')
                injectlarge = np.asarray(data['simpleinject'],dtype='int16').squeeze()
                session = np.asarray(data['sessions'],dtype='int16').squeeze()
                monkeylarge = np.asarray(data['monkey'],dtype='int16').squeeze()
                monkey = np.empty((nses))
                injectfull = np.empty((nses))
                for s in np.unique(session):
                    monkey[s-1] = monkeylarge[session == s][0]
                    injectfull[s-1] = injectlarge[session == s][0]
                monkey = np.asarray(monkey, dtype='int16') - 1
                injectfull = np.asarray(injectfull, dtype='int16')
                uniqueinject, inject = np.unique(injectfull,return_inverse=True) #Remove injection conditions without any samples
                
            median_ter = np.ones((ninject, nmonkey, nses))*np.nan
            tracker = np.zeros((ninject,nmonkey),dtype='int16')
            for s in range(nses):
                median_ter[inject[s],monkey[s],tracker[inject[s],monkey[s]]] = np.quantile(np.reshape(tersamps[s,:,:], (nsamps*nchains)),0.5)
                tracker[inject[s],monkey[s]] += 1
            
            axs[0,2].plot(np.ones(nses)*0 + jitterSPsaline, median_ter[3,0,:], 'v', markersize=markersize, color=precolor)
            axs[0,2].plot(np.ones(nses)*1 + jitterSPsaline, median_ter[4,0,:], 'v', markersize=markersize, color=postcolor)
            axs[0,2].plot(np.ones(nses)*2 + jitterSPsaline, median_ter[5,0,:], 'v', markersize=markersize, color=reccolor)
            axs[2,2].plot(np.ones(nses)*0 + jitterBT, median_ter[0,1,:], 'o', markersize=markersize, color=precolor)
            axs[2,2].plot(np.ones(nses)*1 + jitterBT, median_ter[1,1,:], 'o', markersize=markersize, color=postcolor)
            axs[2,2].plot(np.ones(nses)*2 + jitterBT, median_ter[2,1,:], 'o', markersize=markersize, color=reccolor)
            axs[0,2].plot(np.ones(nses)*0 + jitterSP, median_ter[0,0,:], '^', markersize=markersize, color=precolor)
            axs[0,2].plot(np.ones(nses)*1 + jitterSP, median_ter[1,0,:], '^', markersize=markersize, color=postcolor)
            axs[0,2].plot(np.ones(nses)*2 + jitterSP, median_ter[2,0,:], '^', markersize=markersize, color=reccolor)
            axs[0,2].set_ylabel('Non-decision time (ms)', fontsize=fontsize)
            axs[0,2].set_xticks(np.arange(3))
            axs[0,2].set_xticklabels(('Pre-Injection', 'Post-Injection', 'Recovery'),fontsize=fontsize)
            axs[0,2].tick_params(axis='both', which='major', labelsize=fontsize)
            axs[2,2].set_ylabel('Non-decision time (ms)', fontsize=fontsize)
            axs[2,2].set_xticks(np.arange(3))
            axs[2,2].set_xticklabels(('Pre-Injection', 'Post-Injection', 'Recovery'),fontsize=fontsize)
            axs[2,2].tick_params(axis='both', which='major', labelsize=fontsize)
            # axs[0,2].legend(custom_lines, ['Monkey S', 'Monkey B', 'Saline (Monkey S)'],  loc=4, fontsize=legendsize)
            # plt.savefig((f'../../figures/{dataname}/{modelname}_InitialBiasChange_Session.png'), dpi=300, format='png',bbox_inches="tight")


        if model > 3:
            problapsesamps = samples['problapsehier']*100
        else:
            problapsesamps = samples['problapse']*100
        nchains = problapsesamps.shape[-1]
        nsamps = problapsesamps.shape[-2]
        ninject = problapsesamps.shape[0]
        if len(problapsesamps.shape) > 3:
            nmonkey = problapsesamps.shape[1]
            reshapedsamps = np.reshape(problapsesamps, (ninject, nmonkey, nsamps*nchains))
            median_problapse = np.quantile(reshapedsamps,0.5,axis=2)
            lower_problapse = np.quantile(reshapedsamps,0.025,axis=2)
            upper_problapse = np.quantile(reshapedsamps,0.975,axis=2)
            axs[1,1].plot(np.arange(3), median_problapse[3:6,0], color=SPcolor, lw=1, linestyle='dashdot')
            axs[1,1].plot(np.arange(3), lower_problapse[3:6,0], color=SPcolor, lw=.5, linestyle='dashdot')
            axs[1,1].plot(np.arange(3), upper_problapse[3:6,0], color=SPcolor, lw=.5, linestyle='dashdot')
            axs[3,1].plot(np.arange(3), median_problapse[0:3,1], color=BTcolor, lw=1)
            # axs[1,1].plot(np.arange(3), lower_problapse[0:3,1], color=BTcolor, .5)
            # axs[1,1].plot(np.arange(3), upper_problapse[0:3,1], color=BTcolor, .5)
            axs[3,1].fill_between(np.arange(3), lower_problapse[0:3,1], upper_problapse[0:3,1], color=BTcolor, alpha=0.125)
            axs[1,1].plot(np.arange(3), median_problapse[0:3,0], color=SPcolor, lw=1)
            # axs[1,1].plot(np.arange(3), lower_problapse[0:3,0], color=SPcolor, .5)
            # axs[1,1].plot(np.arange(3), upper_problapse[0:3,0], color=SPcolor, .5)
            axs[1,1].fill_between(np.arange(3), lower_problapse[0:3,0], upper_problapse[0:3,0], color=SPcolor, alpha=0.125)
            # Calculate the probability of change of parameter problapse
            stepsize = .001
            diff_pre_post_BT = reshapedsamps[1,1,:] - reshapedsamps[0,1,:]
            kde_diff_BT = stats.gaussian_kde(diff_pre_post_BT)
            wherecalc_neg_BT = np.arange(np.min(diff_pre_post_BT[diff_pre_post_BT <= 0]),0,step=stepsize)
            prob_muscimol_negchange_BT = np.sum(kde_diff_BT(wherecalc_neg_BT)*stepsize)
            print("The probability of a negative change in lapse proportion for monkey B was %3.2f%%" % (prob_muscimol_negchange_BT*100))
            wherecalc_pos_BT = np.arange(0, np.max(diff_pre_post_BT[diff_pre_post_BT >= 0]),step=stepsize)
            prob_muscimol_poschange_BT = np.sum(kde_diff_BT(wherecalc_pos_BT)*stepsize)
            print("The probability of a positive change in lapse proportion for monkey B was %3.2f%%" % (prob_muscimol_poschange_BT*100))
            diff_pre_post_SP = reshapedsamps[1,0,:] - reshapedsamps[0,0,:]
            kde_diff_SP = stats.gaussian_kde(diff_pre_post_SP)
            wherecalc_neg_SP = np.arange(np.min(diff_pre_post_SP[diff_pre_post_SP <= 0]),0,step=stepsize)
            prob_muscimol_negchange_SP = np.sum(kde_diff_SP(wherecalc_neg_SP)*stepsize)
            print("The probability of a negative change in lapse proportion for monkey S was %3.2f%%" % (prob_muscimol_negchange_SP*100))
            wherecalc_pos_SP = np.arange(0,np.max(diff_pre_post_SP[diff_pre_post_SP >= 0]),step=stepsize)
            prob_muscimol_poschange_SP = np.sum(kde_diff_SP(wherecalc_pos_SP)*stepsize)
            print("The probability of a positive change in lapse proportion for monkey S was %3.2f%%" % (prob_muscimol_poschange_SP*100))
            print("The posterior median of proportion lapse trials for monkey B in Pre-data was %3.4f" % (median_problapse[0,1]))
            print("The posterior median of proportion lapse trials for monkey B in Post-data was %3.4f" % (median_problapse[1,1]))
            print("The posterior median of proportion lapse trials for monkey S in Pre-data was %3.4f" % (median_problapse[0,0]))
            print("The posterior median of proportion lapse trials for monkey S in Post-data was %3.4f" % (median_problapse[1,0]))
        else:
            median_problapse = np.quantile(np.reshape(problapsesamps, (ninject, nsamps*nchains)),0.5, axis=1)
            lower_problapse = np.quantile(np.reshape(problapsesamps, (ninject, nsamps*nchains)),0.025,axis=1)
            upper_problapse = np.quantile(np.reshape(problapsesamps, (ninject, nsamps*nchains)),0.975,axis=1)
            axs[1,1].plot(np.arange(3), median_problapse, color='g', lw=1)
            axs[1,1].plot(np.arange(3), lower_problapse, color='g', lw=.5)
            axs[1,1].plot(np.arange(3), upper_problapse, color='g', lw=.5)
            axs[1,1].fill_between(np.arange(3), lower_problapse, upper_problapse, color='g', alpha=0.125)

        axs[1,1].set_xticks(np.arange(3))
        axs[1,1].set_xticklabels(('Pre-Injection', 'Post-Injection', 'Recovery'),fontsize=fontsize)
        axs[1,1].set_ylabel('Probability of lapse trial', fontsize=fontsize)
        axs[1,1].tick_params(axis='both', which='major', labelsize=fontsize)
        axs[1,1].set_ylim(0,8)
        axs[3,1].set_xticks(np.arange(3))
        axs[3,1].set_xticklabels(('Pre-Injection', 'Post-Injection', 'Recovery'),fontsize=fontsize)
        axs[3,1].set_ylabel('Probability of lapse trial', fontsize=fontsize)
        axs[3,1].tick_params(axis='both', which='major', labelsize=fontsize)
        axs[3,1].set_ylim(0,8)
        # plt.savefig((f'../../figures/{dataname}/{modelname}_LapseChange.png'), dpi=300, format='png',bbox_inches="tight")

        if (model > 3) & (len(problapsesamps.shape) > 3):
            problapsesamps = samples['problapse']*100
            nchains = problapsesamps.shape[-1]
            nsamps = problapsesamps.shape[-2]
            nses = problapsesamps.shape[0]
            
            if dataloaded == False:
                dataloaded = True
                data = sio.loadmat(f'../data/{dataname}.mat')
                injectlarge = np.asarray(data['simpleinject'],dtype='int16').squeeze()
                session = np.asarray(data['sessions'],dtype='int16').squeeze()
                monkeylarge = np.asarray(data['monkey'],dtype='int16').squeeze()
                monkey = np.empty((nses))
                injectfull = np.empty((nses))
                for s in np.unique(session):
                    monkey[s-1] = monkeylarge[session == s][0]
                    injectfull[s-1] = injectlarge[session == s][0]
                monkey = np.asarray(monkey, dtype='int16') - 1
                injectfull = np.asarray(injectfull, dtype='int16')
                uniqueinject, inject = np.unique(injectfull,return_inverse=True) #Remove injection conditions without any samples
                
            median_prob = np.ones((ninject, nmonkey, nses))*np.nan
            tracker = np.zeros((ninject,nmonkey),dtype='int16')
            for s in range(nses):
                median_prob[inject[s],monkey[s],tracker[inject[s],monkey[s]]] = np.quantile(np.reshape(problapsesamps[s,:,:], (nsamps*nchains)),0.5)
                tracker[inject[s],monkey[s]] += 1
            
            axs[1,1].plot(np.ones(nses)*0 + jitterSPsaline, median_prob[3,0,:], 'v', markersize=markersize, color=precolor)
            axs[1,1].plot(np.ones(nses)*1 + jitterSPsaline, median_prob[4,0,:], 'v', markersize=markersize, color=postcolor)
            axs[1,1].plot(np.ones(nses)*2 + jitterSPsaline, median_prob[5,0,:], 'v', markersize=markersize, color=reccolor)
            axs[3,1].plot(np.ones(nses)*0 + jitterBT, median_prob[0,1,:], 'o', markersize=markersize, color=precolor)
            axs[3,1].plot(np.ones(nses)*1 + jitterBT, median_prob[1,1,:], 'o', markersize=markersize, color=postcolor)
            axs[3,1].plot(np.ones(nses)*2 + jitterBT, median_prob[2,1,:], 'o', markersize=markersize, color=reccolor)
            axs[1,1].plot(np.ones(nses)*0 + jitterSP, median_prob[0,0,:], '^', markersize=markersize, color=precolor)
            axs[1,1].plot(np.ones(nses)*1 + jitterSP, median_prob[1,0,:], '^', markersize=markersize, color=postcolor)
            axs[1,1].plot(np.ones(nses)*2 + jitterSP, median_prob[2,0,:], '^', markersize=markersize, color=reccolor)
            axs[1,1].set_ylabel('Probability of lapse trial', fontsize=fontsize)
            axs[1,1].set_ylim(0,8)
            axs[1,1].set_xticks(np.arange(3))
            axs[1,1].set_xticklabels(('Pre-Injection', 'Post-Injection', 'Recovery'),fontsize=fontsize)
            axs[1,1].tick_params(axis='both', which='major', labelsize=fontsize)
            axs[3,1].set_ylabel('Probability of lapse trial', fontsize=fontsize)
            axs[3,1].set_xticks(np.arange(3))
            axs[3,1].set_xticklabels(('Pre-Injection', 'Post-Injection', 'Recovery'),fontsize=fontsize)
            axs[3,1].tick_params(axis='both', which='major', labelsize=fontsize)
            axs[3,1].set_ylim(0,8)
            # axs[0,2].legend(custom_lines, ['Monkey S', 'Monkey B', 'Saline (Monkey S)'],  loc=4, fontsize=legendsize)
            # plt.savefig((f'../../figures/{dataname}/{modelname}_InitialBiasChange_Session.png'), dpi=300, format='png',bbox_inches="tight")


        if (model != 12):
            if model > 3:
                alphasamps = samples['alphahier']
            else:
                alphasamps = samples['alpha']
            nchains = alphasamps.shape[-1]
            nsamps = alphasamps.shape[-2]
            ninject = alphasamps.shape[0]
            # plt.figure(dpi=300)
            if len(alphasamps.shape) > 3:
                nmonkey = alphasamps.shape[1]
                reshapedsamps = np.reshape(alphasamps, (ninject, nmonkey, nsamps*nchains))
                median_alpha = np.quantile(reshapedsamps,0.5,axis=2)
                lower_alpha = np.quantile(reshapedsamps,0.025,axis=2)
                upper_alpha = np.quantile(reshapedsamps,0.975,axis=2)
                axs[1,0].plot(np.arange(3), median_alpha[3:6,0], color=SPcolor, lw=1, linestyle='dashdot')
                axs[1,0].plot(np.arange(3), lower_alpha[3:6,0], color=SPcolor, lw=.5, linestyle='dashdot')
                axs[1,0].plot(np.arange(3), upper_alpha[3:6,0], color=SPcolor, lw=.5, linestyle='dashdot')
                axs[3,0].plot(np.arange(3), median_alpha[0:3,1], color=BTcolor, lw=1)
                # axs[1,0].plot(np.arange(3), lower_alpha[0:3,1], color=BTcolor, .5)
                # axs[1,0].plot(np.arange(3), upper_alpha[0:3,1], color=BTcolor, .5)
                axs[3,0].fill_between(np.arange(3), lower_alpha[0:3,1], upper_alpha[0:3,1], color=BTcolor, alpha=0.125)
                axs[1,0].plot(np.arange(3), median_alpha[0:3,0], color=SPcolor, lw=1)
                # axs[1,0].plot(np.arange(3), lower_alpha[0:3,0], color=SPcolor, .5)
                # axs[1,0].plot(np.arange(3), upper_alpha[0:3,0], color=SPcolor, .5)
                axs[1,0].fill_between(np.arange(3), lower_alpha[0:3,0], upper_alpha[0:3,0], color=SPcolor, alpha=0.125)
                # Calculate the probability of change of parameter problapse
                stepsize = .001
                diff_pre_post_BT = reshapedsamps[1,1,:] - reshapedsamps[0,1,:]
                kde_diff_BT = stats.gaussian_kde(diff_pre_post_BT)
                wherecalc_neg_BT = np.arange(np.min(diff_pre_post_BT[diff_pre_post_BT <= 0]),0,step=stepsize)
                prob_muscimol_negchange_BT = np.sum(kde_diff_BT(wherecalc_neg_BT)*stepsize)
                print("The probability of a negative change in symmetric boundary for monkey B was %3.2f%%" % (prob_muscimol_negchange_BT*100))
                wherecalc_pos_BT = np.arange(0, np.max(diff_pre_post_BT[diff_pre_post_BT >= 0]),step=stepsize)
                prob_muscimol_poschange_BT = np.sum(kde_diff_BT(wherecalc_pos_BT)*stepsize)
                print("The probability of a positive change in symmetric boundary for monkey B was %3.2f%%" % (prob_muscimol_poschange_BT*100))
                diff_pre_post_SP = reshapedsamps[1,0,:] - reshapedsamps[0,0,:]
                kde_diff_SP = stats.gaussian_kde(diff_pre_post_SP)
                wherecalc_neg_SP = np.arange(np.min(diff_pre_post_SP[diff_pre_post_SP <= 0]),0,step=stepsize)
                prob_muscimol_negchange_SP = np.sum(kde_diff_SP(wherecalc_neg_SP)*stepsize)
                print("The probability of a negative change in symmetric boundary for monkey S was %3.2f%%" % (prob_muscimol_negchange_SP*100))
                wherecalc_pos_SP = np.arange(0,np.max(diff_pre_post_SP[diff_pre_post_SP >= 0]),step=stepsize)
                prob_muscimol_poschange_SP = np.sum(kde_diff_SP(wherecalc_pos_SP)*stepsize)
                print("The probability of a positive change in symmetric boundary for monkey S was %3.2f%%" % (prob_muscimol_poschange_SP*100))
                print("The posterior median of symmetric boundary for monkey B in Pre-data was %3.4f" % (median_alpha[0,1]))
                print("The posterior median of symmetric boundary for monkey B in Post-data was %3.4f" % (median_alpha[1,1]))
                print("The posterior median of symmetric boundary for monkey S in Pre-data was %3.4f" % (median_alpha[0,0]))
                print("The posterior median of symmetric boundary for monkey S in Post-data was %3.4f" % (median_alpha[1,0]))
            else:
                median_alpha = np.quantile(np.reshape(alphasamps, (ninject, nsamps*nchains)),0.5, axis=1)
                lower_alpha = np.quantile(np.reshape(alphasamps, (ninject, nsamps*nchains)),0.025,axis=1)
                upper_alpha = np.quantile(np.reshape(alphasamps, (ninject, nsamps*nchains)),0.975,axis=1)
                axs[1,0].plot(np.arange(3), median_alpha, color='g', lw=1)
                axs[1,0].plot(np.arange(3), lower_alpha, color='g', lw=.5)
                axs[1,0].plot(np.arange(3), upper_alpha, color='g', lw=.5)
                axs[1,0].fill_between(np.arange(3), lower_alpha, upper_alpha, color='g', alpha=0.125)
            
            
            axs[1,0].set_xticks(np.arange(3))
            axs[1,0].set_xticklabels(('Pre-Injection', 'Post-Injection', 'Recovery'),fontsize=fontsize)
            axs[1,0].set_ylabel('Symmetric boundary', fontsize=fontsize)
            axs[1,0].set_ylim(1.1,1.9)
            axs[3,0].set_xticks(np.arange(3))
            axs[3,0].set_xticklabels(('Pre-Injection', 'Post-Injection', 'Recovery'),fontsize=fontsize)
            axs[3,0].set_ylabel('Symmetric boundary', fontsize=fontsize)
            axs[3,0].set_ylim(1.1,1.9)
            # plt.tick_params(axis='both', which='major', labelsize=fontsize)
            # plt.savefig((f'../../figures/{dataname}/{modelname}_BoundaryChange.png'), dpi=300, format='png',bbox_inches="tight")

            if (model > 3) & (len(alphasamps.shape) > 3):
                alphasamps = samples['alpha']
                nchains = alphasamps.shape[-1]
                nsamps = alphasamps.shape[-2]
                nses = alphasamps.shape[0]
            
            if dataloaded == False:
                dataloaded = True
                data = sio.loadmat(f'../data/{dataname}.mat')
                injectlarge = np.asarray(data['simpleinject'],dtype='int16').squeeze()
                session = np.asarray(data['sessions'],dtype='int16').squeeze()
                monkeylarge = np.asarray(data['monkey'],dtype='int16').squeeze()
                monkey = np.empty((nses))
                injectfull = np.empty((nses))
                for s in np.unique(session):
                    monkey[s-1] = monkeylarge[session == s][0]
                    injectfull[s-1] = injectlarge[session == s][0]
                monkey = np.asarray(monkey, dtype='int16') - 1
                injectfull = np.asarray(injectfull, dtype='int16')
                uniqueinject, inject = np.unique(injectfull,return_inverse=True) #Remove injection conditions without any samples
                
            median_alpha = np.ones((ninject, nmonkey, nses))*np.nan
            tracker = np.zeros((ninject,nmonkey),dtype='int16')
            for s in range(nses):
                median_alpha[inject[s],monkey[s],tracker[inject[s],monkey[s]]] = np.quantile(np.reshape(alphasamps[s,:,:], (nsamps*nchains)),0.5)
                tracker[inject[s],monkey[s]] += 1
            
            axs[1,0].plot(np.ones(nses)*0 + jitterSPsaline, median_alpha[3,0,:], 'v', markersize=markersize, color=precolor)
            axs[1,0].plot(np.ones(nses)*1 + jitterSPsaline, median_alpha[4,0,:], 'v', markersize=markersize, color=postcolor)
            axs[1,0].plot(np.ones(nses)*2 + jitterSPsaline, median_alpha[5,0,:], 'v', markersize=markersize, color=reccolor)
            axs[3,0].plot(np.ones(nses)*0 + jitterBT, median_alpha[0,1,:], 'o', markersize=markersize, color=precolor)
            axs[3,0].plot(np.ones(nses)*1 + jitterBT, median_alpha[1,1,:], 'o', markersize=markersize, color=postcolor)
            axs[3,0].plot(np.ones(nses)*2 + jitterBT, median_alpha[2,1,:], 'o', markersize=markersize, color=reccolor)
            axs[1,0].plot(np.ones(nses)*0 + jitterSP, median_alpha[0,0,:], '^', markersize=markersize, color=precolor)
            axs[1,0].plot(np.ones(nses)*1 + jitterSP, median_alpha[1,0,:], '^', markersize=markersize, color=postcolor)
            axs[1,0].plot(np.ones(nses)*2 + jitterSP, median_alpha[2,0,:], '^', markersize=markersize, color=reccolor)
            axs[1,0].set_ylabel('Symmetric boundary', fontsize=fontsize)
            axs[1,0].set_xticks(np.arange(3))
            axs[1,0].set_xticklabels(('Pre-Injection', 'Post-Injection', 'Recovery'),fontsize=fontsize)
            axs[1,0].tick_params(axis='both', which='major', labelsize=fontsize)
            axs[3,0].set_ylabel('Symmetric boundary', fontsize=fontsize)
            axs[3,0].set_xticks(np.arange(3))
            axs[3,0].set_xticklabels(('Pre-Injection', 'Post-Injection', 'Recovery'),fontsize=fontsize)
            axs[3,0].tick_params(axis='both', which='major', labelsize=fontsize)
            plt.subplots_adjust(hspace=0.2, wspace=0.3)
            fig.delaxes(axs[1,2])
            fig.delaxes(axs[3,2])
            fig.set_size_inches(mm2inch(218,147),forward=False)
            plt.savefig((f'../../figures/{dataname}/{modelname}_Results.png'), dpi=300, format='png',bbox_inches='tight')
            plt.savefig((f'../../figures/{dataname}/{modelname}_Results.pdf'), dpi=300, format='pdf',bbox_inches='tight')
            plt.savefig((f'../../figures/{dataname}/{modelname}_Results.svg'), dpi=300, format='svg',bbox_inches='tight')
        else:
            varsigmasamps = samples['varsigmahier']
            nchains = varsigmasamps.shape[-1]
            nsamps = varsigmasamps.shape[-2]
            ninject = varsigmasamps.shape[0]
            nmonkey = varsigmasamps.shape[1]
            reshapedsamps = np.reshape(varsigmasamps, (ninject, nmonkey, nsamps*nchains))
            median_varsigma = np.quantile(reshapedsamps,0.5,axis=2)
            lower_varsigma = np.quantile(reshapedsamps,0.025,axis=2)
            upper_varsigma = np.quantile(reshapedsamps,0.975,axis=2)
            axs[1,0].plot(np.arange(3), median_varsigma[3:6,0], color=SPcolor, lw=1, linestyle='dashdot')
            axs[1,0].plot(np.arange(3), lower_varsigma[3:6,0], color=SPcolor, lw=.5, linestyle='dashdot')
            axs[1,0].plot(np.arange(3), upper_varsigma[3:6,0], color=SPcolor, lw=.5, linestyle='dashdot')
            axs[3,0].plot(np.arange(3), median_varsigma[0:3,1], color=BTcolor, lw=1)
            # axs[1,0].plot(np.arange(3), lower_varsigma[0:3,1], color=BTcolor, .5)
            # axs[1,0].plot(np.arange(3), upper_varsigma[0:3,1], color=BTcolor, .5)
            axs[3,0].fill_between(np.arange(3), lower_varsigma[0:3,1], upper_varsigma[0:3,1], color=BTcolor, alpha=0.125)
            axs[1,0].plot(np.arange(3), median_varsigma[0:3,0], color=SPcolor, lw=1)
            # axs[1,0].plot(np.arange(3), lower_varsigma[0:3,0], color=SPcolor, .5)
            # axs[1,0].plot(np.arange(3), upper_varsigma[0:3,0], color=SPcolor, .5)
            axs[1,0].fill_between(np.arange(3), lower_varsigma[0:3,0], upper_varsigma[0:3,0], color=SPcolor, alpha=0.125)
            # Calculate the probability of change of parameter problapse
            stepsize = .001
            diff_pre_post_BT = reshapedsamps[1,1,:] - reshapedsamps[0,1,:]
            kde_diff_BT = stats.gaussian_kde(diff_pre_post_BT)
            wherecalc_neg_BT = np.arange(np.min(diff_pre_post_BT[diff_pre_post_BT <= 0]),0,step=stepsize)
            prob_muscimol_negchange_BT = np.sum(kde_diff_BT(wherecalc_neg_BT)*stepsize)
            print("The probability of a negative change in accumulation noise for monkey B was %3.2f%%" % (prob_muscimol_negchange_BT*100))
            wherecalc_pos_BT = np.arange(0, np.max(diff_pre_post_BT[diff_pre_post_BT >= 0]),step=stepsize)
            prob_muscimol_poschange_BT = np.sum(kde_diff_BT(wherecalc_pos_BT)*stepsize)
            print("The probability of a positive change in accumulation noise for monkey B was %3.2f%%" % (prob_muscimol_poschange_BT*100))
            diff_pre_post_SP = reshapedsamps[1,0,:] - reshapedsamps[0,0,:]
            kde_diff_SP = stats.gaussian_kde(diff_pre_post_SP)
            wherecalc_neg_SP = np.arange(np.min(diff_pre_post_SP[diff_pre_post_SP <= 0]),0,step=stepsize)
            prob_muscimol_negchange_SP = np.sum(kde_diff_SP(wherecalc_neg_SP)*stepsize)
            print("The probability of a negative change in accumulation noise for monkey S was %3.2f%%" % (prob_muscimol_negchange_SP*100))
            wherecalc_pos_SP = np.arange(0,np.max(diff_pre_post_SP[diff_pre_post_SP >= 0]),step=stepsize)
            prob_muscimol_poschange_SP = np.sum(kde_diff_SP(wherecalc_pos_SP)*stepsize)
            print("The probability of a positive change in accumulation noise for monkey S was %3.2f%%" % (prob_muscimol_poschange_SP*100))
            print("The posterior median of accumulation noise for monkey B in Pre-data was %3.4f" % (median_varsigma[0,1]))
            print("The posterior median of accumulation noise for monkey B in Post-data was %3.4f" % (median_varsigma[1,1]))
            print("The posterior median of accumulation noise for monkey S in Pre-data was %3.4f" % (median_varsigma[0,0]))
            print("The posterior median of accumulation noise for monkey S in Post-data was %3.4f" % (median_varsigma[1,0]))
            
            
            axs[1,0].set_xticks(np.arange(3))
            axs[1,0].set_xticklabels(('Pre-Injection', 'Post-Injection', 'Recovery'),fontsize=fontsize)
            axs[1,0].set_ylabel('Accumulation noise', fontsize=fontsize)
            axs[1,0].set_ylim(0.75,1.25)
            axs[3,0].set_xticks(np.arange(3))
            axs[3,0].set_xticklabels(('Pre-Injection', 'Post-Injection', 'Recovery'),fontsize=fontsize)
            axs[3,0].set_ylabel('Accumulation noise', fontsize=fontsize)
            axs[3,0].set_ylim(0.75,1.25)

            varsigmasamps = samples['varsigma']
            nchains = varsigmasamps.shape[-1]
            nsamps = varsigmasamps.shape[-2]
            nses = varsigmasamps.shape[0]
            
            if dataloaded == False:
                dataloaded = True
                data = sio.loadmat(f'../data/{dataname}.mat')
                injectlarge = np.asarray(data['simpleinject'],dtype='int16').squeeze()
                session = np.asarray(data['sessions'],dtype='int16').squeeze()
                monkeylarge = np.asarray(data['monkey'],dtype='int16').squeeze()
                monkey = np.empty((nses))
                injectfull = np.empty((nses))
                for s in np.unique(session):
                    monkey[s-1] = monkeylarge[session == s][0]
                    injectfull[s-1] = injectlarge[session == s][0]
                monkey = np.asarray(monkey, dtype='int16') - 1
                injectfull = np.asarray(injectfull, dtype='int16')
                uniqueinject, inject = np.unique(injectfull,return_inverse=True) #Remove injection conditions without any samples
                
            median_varsigma = np.ones((ninject, nmonkey, nses))*np.nan
            tracker = np.zeros((ninject,nmonkey),dtype='int16')
            for s in range(nses):
                median_varsigma[inject[s],monkey[s],tracker[inject[s],monkey[s]]] = np.quantile(np.reshape(varsigmasamps[s,:,:], (nsamps*nchains)),0.5)
                tracker[inject[s],monkey[s]] += 1
            
            axs[1,0].plot(np.ones(nses)*0 + jitterSPsaline, median_varsigma[3,0,:], 'v', markersize=markersize, color=precolor)
            axs[1,0].plot(np.ones(nses)*1 + jitterSPsaline, median_varsigma[4,0,:], 'v', markersize=markersize, color=postcolor)
            axs[1,0].plot(np.ones(nses)*2 + jitterSPsaline, median_varsigma[5,0,:], 'v', markersize=markersize, color=reccolor)
            axs[3,0].plot(np.ones(nses)*0 + jitterBT, median_varsigma[0,1,:], 'o', markersize=markersize, color=precolor)
            axs[3,0].plot(np.ones(nses)*1 + jitterBT, median_varsigma[1,1,:], 'o', markersize=markersize, color=postcolor)
            axs[3,0].plot(np.ones(nses)*2 + jitterBT, median_varsigma[2,1,:], 'o', markersize=markersize, color=reccolor)
            axs[1,0].plot(np.ones(nses)*0 + jitterSP, median_varsigma[0,0,:], '^', markersize=markersize, color=precolor)
            axs[1,0].plot(np.ones(nses)*1 + jitterSP, median_varsigma[1,0,:], '^', markersize=markersize, color=postcolor)
            axs[1,0].plot(np.ones(nses)*2 + jitterSP, median_varsigma[2,0,:], '^', markersize=markersize, color=reccolor)
            axs[1,0].set_ylabel('Accumulation noise', fontsize=fontsize)
            axs[1,0].set_xticks(np.arange(3))
            axs[1,0].set_xticklabels(('Pre-Injection', 'Post-Injection', 'Recovery'),fontsize=fontsize)
            axs[1,0].tick_params(axis='both', which='major', labelsize=fontsize)
            axs[3,0].set_ylabel('Accumulation noise', fontsize=fontsize)
            axs[3,0].set_xticks(np.arange(3))
            axs[3,0].set_xticklabels(('Pre-Injection', 'Post-Injection', 'Recovery'),fontsize=fontsize)
            axs[3,0].tick_params(axis='both', which='major', labelsize=fontsize)
            plt.subplots_adjust(hspace=0.2, wspace=0.3)
            fig.delaxes(axs[1,2])
            fig.delaxes(axs[3,2])
            fig.set_size_inches(mm2inch(218,147),forward=False)
            plt.savefig((f'../../figures/{dataname}/{modelname}_Results.png'), dpi=300, format='png',bbox_inches='tight')
            plt.savefig((f'../../figures/{dataname}/{modelname}_Results.pdf'), dpi=300, format='pdf',bbox_inches='tight')
            plt.savefig((f'../../figures/{dataname}/{modelname}_Results.svg'), dpi=300, format='svg',bbox_inches='tight')


        if (model != 5) & (model !=12):
            if model > 3:
                deltaaccsamps = samples['deltaacchier']
            else:
                deltaaccsamps = samples['deltaacc']
            nchains = deltaaccsamps.shape[-1]
            nsamps = deltaaccsamps.shape[-2]
            ninject = deltaaccsamps.shape[0]
            plt.figure(dpi=300)
            if len(deltaaccsamps.shape) > 3:
                nmonkey = deltaaccsamps.shape[1]
                median_deltaacc = np.quantile(np.reshape(deltaaccsamps, (ninject, nmonkey, nsamps*nchains)),0.5,axis=2)
                lower_deltaacc = np.quantile(np.reshape(deltaaccsamps, (ninject, nmonkey, nsamps*nchains)),0.025,axis=2)
                upper_deltaacc = np.quantile(np.reshape(deltaaccsamps, (ninject, nmonkey, nsamps*nchains)),0.975,axis=2)
                plt.plot(np.arange(3), median_deltaacc[3:6,0], color=SPcolor, lw=1, linestyle='dashdot')
                plt.plot(np.arange(3), lower_deltaacc[3:6,0], color=SPcolor, lw=.5, linestyle='dashdot')
                plt.plot(np.arange(3), upper_deltaacc[3:6,0], color=SPcolor, lw=.5, linestyle='dashdot')
                plt.plot(np.arange(3), median_deltaacc[0:3,1], color=BTcolor, lw=1)
                # plt.plot(np.arange(3), lower_deltaacc[0:3,1], color=BTcolor, .5)
                # plt.plot(np.arange(3), upper_deltaacc[0:3,1], color=BTcolor, .5)
                plt.fill_between(np.arange(3), lower_deltaacc[0:3,1], upper_deltaacc[0:3,1], color=BTcolor, alpha=0.125)
                plt.plot(np.arange(3), median_deltaacc[0:3,0], color=SPcolor, lw=1)
                # plt.plot(np.arange(3), lower_deltaacc[0:3,0], color=SPcolor, .5)
                # plt.plot(np.arange(3), upper_deltaacc[0:3,0], color=SPcolor, .5)
                plt.fill_between(np.arange(3), lower_deltaacc[0:3,0], upper_deltaacc[0:3,0], color=SPcolor, alpha=0.125)
            else:
                median_deltaacc = np.quantile(np.reshape(deltaaccsamps, (ninject, nsamps*nchains)),0.5, axis=1)
                lower_deltaacc = np.quantile(np.reshape(deltaaccsamps, (ninject, nsamps*nchains)),0.025,axis=1)
                upper_deltaacc = np.quantile(np.reshape(deltaaccsamps, (ninject, nsamps*nchains)),0.975,axis=1)
                plt.plot(np.arange(3), median_deltaacc, color='g', lw=1)
                plt.plot(np.arange(3), lower_deltaacc, color='g', lw=.5)
                plt.plot(np.arange(3), upper_deltaacc, color='g', lw=.5)
                plt.fill_between(np.arange(3), lower_deltaacc, upper_deltaacc, color='g', alpha=0.125)

            plt.xticks(np.arange(3), ('Pre-Injection', 'Post-Injection', 'Recovery'),fontsize=fontsize)
            plt.ylabel('Drift rate towards correct response (evidence units / sec)', fontsize=fontsize)
            plt.tick_params(axis='both', which='major', labelsize=fontsize)
            plt.savefig((f'../../figures/{dataname}/{modelname}_DriftCorrectChange.png'), dpi=300, format='png',bbox_inches="tight")
    
    if model <= 3:
        plt.figure(figsize=(12,6))
        trials = samples['DDMorLapse'].shape[0]
        for n in range(0,nchains):
            plt.plot(np.squeeze(np.mean(samples['DDMorLapse'][:, :, n], axis=1)), 'o', markersize=markersize, label=('chain %d' % (n+1)))
        plt.title('Decision or lapse trial', fontsize=fontsize)
        plt.xlabel('Trial number', fontsize=fontsize)
        plt.ylabel('Trial type', fontsize=fontsize)
        plt.legend()
        plt.savefig((f'../../figures/{dataname}/{modelname}_DDMorLapse.png'), dpi=300, format='png',bbox_inches="tight")

    if 'deltarightsesdiff' in samples.keys():
        plt.figure()
        jellyfish(samples['deltarightsesdiff'])
        plt.xlabel('Evidence units per second', fontsize=fontsize)
        plt.title('Rightward drift rate differences posteriors', fontsize=fontsize)
        plt.savefig((f'../../figures/{dataname}/{modelname}_DeltaRightDiff.png'), dpi=300, format='png',bbox_inches="tight")

    if 'terdiff' in samples.keys():
        plt.figure()
        jellyfish(samples['terdiff']*1000)
        plt.xlabel('Milliseconds', fontsize=fontsize)
        plt.title('Non-decision time differences posteriors', fontsize=fontsize)
        plt.savefig((f'../../figures/{dataname}/{modelname}_NDTDiffPosteriors.png'), dpi=300, format='png',bbox_inches="tight")

    if 'alphadiff' in samples.keys():
        plt.figure()
        jellyfish(samples['alphadiff'])
        plt.xlabel('Evidence units', fontsize=fontsize)
        plt.title('Boundary differences posteriors', fontsize=fontsize)
        plt.savefig((f'../../figures/{dataname}/{modelname}_AlphaDiffPosteriors.png'), dpi=300, format='png',bbox_inches="tight")

    if 'betadiff' in samples.keys():
        plt.figure()
        jellyfish(samples['betadiff'])
        plt.xlabel('Evidence units', fontsize=fontsize)
        plt.title('Initial evidence differences posteriors', fontsize=fontsize)
        plt.savefig((f'../../figures/{dataname}/{modelname}_BetaDiffPosteriors.png'), dpi=300, format='png',bbox_inches="tight")

    if 'deltaaccsesdiff' in samples.keys():
        plt.figure()
        jellyfish(samples['deltaaccsesdiff'])
        plt.xlabel('Evidence units per second', fontsize=fontsize)
        plt.title('Accuracy drift rate differences posteriors', fontsize=fontsize)
        plt.savefig((f'../../figures/{dataname}/{modelname}_DeltaAccDiffPosteriors.png'), dpi=300, format='png',bbox_inches="tight")

if __name__ == '__main__':
    #When mky3_modelresults.py is run as a script, do this
    (samples, diags, filename) = load_data(filename=sys.argv[1])
    save_figures(samples, filename)