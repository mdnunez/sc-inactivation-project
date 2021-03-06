# mky3_insample.py - Evaluates prediction of models fit to RT & choice data RT data during the SCInject experiment
#
# Copyright (C) 2021 Michael D. Nunez, <mdnunez1@uci.edu>
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
# 01/13/20      Michael Nunez                 Original code, reference: pdm3b_crossvaleval.m
# 02/27/20      Michael Nunez                 Add comparison plots of psychophysical curves
# 03/02/20      Michael Nunez     Correct equations for p(rightward saccade | rightward orientation) and p(rightward saccade | leftward orientation)
# 03/03/20      Michael Nunez            Fix calculation of percentage choices to IF in 0 coherence condition
# 03/09/20      Michael Nunez               Add comparison plots of  RT distributions
# 03/12/20      Michael Nunez                Adaptive home directory
# 03/20/20      Michael Nunez              Change color scheme, plot model results as PMF curves
# 04/13/20      Michael Nunez             Corrections for model 'PrePostRecHalf' and data 'fixedexpSCInjectGPRTdata'
# 04/28/20      Michael Nunez          Create function pooled_insample() for direct comparison to ChartR pooled data
# 05/20/20      Michael Nunez           Include model 'DayHalfNDTFree' & separate full prediction results by monkey
# 06/05/20      Michael Nunez          Include model mky3_DayHalfSPBound.py
# 11/30/20      Michael Nunez          Include model mky3_PrePostRecAccum.py
# 03/08/21      Michael Nunez          Include model mky3_PrePostRecLinearDrift.py


# Online referencs:
# https://stackoverflow.com/questions/7986567/matplotlib-how-to-set-the-current-figure

# Imports
import numpy as np
import scipy.io as sio
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats.stats import pearsonr
import os
import sys
from scipy.optimize import curve_fit

def rsquared_pred(trueval,predval):
    """
    RSQUARED_PRED  Calculates R^2_prediction for data and statistics derived from data
    """
    divisor = np.sum(np.isfinite(trueval)) -1
    # Mean squared error of prediction
    MSEP = np.nansum(np.power(trueval - predval,2)) / divisor
    # Variance estimate of the true values
    vartrue = np.nansum(np.power(trueval - np.nanmean(trueval),2)) / divisor
    # R-squared definition
    rsquared = 1 - (MSEP / vartrue)
    return rsquared


def fpsychometric(x, gamma, lam, beta, alpha):
        return (gamma + (1 - gamma - lam)) / (1.0 + np.exp(-beta*(x-alpha))) #From Shen and Richards, 2012, Acoustical Society of America


def plot_psychometric(filename):
    """
    Plot psychometric functions for true and predicted data from muscimol experiment
    """

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
    elif filename.find('PrePostRec2Accum') != -1:
        model = 13
        modelname = 'PrePostRec2Accum'
    elif filename.find('PrePostRecLinearDrift') != -1:
        model = 14
        modelname = 'PrePostRecLinearDrift'

    if (model>=4):

        fontsize = 12
        precolor = 'k'
        precolorlite = '#7e7e7e'
        postcolor = '#6666b3'
        postcolorlite = '#b3b3e6'
        reccolor = '#008080'
        reccolorlite = '#80b3b3'

        dataname = filename[0:filename.find('_')]
        print(f'Loading original data {dataname}...')
        truedata = sio.loadmat(f'../data/{dataname}.mat')
        alltarget = np.asarray(truedata['RF_target'],dtype='int16').squeeze()
        allinjectioncond = np.asarray(truedata['injectioncond'],dtype='int16').squeeze()
        allsimpleinject = np.asarray(truedata['simpleinject'],dtype='int16').squeeze()
        if 'training' in truedata.keys():
            print(f'Training data found in {dataname}')
            target = alltarget[truedata['training'][0,:]==1]
            injectioncond = allinjectioncond[truedata['training'][0,:]==1]
            simpleinject = allsimpleinject[truedata['training'][0,:]==1]
        else:
            target = alltarget
            injectioncond = allinjectioncond
            simpleinject = allsimpleinject


        print(f'Loading {filename} in-sample prediction...')
        fileloc = '../modelfits/'
        inpredloc = fileloc + filename + '_inpred.mat'
        inpred = sio.loadmat(inpredloc)

        nconditions = np.max(np.squeeze(inpred['condition'])) + 1
        nsessions = np.max(np.squeeze(inpred['session'])) + 1

        evidenceright = np.array([-36, -24, -17, -10, -5, -3, 0, 3, 5, 10, 17, 24, 36])

        plt.figure("psychometric")
        plt.figure("rtmeans")

        plt.figure("psychometric")
        plt.plot(np.array([-50, 50]), np.array([50, 50]), color='k', LineStyle='--')
        plt.plot(np.array([0, 0]), np.array([0, 100]), color='k', LineStyle='--')

        psychometric_true = np.ones((nsessions, evidenceright.shape[0]))*np.nan
        psychometric_pred = np.ones((nsessions, evidenceright.shape[0]))*np.nan
        meanRT_correct_true = np.ones((nsessions, evidenceright.shape[0]))*np.nan
        meanRT_incorrect_true = np.ones((nsessions, evidenceright.shape[0]))*np.nan
        meanRT_correct_pred = np.ones((nsessions, evidenceright.shape[0]))*np.nan
        meanRT_incorrect_pred = np.ones((nsessions, evidenceright.shape[0]))*np.nan
        nconditions = np.unique(np.squeeze(inpred['condition'])).shape[0]
        for s in range(np.unique(np.squeeze(inpred['session'])).shape[0]):
            for c in range(np.unique(np.squeeze(inpred['condition'])).shape[0]):
                wherecond = (np.squeeze(inpred['condition']) == c)
                wheresession = (np.squeeze(inpred['session']) == s)
                whereboth = (wherecond) & (wheresession)
                thisinjectioncond = injectioncond[whereboth][0]
                thissimpleinject = simpleinject[whereboth][0]


                chose_right_when_right_true = (inpred['rightwardRT'][0,whereboth] > 0) & (target[whereboth] == 2.)
                chose_right_when_left_true = (inpred['rightwardRT'][0,whereboth] > 0) & (target[whereboth] == 1.)
                chose_left_when_right_true = (inpred['rightwardRT'][0,whereboth] < 0) & (target[whereboth] == 2.)
                chose_left_when_left_true = (inpred['rightwardRT'][0,whereboth] < 0) & (target[whereboth] == 1.)

                if (nconditions == 7):
                    if (c != 0):
                        psychometric_true[s, (6-c)] = ((np.sum(chose_right_when_left_true)) / (np.sum(chose_left_when_left_true) + np.sum(chose_right_when_left_true)))*100
                        psychometric_true[s, (c+6)] = ((np.sum(chose_right_when_right_true)) / (np.sum(chose_left_when_right_true) + np.sum(chose_right_when_right_true)))*100
                        
                        # p(rightward saccade | rightward orientation) = p(Correct) - p(leftward saccade) + .5
                        # p(rightward saccade | leftward orientation) = 1 - 2*p(Correct) + p(rightward saccade | rightward orientation)
                        psychometric_pred[s, (c+6)] = (np.mean(inpred['correctRT_pred'][whereboth,:] > 0,axis=None) - np.mean(inpred['rightwardRT_pred'][whereboth,:] < 0,axis=None) + .5)*100
                        psychometric_pred[s, (6-c)] = (1 - 2*np.mean(inpred['correctRT_pred'][whereboth,:] > 0,axis=None))*100 + psychometric_pred[s, (c+6)]
                        # psychometric_pred[s, (c+6)] = (np.mean(inpred['correctRT'][0,whereboth] > 0) - np.mean(inpred['rightwardRT'][0,whereboth] < 0) + .5)*100
                        # psychometric_pred[s, (6-c)] = (1 - 2*np.mean(inpred['correctRT'][0,whereboth] > 0))*100 + psychometric_pred[s, (c+6)]
                    else:
                        psychometric_true[s, 6] = ((np.sum(chose_right_when_left_true) + np.sum(chose_right_when_right_true)) / (chose_right_when_right_true.shape[0]))*100
                        psychometric_pred[s, 6] = np.mean(inpred['rightwardRT_pred'][whereboth,:] > 0,axis=None)*100
                        # psychometric_pred[s, (c+6)] = np.mean(inpred['rightwardRT'][0,whereboth] > 0)*100
                    if (psychometric_pred[s, (c+6)] > 100):
                            psychometric_pred[s, (c+6)] = 100
                            psychometric_pred[s, (6-c)] = 0
                    elif (psychometric_pred[s, (c+6)] < 0):
                        psychometric_pred[s, (c+6)] = 0
                        psychometric_pred[s, (6-c)] = 100
                elif (nconditions == 13):
                    if (c >= 7):
                        psychometric_true[s, c] = ((np.sum(chose_right_when_left_true) + np.sum(chose_right_when_right_true)) / (chose_right_when_right_true.shape[0]))*100
                        psychometric_pred[s, c] = np.mean(inpred['rightwardRT_pred'][whereboth,:] > 0,axis=None)*100
                        meanRT_correct_true[s, c] = np.mean(np.abs(inpred['rightwardRT'][0,whereboth][inpred['rightwardRT'][0,whereboth] > 0]))
                        meanRT_incorrect_true[s, c] = np.mean(np.abs(inpred['rightwardRT'][0,whereboth][inpred['rightwardRT'][0,whereboth] < 0]))
                        allpredRT = inpred['rightwardRT_pred'][whereboth,:].flatten()
                        meanRT_correct_pred[s, c] = np.mean(np.abs(allpredRT[allpredRT > 0]))
                        meanRT_incorrect_pred[s, c] = np.mean(np.abs(allpredRT[allpredRT < 0]))
                    else:
                        psychometric_true[s, (6-c)] = ((np.sum(chose_right_when_left_true) + np.sum(chose_right_when_right_true)) / (chose_right_when_right_true.shape[0]))*100
                        psychometric_pred[s, (6-c)] = np.mean(inpred['rightwardRT_pred'][whereboth,:] > 0,axis=None)*100
                        meanRT_correct_true[s, (6-c)] = np.mean(np.abs(inpred['rightwardRT'][0,whereboth][inpred['rightwardRT'][0,whereboth] < 0]))
                        meanRT_incorrect_true[s, (6-c)] = np.mean(np.abs(inpred['rightwardRT'][0,whereboth][inpred['rightwardRT'][0,whereboth] > 0]))
                        allpredRT = inpred['rightwardRT_pred'][whereboth,:].flatten()
                        meanRT_correct_pred[s, (6-c)] = np.mean(np.abs(allpredRT[allpredRT < 0]))
                        meanRT_incorrect_pred[s, (6-c)] = np.mean(np.abs(allpredRT[allpredRT > 0]))

            if ( (thissimpleinject == 1) | (thissimpleinject > 3) ):
                thiscolor = precolor
                thislite = precolorlite
            elif (thissimpleinject == 2):
                thiscolor = postcolor
                thislite = postcolorlite
            elif (thissimpleinject == 3):
                thiscolor = reccolor
                thislite = reccolorlite
            # popt_true, pcov = curve_fit(fpsychometric, evidenceright, psychometric_true[s,:]/100)
            plt.figure("psychometric")
            plt.plot(evidenceright,psychometric_true[s,:], 'o',color=thiscolor)
            plt.figure("rtmeans")
            plt.plot(evidenceright,meanRT_correct_true[s,:]*1000, 'o', color=thiscolor)
            plt.plot(evidenceright,meanRT_incorrect_true[s,:]*1000, 's', color=thiscolor)
            plt.plot(evidenceright,meanRT_correct_pred[s,:]*1000, linestyle='-', color=thislite)
            plt.plot(evidenceright,meanRT_incorrect_pred[s,:]*1000, linestyle='--', color=thislite)
            # plt.plot(evidenceright, fpsychometric(evidenceright, *popt_true)*100, color=precolorlite)
            
            # plt.plot(evidenceright,psychometric_pred[s,:],color=thislite)
            plt.figure("psychometric")
            try:
                popt_pred, pcov = curve_fit(fpsychometric, evidenceright, psychometric_pred[s,:]/100)
                plt.plot(evidenceright, fpsychometric(evidenceright, *popt_pred)*100, color=thislite)
            except:
                print('PMF curve could not be fit for session %d!' % (s))
                pass
    
        plt.figure("psychometric")
        plt.xlabel('Coherence (%)', fontsize=fontsize)
        plt.xticks(evidenceright)
        plt.ylim((0, 100))
        plt.xlim((-40,40))
        plt.yticks(np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]))
        plt.ylabel('Choices to IF (%)', fontsize=fontsize)
        plt.tick_params(axis='both', which='major', labelsize=fontsize)
        custom_lines = [Line2D([0], [0], color=precolor, marker='o', markersize=8, linestyle='None'),
                Line2D([0], [0], color=postcolor, marker='o', markersize=8, linestyle='None'),
                Line2D([0], [0], color=reccolor, marker='o', markersize=8, linestyle='None'),
                Line2D([0], [0], color=precolorlite, lw=2)]
        plt.legend(custom_lines, ['Pre-Inject or Saline','Post-Inject','Recovery', 'Fitted HDDM'],  loc=2, fontsize=fontsize)
        plt.savefig((f'../../figures/{dataname}/{modelname}_Insamp_PMFcurves.png'), dpi=300, format='png',bbox_inches="tight")
        plt.figure("rtmeans")
        plt.xlabel('Coherence (%)', fontsize=fontsize)
        plt.xticks(evidenceright)
        plt.ylim((300,2000))
        axes = plt.gca()
        ymin, ymax = axes.get_ylim()
        plt.plot(np.array([0, 0]), np.array([ymin, ymax]), color='k', LineStyle='--')
        plt.ylabel('Mean Response Time (ms)', fontsize=fontsize)
        plt.tick_params(axis='both', which='major', labelsize=fontsize)
        custom_lines = [Line2D([0], [0], color=precolor, marker='o', markersize=8, linestyle='None'),
                Line2D([0], [0], color=postcolor, marker='o', markersize=8, linestyle='None'),
                Line2D([0], [0], color=reccolor, marker='o', markersize=8, linestyle='None'),
                Line2D([0], [0], color=precolor, marker='s', markersize=8, linestyle='None'),
                Line2D([0], [0], color=precolorlite, lw=2, linestyle='-'),
                Line2D([0], [0], color=precolorlite, lw=2, linestyle='--')]
        plt.legend(custom_lines, ['Pre-Inject or Saline','Post-Inject','Recovery', 'Incorrect RT', 'Fitted HDDM Correct', 'Fitted HDDM Incorrect'],  loc=2, fontsize=fontsize)
        plt.savefig((f'../../figures/{dataname}/{modelname}_Insamp_RTcurves.png'), dpi=300, format='png',bbox_inches="tight")



def eval_insample(filename):
    # Load model results
    print(f'Loading {filename} in-sample prediction...')
    fileloc = '../modelfits/'
    inpredloc = fileloc + filename + '_inpred.mat'
    inpred = sio.loadmat(inpredloc)
    nconditions = np.max(np.squeeze(inpred['condition'])) + 1
    nsessions = np.max(np.squeeze(inpred['session'])) + 1
    accuracy_true = np.empty((nsessions, nconditions))
    rightward_true = np.empty((nsessions, nconditions))
    leftward_true = np.empty((nsessions, nconditions))
    accuracy_pred = np.empty((nsessions, nconditions))
    rightward_pred = np.empty((nsessions, nconditions))
    leftward_pred = np.empty((nsessions, nconditions))
    rtmean_true = np.empty((nsessions, nconditions))
    rtmean_pred = np.empty((nsessions, nconditions))
    rtmed_true = np.empty((nsessions, nconditions))
    rtmed_pred = np.empty((nsessions, nconditions))
    rt25_true = np.empty((nsessions, nconditions))
    rt25_pred = np.empty((nsessions, nconditions))
    rt75_true = np.empty((nsessions, nconditions))
    rt75_pred = np.empty((nsessions, nconditions))
    for s in range(np.unique(np.squeeze(inpred['session'])).shape[0]):
        for c in range(np.unique(np.squeeze(inpred['condition'])).shape[0]):
            wherecond = (np.squeeze(inpred['condition']) == c)
            wheresession = (np.squeeze(inpred['session']) == s)
            whereboth = (wherecond) & (wheresession)
            rightward_true[s,c] = np.mean(inpred['rightwardRT'][0,whereboth] > 0)
            rightward_pred[s,c] = np.mean(inpred['rightwardRT_pred'][whereboth,:] > 0,axis=None)
            leftward_true[s,c] = np.mean(inpred['rightwardRT'][0,whereboth] < 0)
            leftward_pred[s,c] = np.mean(inpred['rightwardRT_pred'][whereboth,:] < 0,axis=None)
            accuracy_true[s,c] = np.mean(inpred['correctRT'][0,whereboth] > 0)
            accuracy_pred[s,c] = np.mean(inpred['correctRT_pred'][whereboth,:] > 0,axis=None)
            rtmean_true[s,c] = np.mean(np.abs(inpred['correctRT'][0,whereboth]))
            rtmean_pred[s,c] = np.mean(np.abs(inpred['correctRT_pred'][whereboth,:]),axis=None)
            rtmed_true[s,c] = np.quantile(np.abs(inpred['correctRT'][0,whereboth]),0.5)
            rtmed_pred[s,c] = np.quantile(np.abs(inpred['correctRT_pred'][whereboth,:]),0.5,axis=None)
            rt25_true[s,c] = np.quantile(np.abs(inpred['correctRT'][0,whereboth]),0.25)
            rt25_pred[s,c] = np.quantile(np.abs(inpred['correctRT_pred'][whereboth,:]),0.25,axis=None)
            rt75_true[s,c] = np.quantile(np.abs(inpred['correctRT'][0,whereboth]),0.75)
            rt75_pred[s,c] = np.quantile(np.abs(inpred['correctRT_pred'][whereboth,:]),0.75,axis=None)
    print('The R^2_pred of rightward choice percentage is %.3f' % (rsquared_pred(rightward_true.flatten(),rightward_pred.flatten())))
    print('The R^2_pred of leftward choice percentage is %.3f' % (rsquared_pred(leftward_true.flatten(),leftward_pred.flatten())))
    print('The R^2_pred of accuracy is %.3f' % (rsquared_pred(accuracy_true.flatten(),accuracy_pred.flatten())))
    print('The R^2_pred of reaction time means is %.3f' % (rsquared_pred(rtmean_true.flatten(),rtmean_pred.flatten())))
    print('The R^2_pred of reaction time 25th percentiles is %.3f' % (rsquared_pred(rt25_true.flatten(),rt25_pred.flatten())))
    print('The R^2_pred of reaction time medians is %.3f' % (rsquared_pred(rtmed_true.flatten(),rtmed_pred.flatten())))
    print('The R^2_pred of reaction time 75th percentiles is %.3f' % (rsquared_pred(rt75_true.flatten(),rt75_pred.flatten())))
    return rightward_pred, rightward_true


def pooled_insample(filename):
    # Load model results
    dataname = filename[0:filename.find('_')]
    print(f'Loading original data {dataname}...')
    truedata = sio.loadmat(f'../data/{dataname}.mat')
    allsimpleinject = np.asarray(truedata['simpleinject'],dtype='int16').squeeze()
    allmonkeylarge = np.asarray(truedata['monkey'],dtype='int16').squeeze()
    if 'training' in truedata.keys():
        print(f'Training data found in {dataname}')
        simpleinject = allsimpleinject[truedata['training'][0,:]==1]
        monkeylarge = allmonkeylarge[truedata['training'][0,:]==1]
    else:
        simpleinject = allsimpleinject
        monkeylarge = allmonkeylarge


    print(f'Loading {filename} in-sample prediction...')
    fileloc = '../modelfits/'
    inpredloc = fileloc + filename + '_inpred.mat'
    inpred = sio.loadmat(inpredloc)

    nconditions = np.max(np.squeeze(inpred['condition'])) + 1
    nsessions = np.max(np.squeeze(inpred['session'])) + 1

    nmonkeys = np.max(monkeylarge)
    accuracy_true_pooled = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    rightward_true_pooled = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    leftward_true_pooled = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    accuracy_pred_pooled = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    rightward_pred_pooled = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    leftward_pred_pooled = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    rtmean_true_pooled = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    rtmean_pred_pooled = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    rtmed_true_pooled = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    rtmed_pred_pooled = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    rt25_true_pooled = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    rt25_pred_pooled = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    rt75_true_pooled = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    rt75_pred_pooled = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    accuracy_true_presubset = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    rightward_true_presubset = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    leftward_true_presubset = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    accuracy_pred_presubset = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    rightward_pred_presubset = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    leftward_pred_presubset = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    rtmean_true_presubset = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    rtmean_pred_presubset = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    rtmed_true_presubset = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    rtmed_pred_presubset = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    rt25_true_presubset = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    rt25_pred_presubset = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    rt75_true_presubset = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    rt75_pred_presubset = np.ones((nmonkeys, nsessions, nconditions))*np.nan
    for m in range(nmonkeys):
        if m == 0:
            monkeylabel = 'Monkey S'
        elif m == 1:
            monkeylabel = 'Monkey B'
        for s in range(np.unique(np.squeeze(inpred['session'])).shape[0]):
            for c in range(np.unique(np.squeeze(inpred['condition'])).shape[0]):
                wherecond = (np.squeeze(inpred['condition']) == c)
                wheresession = (np.squeeze(inpred['session']) == s)
                wheremonkey = (monkeylarge == (m+1))
                whereinjection = (simpleinject == 1) # Use only the Pre- data
                whereall = (wherecond) & (wheremonkey) & (whereinjection) & (wheresession)
                if np.sum(whereall) > 0:
                    rightward_true_presubset[m,s,c] = np.mean(inpred['rightwardRT'][0,whereall] > 0)
                    rightward_pred_presubset[m,s,c] = np.mean(inpred['rightwardRT_pred'][whereall,:] > 0,axis=None)
                    leftward_true_presubset[m,s,c] = np.mean(inpred['rightwardRT'][0,whereall] < 0)
                    leftward_pred_presubset[m,s,c] = np.mean(inpred['rightwardRT_pred'][whereall,:] < 0,axis=None)
                    accuracy_true_presubset[m,s,c] = np.mean(inpred['correctRT'][0,whereall] > 0)
                    accuracy_pred_presubset[m,s,c] = np.mean(inpred['correctRT_pred'][whereall,:] > 0,axis=None)
                    rtmean_true_presubset[m,s,c] = np.mean(np.abs(inpred['correctRT'][0,whereall]))
                    rtmean_pred_presubset[m,s,c] = np.mean(np.abs(inpred['correctRT_pred'][whereall,:]),axis=None)
                    rtmed_true_presubset[m,s,c] = np.quantile(np.abs(inpred['correctRT'][0,whereall]),0.5)
                    rtmed_pred_presubset[m,s,c] = np.quantile(np.abs(inpred['correctRT_pred'][whereall,:]),0.5,axis=None)
                    rt25_true_presubset[m,s,c] = np.quantile(np.abs(inpred['correctRT'][0,whereall]),0.25)
                    rt25_pred_presubset[m,s,c] = np.quantile(np.abs(inpred['correctRT_pred'][whereall,:]),0.25,axis=None)
                    rt75_true_presubset[m,s,c] = np.quantile(np.abs(inpred['correctRT'][0,whereall]),0.75)
                    rt75_pred_presubset[m,s,c] = np.quantile(np.abs(inpred['correctRT_pred'][whereall,:]),0.75,axis=None)
        theseconds = np.array([0,1,2,3,4,7,8,9,10])
        print('The R^2_pred of rightward choice percentage for %s in Pre-subset is %.3f' % (monkeylabel, rsquared_pred(rightward_true_presubset[m,:,theseconds].flatten(),rightward_pred_presubset[m,:,theseconds].flatten())))
        print('The R^2_pred of leftward choice percentage for %s in Pre-subset is %.3f' % (monkeylabel, rsquared_pred(leftward_true_presubset[m,:,theseconds].flatten(),leftward_pred_presubset[m,:,theseconds].flatten())))
        print('The R^2_pred of accuracy for %s in Pre-subset is %.3f' % (monkeylabel, rsquared_pred(accuracy_true_presubset[m,:,theseconds].flatten(),accuracy_pred_presubset[m,:,theseconds].flatten())))
        print('The R^2_pred of reaction time means for %s in Pre-subset is %.3f' % (monkeylabel, rsquared_pred(rtmean_true_presubset[m,:,theseconds].flatten(),rtmean_pred_presubset[m,:,theseconds].flatten())))
        print('The R^2_pred of reaction time 25th percentiles for %s in Pre-subset is %.3f' % (monkeylabel, rsquared_pred(rt25_true_presubset[m,:,theseconds].flatten(),rt25_pred_presubset[m,:,theseconds].flatten())))
        print('The R^2_pred of reaction time medians for %s in Pre-subset is %.3f' % (monkeylabel, rsquared_pred(rtmed_true_presubset[m,:,theseconds].flatten(),rtmed_pred_presubset[m,:,theseconds].flatten())))
        print('The R^2_pred of reaction time 75th percentiles for %s in Pre-subset is %.3f' % (monkeylabel, rsquared_pred(rt75_true_presubset[m,:,theseconds].flatten(),rt75_pred_presubset[m,:,theseconds].flatten())))
        for s in range(np.unique(np.squeeze(inpred['session'])).shape[0]):
            for c in range(np.unique(np.squeeze(inpred['condition'])).shape[0]):
                wherecond = (np.squeeze(inpred['condition']) == c)
                wheresession = (np.squeeze(inpred['session']) == s)
                wheremonkey = (monkeylarge == (m+1))
                whereall = (wherecond) & (wheremonkey) & (wheresession)
                if np.sum(whereall) > 0:
                    rightward_true_pooled[m,s,c] = np.mean(inpred['rightwardRT'][0,whereall] > 0)
                    rightward_pred_pooled[m,s,c] = np.mean(inpred['rightwardRT_pred'][whereall,:] > 0,axis=None)
                    leftward_true_pooled[m,s,c] = np.mean(inpred['rightwardRT'][0,whereall] < 0)
                    leftward_pred_pooled[m,s,c] = np.mean(inpred['rightwardRT_pred'][whereall,:] < 0,axis=None)
                    accuracy_true_pooled[m,s,c] = np.mean(inpred['correctRT'][0,whereall] > 0)
                    accuracy_pred_pooled[m,s,c] = np.mean(inpred['correctRT_pred'][whereall,:] > 0,axis=None)
                    rtmean_true_pooled[m,s,c] = np.mean(np.abs(inpred['correctRT'][0,whereall]))
                    rtmean_pred_pooled[m,s,c] = np.mean(np.abs(inpred['correctRT_pred'][whereall,:]),axis=None)
                    rtmed_true_pooled[m,s,c] = np.quantile(np.abs(inpred['correctRT'][0,whereall]),0.5)
                    rtmed_pred_pooled[m,s,c] = np.quantile(np.abs(inpred['correctRT_pred'][whereall,:]),0.5,axis=None)
                    rt25_true_pooled[m,s,c] = np.quantile(np.abs(inpred['correctRT'][0,whereall]),0.25)
                    rt25_pred_pooled[m,s,c] = np.quantile(np.abs(inpred['correctRT_pred'][whereall,:]),0.25,axis=None)
                    rt75_true_pooled[m,s,c] = np.quantile(np.abs(inpred['correctRT'][0,whereall]),0.75)
                    rt75_pred_pooled[m,s,c] = np.quantile(np.abs(inpred['correctRT_pred'][whereall,:]),0.75,axis=None)
        print('The R^2_pred of rightward choice percentage for %s is %.3f' % (monkeylabel, rsquared_pred(rightward_true_pooled[m,:,:].flatten(),rightward_pred_pooled[m,:,:].flatten())))
        print('The R^2_pred of leftward choice percentage for %s is %.3f' % (monkeylabel, rsquared_pred(leftward_true_pooled[m,:,:].flatten(),leftward_pred_pooled[m,:,:].flatten())))
        print('The R^2_pred of accuracy for %s is %.3f' % (monkeylabel, rsquared_pred(accuracy_true_pooled[m,:,:].flatten(),accuracy_pred_pooled[m,:,:].flatten())))
        print('The R^2_pred of reaction time means for %s is %.3f' % (monkeylabel, rsquared_pred(rtmean_true_pooled[m,:,:].flatten(),rtmean_pred_pooled[m,:,:].flatten())))
        print('The R^2_pred of reaction time 25th percentiles for %s is %.3f' % (monkeylabel, rsquared_pred(rt25_true_pooled[m,:,:].flatten(),rt25_pred_pooled[m,:,:].flatten())))
        print('The R^2_pred of reaction time medians for %s is %.3f' % (monkeylabel, rsquared_pred(rtmed_true_pooled[m,:,:].flatten(),rtmed_pred_pooled[m,:,:].flatten())))
        print('The R^2_pred of reaction time 75th percentiles for %s is %.3f' % (monkeylabel, rsquared_pred(rt75_true_pooled[m,:,:].flatten(),rt75_pred_pooled[m,:,:].flatten())))
    return inpred, monkeylarge, simpleinject, rightward_pred_pooled, rightward_true_pooled
    


if __name__ == '__main__':
    #When mky3_insample.py is run as a script, do this
    rightward_pred, rightward_true = eval_insample(filename=sys.argv[1])
    inpred, monkeylarge, simpleinject, rightward_pred_pooled, rightward_true_pooled = pooled_insample(filename=sys.argv[1])
    # plot_psychometric(filename=sys.argv[1])