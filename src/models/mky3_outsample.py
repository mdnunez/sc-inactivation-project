# mky3_outsample.py - Evaluates prediction of models fit to RT & choice data RT data during the SCInject experiment with data not used for model training
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
# 04/13/20     Michael Nunez                        Converted from mky3_insample.py
# 05/20/20      Michael Nunez           Include model 'DayHalfNDTFree' & separate full prediction results by monkey
# 06/05/20      Michael Nunez          Include model mky3_DayHalfSPBound.py
# 11/30/20      Michael Nunez          Include model mky3_PrePostRecAccum.py

# Online referencs:
# https://stackoverflow.com/questions/7986567/matplotlib-how-to-set-the-current-figure
# https://stackoverflow.com/questions/179369/how-do-i-abort-the-execution-of-a-python-script

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
        modename = 'DayHalfNDTFree'
    elif filename.find('DayHalfSPBound') != -1:
        model = 11
        modelname = 'DayHalfSPBound'
    elif filename.find('PrePostRec2Accum') != -1:
        model = 13
        modelname = 'PrePostRec2Accum'

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
        allrightwardRT = ((truedata['RF_choice']*2 - 3)*truedata['responsetime']).squeeze()
        allcorrectRT = ((truedata['correct']*2 - 1)*truedata['responsetime']).squeeze()
        allcondition = np.asarray(truedata['condition'],dtype='int16').squeeze() - 1
        allsession = np.asarray(truedata['sessions'],dtype='int16').squeeze() - 1
        if 'training' in truedata.keys():
            print(f'Test data found in {dataname}')
            target = alltarget[truedata['training'][0,:]==0]
            injectioncond = allinjectioncond[truedata['training'][0,:]==0]
            simpleinject = allsimpleinject[truedata['training'][0,:]==0]
            rightwardRT = allrightwardRT[truedata['training'][0,:]==0]
            correctRT = allcorrectRT[truedata['training'][0,:]==0]
            condition = allcondition[truedata['training'][0,:]==0]
            session = allsession[truedata['training'][0,:]==0]
        else:
            sys.exit(f'No test data found in {dataname}! Exiting!')




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
        for s in range(nsessions):
            for c in range(nconditions):
                wherecond = (condition == c)
                wheresession = (session == s)
                whereboth = (wherecond) & (wheresession)
                predcond = (np.squeeze(inpred['condition']) == c)
                predsession = (np.squeeze(inpred['session']) == s)
                predboth = (predcond) & (predsession)
                thisinjectioncond = injectioncond[whereboth][0]
                thissimpleinject = simpleinject[whereboth][0]


                chose_right_when_right_true = (rightwardRT[whereboth] > 0) & (target[whereboth] == 2.)
                chose_right_when_left_true = (rightwardRT[whereboth] > 0) & (target[whereboth] == 1.)
                chose_left_when_right_true = (rightwardRT[whereboth] < 0) & (target[whereboth] == 2.)
                chose_left_when_left_true = (rightwardRT[whereboth] < 0) & (target[whereboth] == 1.)

                if (nconditions == 7):
                    if (c != 0):
                        psychometric_true[s, (6-c)] = ((np.sum(chose_right_when_left_true)) / (np.sum(chose_left_when_left_true) + np.sum(chose_right_when_left_true)))*100
                        psychometric_true[s, (c+6)] = ((np.sum(chose_right_when_right_true)) / (np.sum(chose_left_when_right_true) + np.sum(chose_right_when_right_true)))*100
                        
                        # p(rightward saccade | rightward orientation) = p(Correct) - p(leftward saccade) + .5
                        # p(rightward saccade | leftward orientation) = 1 - 2*p(Correct) + p(rightward saccade | rightward orientation)
                        psychometric_pred[s, (c+6)] = (np.mean(inpred['correctRT_pred'][predboth,:] > 0,axis=None) - np.mean(inpred['rightwardRT_pred'][predboth,:] < 0,axis=None) + .5)*100
                        psychometric_pred[s, (6-c)] = (1 - 2*np.mean(inpred['correctRT_pred'][predboth,:] > 0,axis=None))*100 + psychometric_pred[s, (c+6)]
                        # psychometric_pred[s, (c+6)] = (np.mean(inpred['correctRT'][0,whereboth] > 0) - np.mean(inpred['rightwardRT'][0,whereboth] < 0) + .5)*100
                        # psychometric_pred[s, (6-c)] = (1 - 2*np.mean(inpred['correctRT'][0,whereboth] > 0))*100 + psychometric_pred[s, (c+6)]
                    else:
                        psychometric_true[s, 6] = ((np.sum(chose_right_when_left_true) + np.sum(chose_right_when_right_true)) / (chose_right_when_right_true.shape[0]))*100
                        psychometric_pred[s, 6] = np.mean(inpred['rightwardRT_pred'][predboth,:] > 0,axis=None)*100
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
                        psychometric_pred[s, c] = np.mean(inpred['rightwardRT_pred'][predboth,:] > 0,axis=None)*100
                        meanRT_correct_true[s, c] = np.mean(np.abs(rightwardRT[whereboth][rightwardRT[whereboth] > 0]))
                        meanRT_incorrect_true[s, c] = np.mean(np.abs(rightwardRT[whereboth][rightwardRT[whereboth] < 0]))
                        allpredRT = inpred['rightwardRT_pred'][predboth,:].flatten()
                        meanRT_correct_pred[s, c] = np.mean(np.abs(allpredRT[allpredRT > 0]))
                        meanRT_incorrect_pred[s, c] = np.mean(np.abs(allpredRT[allpredRT < 0]))
                    else:
                        psychometric_true[s, (6-c)] = ((np.sum(chose_right_when_left_true) + np.sum(chose_right_when_right_true)) / (chose_right_when_right_true.shape[0]))*100
                        psychometric_pred[s, (6-c)] = np.mean(inpred['rightwardRT_pred'][predboth,:] > 0,axis=None)*100
                        meanRT_correct_true[s, (6-c)] = np.mean(np.abs(rightwardRT[whereboth][rightwardRT[whereboth] < 0]))
                        meanRT_incorrect_true[s, (6-c)] = np.mean(np.abs(rightwardRT[whereboth][rightwardRT[whereboth] > 0]))
                        allpredRT = inpred['rightwardRT_pred'][predboth,:].flatten()
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
        plt.savefig((f'../../figures/{dataname}/{modelname}_Outsamp_PMFcurves.png'), dpi=300, format='png',bbox_inches="tight")
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
        plt.savefig((f'../../figures/{dataname}/{modelname}_Outsamp_RTcurves.png'), dpi=300, format='png',bbox_inches="tight")



def eval_outsample(filename):
    # Load model results
    print(f'Loading {filename} in-sample prediction...')
    fileloc = '../modelfits/'
    inpredloc = fileloc + filename + '_inpred.mat'
    inpred = sio.loadmat(inpredloc)
    

    dataname = filename[0:filename.find('_')]
    print(f'Loading original data {dataname}...')
    truedata = sio.loadmat(f'../data/{dataname}.mat')
    alltarget = np.asarray(truedata['RF_target'],dtype='int16').squeeze()
    allinjectioncond = np.asarray(truedata['injectioncond'],dtype='int16').squeeze()
    allsimpleinject = np.asarray(truedata['simpleinject'],dtype='int16').squeeze()
    allrightwardRT = ((truedata['RF_choice']*2 - 3)*truedata['responsetime']).squeeze()
    allcorrectRT = ((truedata['correct']*2 - 1)*truedata['responsetime']).squeeze()
    allcondition = np.asarray(truedata['condition'],dtype='int16').squeeze()
    allsession = np.asarray(truedata['sessions'],dtype='int16').squeeze()
    if 'training' in truedata.keys():
        print(f'Test data found in {dataname}')
        target = alltarget[truedata['training'][0,:]==0]
        injectioncond = allinjectioncond[truedata['training'][0,:]==0]
        simpleinject = allsimpleinject[truedata['training'][0,:]==0]
        rightwardRT = allrightwardRT[truedata['training'][0,:]==0]
        correctRT = allcorrectRT[truedata['training'][0,:]==0]
        condition = allcondition[truedata['training'][0,:]==0]
        session = allsession[truedata['training'][0,:]==0]
    else:
        sys.exit(f'No test data found in {dataname}! Exiting!')

    nconditions = np.max(condition)
    nsessions = np.max(session)

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
    for s in range(nsessions):
        for c in range(nconditions):
            wherecond = (condition == (c+1))
            wheresession = (session == (s+1))
            whereboth = (wherecond) & (wheresession)
            predcond = (np.squeeze(inpred['condition']) == c)
            predsession = (np.squeeze(inpred['session']) == s)
            predboth = (predcond) & (predsession)
            rightward_true[s,c] = np.mean(rightwardRT[whereboth] > 0)
            rightward_pred[s,c] = np.mean(inpred['rightwardRT_pred'][predboth,:] > 0,axis=None)
            leftward_true[s,c] = np.mean(rightwardRT[whereboth] < 0)
            leftward_pred[s,c] = np.mean(inpred['rightwardRT_pred'][predboth,:] < 0,axis=None)
            accuracy_true[s,c] = np.mean(correctRT[whereboth] > 0)
            accuracy_pred[s,c] = np.mean(inpred['correctRT_pred'][predboth,:] > 0,axis=None)
            rtmean_true[s,c] = np.mean(np.abs(correctRT[whereboth]))
            rtmean_pred[s,c] = np.mean(np.abs(inpred['correctRT_pred'][predboth,:]),axis=None)
            rtmed_true[s,c] = np.quantile(np.abs(correctRT[whereboth]),0.5)
            rtmed_pred[s,c] = np.quantile(np.abs(inpred['correctRT_pred'][predboth,:]),0.5,axis=None)
            rt25_true[s,c] = np.quantile(np.abs(correctRT[whereboth]),0.25)
            rt25_pred[s,c] = np.quantile(np.abs(inpred['correctRT_pred'][predboth,:]),0.25,axis=None)
            rt75_true[s,c] = np.quantile(np.abs(correctRT[whereboth]),0.75)
            rt75_pred[s,c] = np.quantile(np.abs(inpred['correctRT_pred'][predboth,:]),0.75,axis=None)
    print('The R^2_pred of rightward choice percentage is %.3f' % (rsquared_pred(rightward_true.flatten(),rightward_pred.flatten())))
    print('The R^2_pred of leftward choice percentage is %.3f' % (rsquared_pred(leftward_true.flatten(),leftward_pred.flatten())))
    if (nconditions != 13):
            print('The R^2_pred of accuracy is %.3f' % (rsquared_pred(accuracy_true.flatten(),accuracy_pred.flatten())))
    print('The R^2_pred of reaction time means is %.3f' % (rsquared_pred(rtmean_true.flatten(),rtmean_pred.flatten())))
    print('The R^2_pred of reaction time medians is %.3f' % (rsquared_pred(rtmed_true.flatten(),rtmed_pred.flatten())))
    print('The R^2_pred of reaction time 25th percentiles is %.3f' % (rsquared_pred(rt25_true.flatten(),rt25_pred.flatten())))
    print('The R^2_pred of reaction time 75th percentiles is %.3f' % (rsquared_pred(rt75_true.flatten(),rt75_pred.flatten())))
    return rightward_pred, rightward_true, accuracy_pred, accuracy_true


def pooled_outsample(filename):
    # Load model results
    print(f'Loading {filename} in-sample prediction...')
    fileloc = '../modelfits/'
    inpredloc = fileloc + filename + '_inpred.mat'
    inpred = sio.loadmat(inpredloc)
    

    dataname = filename[0:filename.find('_')]
    print(f'Loading original data {dataname}...')
    truedata = sio.loadmat(f'../data/{dataname}.mat')
    allsimpleinject = np.asarray(truedata['simpleinject'],dtype='int16').squeeze()
    allmonkeylarge = np.asarray(truedata['monkey'],dtype='int16').squeeze()
    allrightwardRT = ((truedata['RF_choice']*2 - 3)*truedata['responsetime']).squeeze()
    allcorrectRT = ((truedata['correct']*2 - 1)*truedata['responsetime']).squeeze()
    allcondition = np.asarray(truedata['condition'],dtype='int16').squeeze()
    allsession = np.asarray(truedata['sessions'],dtype='int16').squeeze()
    if 'training' in truedata.keys():
        print(f'Test data found in {dataname}')
        rightwardRT = allrightwardRT[truedata['training'][0,:]==0]
        correctRT = allcorrectRT[truedata['training'][0,:]==0]
        condition = allcondition[truedata['training'][0,:]==0]
        session = allsession[truedata['training'][0,:]==0]
        monkeylarge_out = allmonkeylarge[truedata['training'][0,:]==0]
        monkeylarge_in = allmonkeylarge[truedata['training'][0,:]==1]
        simpleinject_out = allsimpleinject[truedata['training'][0,:]==0]
        simpleinject_in = allsimpleinject[truedata['training'][0,:]==1]
    else:
        sys.exit(f'No test data found in {dataname}! Exiting!')

    nconditions = np.max(condition)
    nsessions = np.max(session)

    nmonkeys = np.max(monkeylarge_in)
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
        for s in range(nsessions):
            for c in range(np.unique(np.squeeze(inpred['condition'])).shape[0]):
                predcond = (np.squeeze(inpred['condition']) == c)
                predsession = (np.squeeze(inpred['session']) == s)
                predmonkey = (monkeylarge_in == (m+1))
                predinjection = (simpleinject_in == 1) # Use only the Pre- data
                predall = (predcond) & (predmonkey) & (predinjection) & (predsession)
                wherecond = (condition == (c+1))
                wheresession = (session == (s+1))
                wheremonkey = (monkeylarge_out == (m+1))
                whereinjection = (simpleinject_out == 1) # Use only the Pre- data
                whereall = (wherecond) & (wheremonkey) & (whereinjection) & (wheresession)
                if np.sum(whereall) > 0:
                    rightward_true_presubset[m,s,c] = np.mean(rightwardRT[whereall] > 0)
                    rightward_pred_presubset[m,s,c] = np.mean(inpred['rightwardRT_pred'][predall,:] > 0,axis=None)
                    leftward_true_presubset[m,s,c] = np.mean(rightwardRT[whereall] < 0)
                    leftward_pred_presubset[m,s,c] = np.mean(inpred['rightwardRT_pred'][predall,:] < 0,axis=None)
                    accuracy_true_presubset[m,s,c] = np.mean(correctRT[whereall] > 0)
                    accuracy_pred_presubset[m,s,c] = np.mean(inpred['correctRT_pred'][predall,:] > 0,axis=None)
                    rtmean_true_presubset[m,s,c] = np.mean(np.abs(correctRT[whereall]))
                    rtmean_pred_presubset[m,s,c] = np.mean(np.abs(inpred['correctRT_pred'][predall,:]),axis=None)
                    rtmed_true_presubset[m,s,c] = np.quantile(np.abs(correctRT[whereall]),0.5)
                    rtmed_pred_presubset[m,s,c] = np.quantile(np.abs(inpred['correctRT_pred'][predall,:]),0.5,axis=None)
                    rt25_true_presubset[m,s,c] = np.quantile(np.abs(correctRT[whereall]),0.25)
                    rt25_pred_presubset[m,s,c] = np.quantile(np.abs(inpred['correctRT_pred'][predall,:]),0.25,axis=None)
                    rt75_true_presubset[m,s,c] = np.quantile(np.abs(correctRT[whereall]),0.75)
                    rt75_pred_presubset[m,s,c] = np.quantile(np.abs(inpred['correctRT_pred'][predall,:]),0.75,axis=None)
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
                predcond = (np.squeeze(inpred['condition']) == c)
                predsession = (np.squeeze(inpred['session']) == s)
                predmonkey = (monkeylarge_in == (m+1))
                predall = (predcond) & (predmonkey) & (predinjection)
                wherecond = (condition == (c+1))
                wheresession = (session == (s+1))
                wheremonkey = (monkeylarge_out == (m+1))
                whereall = (wherecond) & (wheremonkey) & (whereinjection)
                if np.sum(whereall) > 0:
                    rightward_true_pooled[m,s,c] = np.mean(rightwardRT[whereall] > 0)
                    rightward_pred_pooled[m,s,c] = np.mean(inpred['rightwardRT_pred'][predall,:] > 0,axis=None)
                    leftward_true_pooled[m,s,c] = np.mean(rightwardRT[whereall] < 0)
                    leftward_pred_pooled[m,s,c] = np.mean(inpred['rightwardRT_pred'][predall,:] < 0,axis=None)
                    accuracy_true_pooled[m,s,c] = np.mean(correctRT[whereall] > 0)
                    accuracy_pred_pooled[m,s,c] = np.mean(inpred['correctRT_pred'][predall,:] > 0,axis=None)
                    rtmean_true_pooled[m,s,c] = np.mean(np.abs(correctRT[whereall]))
                    rtmean_pred_pooled[m,s,c] = np.mean(np.abs(inpred['correctRT_pred'][predall,:]),axis=None)
                    rtmed_true_pooled[m,s,c] = np.quantile(np.abs(correctRT[whereall]),0.5)
                    rtmed_pred_pooled[m,s,c] = np.quantile(np.abs(inpred['correctRT_pred'][predall,:]),0.5,axis=None)
                    rt25_true_pooled[m,s,c] = np.quantile(np.abs(correctRT[whereall]),0.25)
                    rt25_pred_pooled[m,s,c] = np.quantile(np.abs(inpred['correctRT_pred'][predall,:]),0.25,axis=None)
                    rt75_true_pooled[m,s,c] = np.quantile(np.abs(correctRT[whereall]),0.75)
                    rt75_pred_pooled[m,s,c] = np.quantile(np.abs(inpred['correctRT_pred'][predall,:]),0.75,axis=None)
        print('The R^2_pred of rightward choice percentage for %s is %.3f' % (monkeylabel, rsquared_pred(rightward_true_pooled[m,:,:].flatten(),rightward_pred_pooled[m,:,:].flatten())))
        print('The R^2_pred of leftward choice percentage for %s is %.3f' % (monkeylabel, rsquared_pred(leftward_true_pooled[m,:,:].flatten(),leftward_pred_pooled[m,:,:].flatten())))
        print('The R^2_pred of accuracy for %s is %.3f' % (monkeylabel, rsquared_pred(accuracy_true_pooled[m,:,:].flatten(),accuracy_pred_pooled[m,:,:].flatten())))
        print('The R^2_pred of reaction time means for %s is %.3f' % (monkeylabel, rsquared_pred(rtmean_true_pooled[m,:,:].flatten(),rtmean_pred_pooled[m,:,:].flatten())))
        print('The R^2_pred of reaction time 25th percentiles for %s is %.3f' % (monkeylabel, rsquared_pred(rt25_true_pooled[m,:,:].flatten(),rt25_pred_pooled[m,:,:].flatten())))
        print('The R^2_pred of reaction time medians for %s is %.3f' % (monkeylabel, rsquared_pred(rtmed_true_pooled[m,:,:].flatten(),rtmed_pred_pooled[m,:,:].flatten())))
        print('The R^2_pred of reaction time 75th percentiles for %s is %.3f' % (monkeylabel, rsquared_pred(rt75_true_pooled[m,:,:].flatten(),rt75_pred_pooled[m,:,:].flatten())))
    return inpred, rightward_pred_presubset, rightward_true_presubset, rightward_pred_pooled, rightward_true_pooled
    

    


if __name__ == '__main__':
    #When mky3_outsample.py is run as a script, do this
    rightward_pred, rightward_true, accuracy_pred, accuracy_true = eval_outsample(filename=sys.argv[1])
    inpred, rightward_pred_presubset, rightward_true_presubset, rightward_pred_pooled, rightward_true_pooled = pooled_outsample(filename=sys.argv[1])
    # plot_psychometric(filename=sys.argv[1])