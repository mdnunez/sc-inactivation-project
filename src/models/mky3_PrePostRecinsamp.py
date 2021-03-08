# mky3_PrePostRecinsamp.py - Generates behavioral data from PrePostRecmodel1 then compares to real data
#
# Copyright (C) 2021 Michael D. Nunez, <mdnunez@mednet.ucla.edu>
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
# 01/13/20      Michael Nunez                             Original code
# 01/14/20      Michael Nunez                     Import mky3_insample
# 01/15/20      Michael Nunez             Adaptive for different sizes of problapse (model1, model2, model3)
# 02/27/20      Michael Nunez        Include models PrePostRecHalf, DayHalfDriftFree, DayHalfSPFree
# 02/28/20      Michael Nunez                      Bug fix
# 03/02/20      Michael Nunez                   uselapse flag
# 03/05/20      Michael Nunez          Include models 'DayFullDriftFree' and 'Day FullSPFree'
# 03/12/20      Michael Nunez                Adaptive home directory
# 03/23/20      Michael Nunez           Run mky3_insample.plot_psychometric
# 05/19/20      Michael Nunez           Include model 'DayHalfNDTFree'
# 06/05/20      Michael Nunez      Correct generation posterior predictives for 'DayHalfNDTFree', Include model mky3_DayHalfSPBound.py
# 06/08/20      Michael Nunez      Correct generation posterior predictives for 'DayHalfSPBound'
# 11/30/20      Michael Nunez          Include model mky3_PrePostRecAccum.py
# 03/08/21      Michael Nunez          Include model mky3_PrePostRecLinearDrift.py

import numpy as np
import scipy.io as sio
import scipy as sp
import os
import sys
import warnings
import random
import pandas as pd
import mky3_insample

# Simulate diffusion models
def simuldiff(N=100,Alpha=1,Ter=.4,Nu=1,Zeta=None,rangeTer=0,rangeZeta=0,Eta=.3,Varsigma=1):
    """
    SIMULDIFF  Generates data according to a diffusion model


    Reference:
    Tuerlinckx, F., Maris, E.,
    Ratcliff, R., & De Boeck, P. (2001). A comparison of four methods for
    simulating the diffusion process. Behavior Research Methods,
    Instruments, & Computers, 33, 443-456.

    Parameters
    ----------
    N: a integer denoting the size of the output vector
    (defaults to 100 experimental trials)

    Alpha: the mean boundary separation across trials  in evidence units
    (defaults to 1 evidence unit)

    Ter: the mean non-decision time across trials in seconds
    (defaults to .4 seconds)

    Nu: the mean drift rate across trials in evidence units per second
    (defaults to 1 evidence units per second, restricted to -5 to 5 units)

    Zeta: the initial bias in the evidence process for choice A
    (defaults to 50% of total evidence units given by Alpha)

    rangeTer: Non-decision time across trials is generated from a uniform
    distribution of Ter - rangeTer/2 to  Ter + rangeTer/2 across trials
    (defaults to 0 seconds)

    rangeZeta: Bias across trials is generated from a uniform distribution
    of Zeta - rangeZeta/2 to Zeta + rangeZeta/2 across trials
    (defaults to 0 evidence units)

    Eta: Standard deviation of the drift rate across trials
    (defaults to 3 evidence units per second, restricted to less than 3 evidence units)

    Varsigma: The diffusion coefficient, the standard deviation of the
    evidence accumulation process within one trial. It is recommended that
    this parameter be kept fixed unless you have reason to explore this parameter
    (defaults to 1 evidence unit per second)

    Returns
    -------
    Numpy complex vector with 1) Real component ( np.real(x) ): reaction times (in seconds) multiplied by the response vector
    such that negative reaction times encode response B and positive reaction times
    encode response A  and 2) Imaginary component ( np.imag(x) ): N200 peak-latencies in seconds
    
    
    Converted from simuldiff.m MATLAB script by Joachim Vandekerckhove
    See also http://ppw.kuleuven.be/okp/dmatoolbox.
    """

    if Zeta is None:
        Zeta = .5*Alpha

    if (Nu < -5) or (Nu > 5):
        Nu = np.sign(Nu)*5
        warnings.warn('Nu is not in the range [-5 5], bounding drift rate to %.1f...' % (Nu))

    if (Eta > 3):
        warning.warn('Standard deviation of drift rate is out of bounds, bounding drift rate to 3')
        eta = 3

    if (Eta == 0):
        Eta = 1e-16

    #Initialize output vectors
    result = np.zeros(N)
    T = np.zeros(N)
    XX = np.zeros(N)
    N200 = np.zeros(N)

    #Called sigma in 2001 paper
    D = np.power(Varsigma,2)/2

    #Program specifications
    eps = 2.220446049250313e-16 #precision from 1.0 to next double-precision number
    delta=eps

    for n in range(0,N):
        r1 = np.random.normal()
        mu = Nu + r1*Eta
        zz = Zeta - rangeZeta/2 + rangeZeta*np.random.uniform()
        finish = 0
        totaltime = 0
        startpos = 0
        Aupper = Alpha - zz
        Alower = -zz
        radius = np.min(np.array([np.abs(Aupper), np.abs(Alower)]))
        while (finish==0):
            lambda_ = 0.25*np.power(mu,2)/D + 0.25*D*np.power(np.pi,2)/np.power(radius,2)
            # eq. formula (13) in 2001 paper with D = sigma^2/2 and radius = Alpha/2
            F = D*np.pi/(radius*mu)
            F = np.power(F,2)/(1 + np.power(F,2) )
            # formula p447 in 2001 paper
            prob = np.exp(radius*mu/D)
            prob = prob/(1 + prob)
            dir_ = 2*(np.random.uniform() < prob) - 1
            l = -1
            s2 = 0
            while (s2>l):
                s2=np.random.uniform()
                s1=np.random.uniform()
                tnew=0
                told=0
                uu=0
                while (np.abs(tnew-told)>eps) or (uu==0):
                    told=tnew
                    uu=uu+1
                    tnew = told + (2*uu+1) * np.power(-1,uu) * np.power(s1,(F*np.power(2*uu+1,2)));
                    # infinite sum in formula (16) in BRMIC,2001
                l = 1 + np.power(s1,(-F)) * tnew;
            # rest of formula (16)
            t = np.abs(np.log(s1))/lambda_;
            # is the negative of t* in (14) in BRMIC,2001
            totaltime=totaltime+t
            dir_=startpos+dir_*radius
            NDTime = Ter - rangeTer/2 + rangeTer*np.random.uniform()
            if ( (dir_ + delta) > Aupper):
                T[n]=NDTime+totaltime
                XX[n]=1
                finish=1
            elif ( (dir_-delta) < Alower ):
                T[n]=NDTime+totaltime
                XX[n]=-1
                finish=1
            else:
                startpos=dir_
                radius=np.min(np.abs([Aupper, Alower]-startpos))

    result = T*XX
    return result


def generate_data(filename,seed=11, Niter=100):
    #Load posterior samples
    print(f'Loading {filename}...')
    fileloc = '../modelfits/'
    modelloc = fileloc + filename + '.mat'
    samples = sio.loadmat(modelloc)
    # Number of chains
    nchains = samples['ter'].shape[-1]

    # Number of samples per chain
    nsamps = samples['ter'].shape[-2]

    #Default trial-to-trial variability parameters
    correctEta = 0
    rightwardEta = 0

    # Use lapse process in the simulation
    uselapse = True

    # Number of samples in total
    ntotalsamps = nchains*nsamps

    #Which model?
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
    elif filename.find('PrePostRecLinearDrift') != -1:
        model = 14
        modelname = 'PrePostRecLinearDrift'


    #Load original data
    dataname = filename[0:filename.find('_')]
    print(f'Loading {dataname}...')
    if model < 4:
        data = pd.read_csv(f'../data/{dataname}.csv')
        correctRT = (data['Correct']*2 - 1)*data['ReactionTime']
        savecorrectRT = np.array(np.copy(correctRT))
        rightwardRT = (data['Rightward']*2 - 1)*data['ReactionTime']
        condition = np.array(data['Condition'],dtype='int16') - 1
        session = np.array(data['File'],dtype='int16') - 1
        nconditions = np.max(condition) + 1
        nsessions = np.max(session) + 1
        maxrt = np.empty((nsessions))
        for s in np.unique(session):
            maxrt[s] = np.max(data['ReactionTime'][session==s])
    else:
        data = sio.loadmat(f'../data/{dataname}.mat')
        allrightwardRT = ((data['RF_choice']*2 - 3)*data['responsetime']).squeeze()
        rightwardRT = allrightwardRT[data['training'][0,:]==1]
        if (model==4) | ((model > 7) & (model < 10)):
            allcorrectRT = ((data['correct']*2 - 1)*data['responsetime']).squeeze()
            correctRT = correctRT[data['training'][0,:]==1]
            savecorrectRT = np.array(np.copy(correctRT))
        allinjectlarge = np.asarray(data['simpleinject'],dtype='int16').squeeze()
        injectlarge = allinjectlarge[data['training'][0,:]==1]
        allcondition = np.asarray(data['condition'],dtype='int16').squeeze() - 1
        condition = allcondition[data['training'][0,:]==1]
        allsession = np.asarray(data['sessions'],dtype='int16').squeeze() - 1
        session = allsession[data['training'][0,:]==1]
        uniquedate, allday = np.unique(data['date'], return_inverse=True)
        day = allday[data['training'][0,:]==1]
        nconditions = np.max(condition) + 1
        nsessions = np.max(session) + 1
        ndays = np.max(day) + 1

    saverightwardRT = np.array(np.copy(rightwardRT))

    N = np.size(rightwardRT)

    #Set the random seed for replicability 
    random.seed(seed)

    #Generate all N trials Niter times for a total of N*Niter samples
    print(f'Generating {N} * {Niter} = {N*Niter} samples...')
    if (model < 5) | ((model > 7) & (model < 10)):
        correctRT_pred = np.empty((N, Niter))
    rightwardRT_pred = np.empty((N, Niter))
    for n in range(N):
        # Indexing
        randsamps = np.random.randint(ntotalsamps,size=Niter)
        whichchain = np.array(np.floor(randsamps/nsamps),dtype=int)
        whichsamp = np.mod(randsamps,nsamps)
        # Find which samples are from the lapse process and which are from the diffusion process
        if model == 1:
            problapse = samples['problapse'][session[n],whichsamp,whichchain]
        elif (model==2) | (model==3):
            problapse = samples['problapse'][0, whichsamp,whichchain]
        elif (model==4) | (model==5) | (model==13):
            problapse = samples['problapse'][session[n], whichsamp,whichchain]
        elif (model>5):
            problapse = samples['problapse'][day[n], whichsamp,whichchain]
        for i in range(Niter):
            if model == 1:
                Alpha = samples['alpha'][session[n],whichsamp[i],whichchain[i]]
                Ter = samples['ter'][session[n],whichsamp[i],whichchain[i]]
                Deltaacc = samples['deltaacc'][session[n],condition[n],whichsamp[i],whichchain[i]]
                Deltaright = samples['deltaright'][session[n],condition[n],whichsamp[i],whichchain[i]]
                Zeta = samples['beta'][session[n],whichsamp[i],whichchain[i]]*samples['alpha'][session[n],whichsamp[i],whichchain[i]]
            elif model ==2:
                Alpha = samples['alpha'][0,whichsamp[i],whichchain[i]]
                Ter = samples['ter'][0,whichsamp[i],whichchain[i]]
                Deltaacc = samples['deltaacc'][session[n],condition[n],whichsamp[i],whichchain[i]]
                Deltaright = samples['deltaright'][session[n],condition[n],whichsamp[i],whichchain[i]]
                Zeta = samples['beta'][0,whichsamp[i],whichchain[i]]*samples['alpha'][0,whichsamp[i],whichchain[i]]
            elif model ==3:
                Alpha = samples['alpha'][0,whichsamp[i],whichchain[i]]
                Ter = samples['ter'][0,whichsamp[i],whichchain[i]]
                Deltaacc = samples['deltaacc'][condition[n],whichsamp[i],whichchain[i]]
                Deltaright = samples['deltaright'][condition[n],whichsamp[i],whichchain[i]]
                Zeta = samples['beta'][session[n],whichsamp[i],whichchain[i]]*samples['alpha'][0,whichsamp[i],whichchain[i]]
            elif model ==4:
                Alpha = samples['alpha'][session[n],whichsamp[i],whichchain[i]]
                Ter = samples['ter'][session[n],whichsamp[i],whichchain[i]]
                Deltaacc = samples['deltaacc'][session[n],condition[n],whichsamp[i],whichchain[i]]
                Deltaright = samples['deltaright'][session[n],condition[n],whichsamp[i],whichchain[i]]
                Zeta = samples['beta'][session[n],whichsamp[i],whichchain[i]]*samples['alpha'][session[n],whichsamp[i],whichchain[i]]
            elif ((model ==5) | (model==13) | (model==14)): #Note that deltaright in models 'PrePostRec2Accum' and 'PrePostRecLinearDrift' is directly derived from other variables and is tracked
                Alpha = samples['alpha'][session[n],whichsamp[i],whichchain[i]]
                Ter = samples['ter'][session[n],whichsamp[i],whichchain[i]]
                Deltaright = samples['deltaright'][session[n],condition[n],whichsamp[i],whichchain[i]]
                Zeta = samples['beta'][session[n],whichsamp[i],whichchain[i]]*samples['alpha'][session[n],whichsamp[i],whichchain[i]]
            elif model ==6:
                Alpha = samples['alpha'][day[n],whichsamp[i],whichchain[i]]
                Ter = samples['ter'][day[n],whichsamp[i],whichchain[i]]
                Deltaright = samples['deltaright'][session[n],condition[n],whichsamp[i],whichchain[i]]
                Zeta = samples['beta'][day[n],whichsamp[i],whichchain[i]]*samples['alpha'][day[n],whichsamp[i],whichchain[i]]
            elif model ==7:
                Alpha = samples['alpha'][day[n],whichsamp[i],whichchain[i]]
                Ter = samples['ter'][day[n],whichsamp[i],whichchain[i]]
                Deltaright = samples['deltaright'][day[n],condition[n],whichsamp[i],whichchain[i]]
                Zeta = samples['beta'][session[n],whichsamp[i],whichchain[i]]*samples['alpha'][day[n],whichsamp[i],whichchain[i]]
            elif model ==8:
                Alpha = samples['alpha'][day[n],whichsamp[i],whichchain[i]]
                Ter = samples['ter'][day[n],whichsamp[i],whichchain[i]]
                Deltaacc = samples['deltaacc'][session[n],condition[n],whichsamp[i],whichchain[i]]
                Deltaright = samples['deltaright'][session[n],condition[n],whichsamp[i],whichchain[i]]
                Zeta = samples['beta'][day[n],whichsamp[i],whichchain[i]]*samples['alpha'][day[n],whichsamp[i],whichchain[i]]
            elif model ==9:
                Alpha = samples['alpha'][day[n],whichsamp[i],whichchain[i]]
                Ter = samples['ter'][day[n],whichsamp[i],whichchain[i]]
                Deltaacc = samples['deltaacc'][day[n],condition[n],whichsamp[i],whichchain[i]]
                Deltaright = samples['deltaright'][day[n],condition[n],whichsamp[i],whichchain[i]]
                Zeta = samples['beta'][session[n],whichsamp[i],whichchain[i]]*samples['alpha'][day[n],whichsamp[i],whichchain[i]]
            elif model==10:
                Alpha = samples['alpha'][day[n],whichsamp[i],whichchain[i]]
                Ter = samples['ter'][session[n],condition[n],whichsamp[i],whichchain[i]]
                Deltaright = samples['deltaright'][day[n],condition[n],whichsamp[i],whichchain[i]]
                Zeta = samples['beta'][day[n],whichsamp[i],whichchain[i]]*samples['alpha'][day[n],whichsamp[i],whichchain[i]]
            elif model==11:
                Alpha = samples['alpha'][session[n],whichsamp[i],whichchain[i]]
                Ter = samples['ter'][day[n],whichsamp[i],whichchain[i]]
                Deltaright = samples['deltaright'][day[n],condition[n],whichsamp[i],whichchain[i]]
                Zeta = samples['beta'][session[n],whichsamp[i],whichchain[i]]*samples['alpha'][session[n],whichsamp[i],whichchain[i]]
                



            

            notlapse = (problapse[i] < np.random.uniform())
            if (notlapse or not uselapse):
                #Even though these are not independent processes, treat them as independent processes for easy of model evaluation using in-sample prediction
                if (model < 5) | ((model > 7) & (model < 10)):
                    correctRT_pred[n, i] = simuldiff(N=1, Alpha=Alpha,Ter=Ter, Zeta=None, Nu=Deltaacc,Eta=correctEta)
                rightwardRT_pred[n, i] = simuldiff(N=1, Alpha=Alpha, Ter=Ter, Zeta=Zeta, Nu=Deltaright,Eta=rightwardEta)
            else:
                if model < 4:
                    correctRT_pred[n, i] = np.random.uniform(-maxrt[session[n]], maxrt[session[n]])
                    rightwardRT_pred[n, i] = np.random.uniform(-maxrt[session[n]], maxrt[session[n]])
                elif (model < 5) | ((model > 7) & (model < 10)):
                    correctRT_pred[n, i] = np.random.uniform(-3, 3)
                    rightwardRT_pred[n, i] = np.random.uniform(-3, 3)
                else:
                    rightwardRT_pred[n, i] = np.random.uniform(-3, 3)


    if model < 5:
        inpred = dict(correctRT=savecorrectRT, rightwardRT=saverightwardRT, session=np.squeeze(session), condition=np.squeeze(condition), 
            correctRT_pred=correctRT_pred, rightwardRT_pred=rightwardRT_pred, correctEta=correctEta, rightwardEta=rightwardEta, uselapse=uselapse)
    elif (model==5) | (model==13) | (model==14):
        inpred = dict(correctRT=saverightwardRT, rightwardRT=saverightwardRT, session=np.squeeze(session), condition=np.squeeze(condition), 
            correctRT_pred=rightwardRT_pred, rightwardRT_pred=rightwardRT_pred, correctEta=rightwardEta, rightwardEta=rightwardEta, uselapse=uselapse)
    elif (model==6) | (model==7) | (model==10) | (model==11):
        inpred = dict(correctRT=saverightwardRT, rightwardRT=saverightwardRT, session=np.squeeze(session), day=np.squeeze(day), condition=np.squeeze(condition), 
            correctRT_pred=rightwardRT_pred, rightwardRT_pred=rightwardRT_pred, correctEta=rightwardEta, rightwardEta=rightwardEta, uselapse=uselapse)
    elif (model==8) | (model==9):
        inpred = dict(correctRT=savecorrectRT, rightwardRT=saverightwardRT, session=np.squeeze(session), day=np.squeeze(day), condition=np.squeeze(condition), 
            correctRT_pred=correctRT_pred, rightwardRT_pred=rightwardRT_pred, correctEta=rightwardEta, rightwardEta=rightwardEta, uselapse=uselapse)
    savestring = filename + "_inpred.mat"
    print('Saving results to: \n %s' % (savestring))
    sio.savemat(fileloc + savestring, inpred)




if __name__ == '__main__':
    #When mky3_modelresults.py is run as a script, do this
    generate_data(filename=sys.argv[1],seed=11)
    mky3_insample.eval_insample(filename=sys.argv[1])
    mky3_insample.plot_psychometric(filename=sys.argv[1])