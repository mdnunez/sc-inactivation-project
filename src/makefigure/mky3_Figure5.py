# mky3_Figure5.py - Generates modeling results plots in subfigures
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
# 05/21/20     Michael Nunez                     Converted from mky3_modelingsubplots.py
# 05/26/20     Michael Nunez          Use hierarchical parameters for the first column
# 05/27/20     Michael Nunez                 Fix condition index for 0% coherence
# 07/20/20     Michael Nunez          Change colors and drift-rate vector line thickness
# 07/21/20     MIchael Nunez              Change subplot spacing after making 4 columns

#References:
#http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=4
#https://stackoverflow.com/questions/28572731/matlab-ksdensity-equivalent-in-python
#https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.axes.Axes.arrow.html
#https://stackoverflow.com/questions/35781612/matplotlib-control-which-plot-is-on-top
#https://matplotlib.org/api/_as_gen/matplotlib.pyplot.subplots.html
#https://stackoverflow.com/questions/11244514/modify-tick-label-text
#https://stackoverflow.com/questions/35692507/plot-several-image-files-in-matplotlib-subplots/35692695
#https://stackoverflow.com/questions/14708695/specify-figure-size-in-centimeter-in-matplotlib
#https://www.nature.com/nature/for-authors/formatting-guide
#https://stackoverflow.com/questions/20057260/how-to-remove-gaps-between-subplots-in-matplotlib
#https://matplotlib.org/3.1.1/tutorials/text/text_props.html
#https://stackoverflow.com/questions/34212241/control-the-number-of-rows-within-a-legend/34218288
#https://stackoverflow.com/questions/925024/how-can-i-remove-the-top-and-right-axis-in-matplotlib
#https://matplotlib.org/3.2.1/tutorials/intermediate/tight_layout_guide.html
#https://stackoverflow.com/questions/2176424/hiding-axis-text-in-matplotlib-plots

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import scipy.io as sio
import scipy as sp
from scipy import stats
from scipy.optimize import curve_fit
from math import atan2, cos, sin, degrees
from matplotlib import rc

rc('font', **{'family':'serif','serif':['Arial']})
rc('text', usetex=False)

# Set seed
np.random.seed(seed=2021)

# Generate samples from a diffusion process without shortcut to save wiener process

#To simulate a Wiener process use:
#X_t = mu*t + sigma*W_t
#which implies X_t/dt = mu + sigma*(W_t/dt)
#Approximate by X_t(t+dt) = X_t(t) + mu*dt + sigma*sqrt(dt)*randn

#alpha - amount of evidence units necessary to make a decision (evidence units)
#beta - initial evidence bias parameter as a proportion of alpha
#drift - average rate of evidence accumulation within a trial (evidence units per second)
#dcoef - diffusion coefficient, standard deviation of each sample of the evidence accumulation process  (evidence units per second)
#eta - standard deviation of the drift rate distribution across trials
#dt - length of time per sample (seconds)
#n - number of obervations of RT and choice
#maxsample - maximum number of samples allowed in the simulation

def ratcliffDDM(beta,drift,alpha,dcoef,dt,n,maxsample,eta=0) :
    resp = []
    rt = []
    data = []
    allevidence = np.empty((n,maxsample))
    for i in range(n):
        evidence = []
        x=beta*alpha
        sample=0
        thisdrift = drift + eta*np.random.normal()
        while (sample<maxsample):
            sample=sample+1
            evidence.append(x)
            x = x + dt*thisdrift + dcoef*np.sqrt(dt)*np.random.normal()
            if (x>=alpha):
                resp.append(float(1.0))
                break
            if (x<=0):
                resp.append(float(-1.0))
                break
            if (sample==maxsample):
                resp.append(np.nan)
                break
        allevidence[i,0:(sample)] = np.asarray(evidence)
        allevidence[i,(sample):maxsample] = np.asarray(evidence)[-1] #The evidence trace no longer accumulates after the last sample
        #Return the response time after the boundary is crossed with a correction to find the true RT after the crossing
        number=((float(sample))*dt)-(dt/(float(2.0))) 
        rt.append(number)
    for i in range(n):
        temp=resp[i]*rt[i]
        data.append(temp)
    alldata = np.asarray(data)
    return (alldata, allevidence)


def fpsychometric(x, gamma, lam, beta, alpha):
        return (gamma + (1 - gamma - lam)) / (1.0 + np.exp(-beta*(x-alpha))) #From Shen and Richards, 2012, Acoustical Society of America


def mm2inch(*tupl):
    mmperinch = 25.4
    if isinstance(tupl[0], tuple):
        return tuple(i/mmperinch for i in tupl[0])
    else:
        return tuple(i/mmperinch for i in tupl)


#Make the figure
fontsize = 7
panelsize = 8
markersize = 2
legendsize = 5
precolor = 'k'
precolorlite = '#7e7e7e'
postcolor = '#e86800'
postcolorlite = '#ffb160'
reccolor = '#008080'
reccolorlite = '#80b3b3'


maxx = 2250 #Maximum number of ms to plot


dt = .001 #Seconds
lastsec = 2
time = np.arange(np.round(lastsec/float(dt)))*dt*1000


#Plot the simulated psychometric curves
evidenceright = np.array([-36, -24, -17, -10, -5, -3, 0, 3, 5, 10, 17, 24, 36])
psychometric_xtick = np.array([-36, -17, -5, 0, 5, 17, 36])


fig, axs = plt.subplots(2, 4,figsize=mm2inch(183,73), sharex='col')



##### Load model results #####
# plt.figure()
filename = 'fixedexpSCInjectGPRTdata_PrePostRecHalfApr_10_20_15_50'
dataname = filename[0:filename.find('_')]

modelloc = (f'../modelfits/{filename}.mat')
print(f'Loading model fit {filename}...')
samples = sio.loadmat(modelloc)
deltarightsamps = samples['deltaright']
nchains = deltarightsamps.shape[-1]
nsamps = deltarightsamps.shape[-2]
nses = deltarightsamps.shape[0]
nconds = deltarightsamps.shape[1]

print(f'Loading original data {dataname}...')
truedata = sio.loadmat(f'../data/{dataname}.mat')
responsetime = np.asarray(truedata['responsetime']).squeeze()
RF_choice = np.asarray(truedata['RF_choice']).squeeze()
rightwardRT = (RF_choice*2 - 3)*responsetime
target = np.asarray(truedata['RF_target'],dtype='int16').squeeze()
injectlarge = np.asarray(truedata['simpleinject'],dtype='int16').squeeze()
exactinjectlarge = np.asarray(truedata['injectioncond'],dtype='int16').squeeze()
session = np.asarray(truedata['sessions'],dtype='int16').squeeze() - 1
monkeylarge = np.asarray(truedata['monkey'],dtype='int16').squeeze()
condition = np.asarray(truedata['condition'],dtype='int16').squeeze() - 1
monkey = np.empty((nses))
injectfull = np.empty((nses))
exactinject = np.empty((nses))
for s in np.unique(session):
    monkey[s] = monkeylarge[session == s][0]
    injectfull[s] = injectlarge[session == s][0]
    exactinject[s] = exactinjectlarge[session == s][0]
monkey = np.asarray(monkey, dtype='int16')
injectfull = np.asarray(injectfull, dtype='int16')
exactinject = np.asarray(exactinject, dtype='int16')

pre_drift_SP = np.zeros((nconds))
post_drift_SP = np.zeros((nconds))
pre_beta_SP = np.zeros((nconds))
post_beta_SP = np.zeros((nconds))
pre_ndt_SP = np.zeros((nconds))
post_ndt_SP = np.zeros((nconds))
pre_alpha_SP = np.zeros((nconds))
post_alpha_SP = np.zeros((nconds))
pre_lapse_SP = np.zeros((nconds))
post_lapse_SP = np.zeros((nconds))

pre_drift_BT = np.zeros((nconds))
post_drift_BT = np.zeros((nconds))
pre_beta_BT = np.zeros((nconds))
post_beta_BT = np.zeros((nconds))
pre_ndt_BT = np.zeros((nconds))
post_ndt_BT = np.zeros((nconds))
pre_alpha_BT = np.zeros((nconds))
post_alpha_BT = np.zeros((nconds))
pre_lapse_BT = np.zeros((nconds))
post_lapse_BT = np.zeros((nconds))
for c in range(nconds):
    if (c >= 7):
        thisindex = c
    else:
        thisindex = (6-c)
    pre_drift_SP[thisindex] = np.median(deltarightsamps[((monkey==1) & (injectfull==1)),c,:,:])
    post_drift_SP[thisindex] = np.median(deltarightsamps[((monkey==1) & (injectfull==2)),c,:,:])
    pre_drift_BT[thisindex] = np.median(deltarightsamps[((monkey==2) & (injectfull==1)),c,:,:])
    post_drift_BT[thisindex] = np.median(deltarightsamps[((monkey==2) & (injectfull==2)),c,:,:])

pre_beta_SP = np.median(samples['betahier'][0,0,:,:])
print('Pre- hierarchical beta for monkey SP was %.3f' % (pre_beta_SP))
post_beta_SP = np.median(samples['betahier'][1,0,:,:])
print('Post- hierarchical beta for monkey SP was %.3f' % (post_beta_SP))
pre_ndt_SP = np.median(samples['terhier'][0,0,:,:])*1000
print('Pre- hierarchical non-decision time for monkey SP was %.3f' % (pre_ndt_SP))
post_ndt_SP = np.median(samples['terhier'][1,0,:,:])*1000
print('Post- hierarchical non-decision time for monkey SP was %.3f' % (post_ndt_SP))
pre_alpha_SP = np.median(samples['alphahier'][0,0,:,:])
print('Pre- hierarchical symmetric boundary for monkey SP was %.3f' % (pre_alpha_SP))
post_alpha_SP = np.median(samples['alphahier'][1,0,:,:])
print('Post- hierarchical symmetric boundary for monkey SP was %.3f' % (post_alpha_SP))
pre_lapse_SP = np.median(samples['problapsehier'][0,0,:,:])
print('Pre- hierarchical lapse proportion for monkey SP was %.3f' % (pre_lapse_SP))
post_lapse_SP = np.median(samples['problapsehier'][1,0,:,:])
print('Post- hierarchical lapse proportion for monkey SP was %.3f' % (post_lapse_SP))
pre_Delta_SP = np.median(samples['deltarighthier'][0,0,:,:])
print('Pre- hierarchical drift-criterion for monkey SP was %.3f' % (pre_Delta_SP))
post_Delta_SP = np.median(samples['deltarighthier'][1,0,:,:])
print('Post- hierarchical drift-criterion for monkey SP was %.3f' % (post_Delta_SP))


pre_beta_BT = np.median(samples['betahier'][0,1,:,:])
print('Pre- hierarchical beta for monkey BT was %.3f' % (pre_beta_BT))
post_beta_BT = np.median(samples['betahier'][1,1,:,:])
print('Post- hierarchical beta for monkey BT was %.3f' % (post_beta_BT))
pre_ndt_BT = np.median(samples['terhier'][0,1,:,:])*1000
print('Pre- hierarchical non-decision time for monkey BT was %.3f' % (pre_ndt_BT))
post_ndt_BT = np.median(samples['terhier'][1,1,:,:])*1000
print('Post- hierarchical non-decision time for monkey BT was %.3f' % (post_ndt_BT))
pre_alpha_BT = np.median(samples['alphahier'][0,1,:,:])
print('Pre- hierarchical symmetric boundary for monkey BT was %.3f' % (pre_alpha_BT))
post_alpha_BT = np.median(samples['alphahier'][1,1,:,:])
print('Post- hierarchical symmetric boundary for monkey BT was %.3f' % (post_alpha_BT))
pre_lapse_BT = np.median(samples['problapsehier'][0,1,:,:])
print('Pre- hierarchical lapse proportion for monkey BT was %.3f' % (pre_lapse_BT))
post_lapse_BT = np.median(samples['problapsehier'][1,1,:,:])
print('Post- hierarchical lapse proportion for monkey BT was %.3f' % (post_lapse_BT))
pre_Delta_BT = np.median(samples['deltarighthier'][0,1,:,:])
print('Pre- hierarchical drift-criterion for monkey BT was %.3f' % (pre_Delta_BT))
post_Delta_BT = np.median(samples['deltarighthier'][1,1,:,:])
print('Post- hierarchical drift-criterion for monkey BT was %.3f' % (post_Delta_BT))



#Plot boundaries
axs[0,0].plot(np.array([0, maxx]), np.array([-pre_alpha_SP/2, -pre_alpha_SP/2]),linewidth=1,color=precolor,linestyle='-', zorder=13)
axs[0,0].plot(np.array([0, maxx]), np.array([-post_alpha_SP/2, -post_alpha_SP/2]),linewidth=1,color=postcolor,linestyle='-', zorder=14)
axs[0,0].plot(np.array([0, maxx]), np.array([pre_alpha_SP/2, pre_alpha_SP/2]),linewidth=1,color=precolor,linestyle='-', zorder=18)
axs[0,0].plot(np.array([0, maxx]), np.array([post_alpha_SP/2, post_alpha_SP/2]),linewidth=1,color=postcolor,linestyle='-', zorder=19)

#Plot non-decision times
axs[0,0].plot(np.array([0, pre_ndt_SP]), np.array([pre_beta_SP-.5, pre_beta_SP-.5]),linewidth=.5,color=precolorlite)
axs[0,0].plot(np.array([0, post_ndt_SP]), np.array([post_beta_SP-.5, post_beta_SP-.5]),linewidth=.5,color=postcolorlite)



angles_radians_pre = np.empty((pre_drift_SP.shape[0]))
angles_radians_post = np.empty((post_drift_SP.shape[0]))
angletrack = 0
arrowlength = 1000
for c in (0, 3, 6, -4, -1):
    angle_pre = atan2( pre_drift_SP[c] , 1000 ) #Find the angle from start point line, drift rate is the average rate of accumulation per second
    angles_radians_pre[angletrack] = angle_pre
    xlen_pre = arrowlength*cos(angle_pre)
    ylen_pre = arrowlength*sin(angle_pre)
    angle_post = atan2( post_drift_SP[c] , 1000 ) #Find the angle from start point line, drift rate is the average rate of accumulation per second
    angles_radians_post[angletrack] = angle_post
    xlen_post = arrowlength*cos(angle_post)
    ylen_post = arrowlength*sin(angle_post)
    angletrack += 1
    axs[0,0].quiver(pre_ndt_SP, pre_beta_SP-.5, xlen_pre, ylen_pre, color=precolor, alpha=1, angles='xy', scale=8000, width=0.005, zorder=20)
    axs[0,0].quiver(post_ndt_SP, post_beta_SP-.5, xlen_post, ylen_post, color=postcolor, alpha=1, angles='xy', scale=8000, width=0.005, zorder=21)

#Plot hierarchical drift rate (mean across sessions and conditions)
arrowlength = 1000
angle_pre = atan2( pre_Delta_SP , 1000 ) #Find the angle from start point line, drift rate is the average rate of accumulation per second
xlen_pre = arrowlength*cos(angle_pre)
ylen_pre = arrowlength*sin(angle_pre)
angle_post = atan2( post_Delta_SP , 1000 ) #Find the angle from start point line, drift rate is the average rate of accumulation per second
xlen_post = arrowlength*cos(angle_post)
ylen_post = arrowlength*sin(angle_post)
angletrack += 1
axs[0,0].quiver(pre_ndt_SP, pre_beta_SP-.5, xlen_pre, ylen_pre, color=precolor, alpha=1, angles='xy', scale=5000, zorder=20)
axs[0,0].quiver(post_ndt_SP, post_beta_SP-.5, xlen_post, ylen_post, color=postcolor, alpha=1, angles='xy', scale=5000, zorder=20)


#Plot and shade RT distributions
allRT_pre = responsetime[((injectlarge==1) & (monkeylarge==1) & (condition==0))]
allchoice_pre = RF_choice[((injectlarge==1) & (monkeylarge==1) & (condition==0))]
alldata_pre = (allchoice_pre*2 - 3)*allRT_pre
kde_upper_pre = stats.gaussian_kde(alldata_pre[alldata_pre>0]*1000)
axs[0,0].plot(time, kde_upper_pre(time)*400 + pre_alpha_SP/2, linewidth=1,color=precolor,linestyle='-', zorder=9)
kde_lower_pre = stats.gaussian_kde(np.abs(alldata_pre[alldata_pre<0])*1000)
axs[0,0].plot(time, kde_lower_pre(time)*-400 - pre_alpha_SP/2, linewidth=1,color=precolor,linestyle='-', zorder=10)
axs[0,0].fill_between(time, (pre_alpha_SP/2)*np.ones((time.shape[0])), kde_upper_pre(time)*400 + pre_alpha_SP/2, color=precolorlite, alpha=0.5, zorder=5)
axs[0,0].fill_between(time, -(pre_alpha_SP/2)*np.ones((time.shape[0])), -kde_lower_pre(time)*400 - pre_alpha_SP/2, color=precolorlite, alpha=0.5, zorder=6)

allRT_post = responsetime[((injectlarge==2) & (monkeylarge==1) & (condition==0))]
allchoice_post = RF_choice[((injectlarge==2) & (monkeylarge==1) & (condition==0))]
alldata_post = (allchoice_post*2 - 3)*allRT_post
kde_upper_post = stats.gaussian_kde(alldata_post[alldata_post>0]*1000)
axs[0,0].plot(time, kde_upper_post(time)*400 + post_alpha_SP/2, linewidth=1,color=postcolor,linestyle='-', zorder=11)
kde_lower_post = stats.gaussian_kde(np.abs(alldata_post[alldata_post<0])*1000)
axs[0,0].plot(time, -kde_lower_post(time)*400 - post_alpha_SP/2, linewidth=1,color=postcolor,linestyle='-', zorder=12)
axs[0,0].fill_between(time, (post_alpha_SP/2)*np.ones((time.shape[0])), kde_upper_post(time)*400 + post_alpha_SP/2, color=postcolorlite, alpha=0.5, zorder=7)
axs[0,0].fill_between(time, -(post_alpha_SP/2)*np.ones((time.shape[0])), -kde_lower_post(time)*400 - post_alpha_SP/2, color=postcolorlite, alpha=0.5, zorder=8)

#Appropriately label plots
# axs[0,0].set_xlabel('Time after Glass Pattern onset (ms)', fontsize=fontsize)
# axs[0,0].set_ylabel('Relative evidence for a choice to the IF', fontsize=fontsize)
axs[0,0].set_xlim([250,maxx]) #X limits in ms
axs[0,0].set_yticks(np.array([-1.0, 0, 1.0]))
custom_lines = [Line2D([0], [0], color=precolor, linestyle='-', lw=1),
                Line2D([0], [0], color=postcolor, lw=1)]
axs[0,0].legend(custom_lines, ['Monkey S Pre-','Monkey S Post-'],  loc=7, fontsize=legendsize)
axs[0,0].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[0,0].spines['right'].set_visible(False)
axs[0,0].spines['top'].set_visible(False)
axs[0,0].get_xaxis().set_visible(False)
axs[0,0].yaxis.set_ticks_position('left')
axs[0,0].xaxis.set_ticks_position('bottom')
axs[0,0].set_ylabel('Relative evidence for a choice to the IF', fontsize=fontsize)
# axs[0,0].get_yaxis().set_label_coords(-0.3,0.5)
axs[0,0].text(0.02, .98, 'q', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[0,0].transAxes)
# plt.savefig((f'../../figures/diffusion_figures/Diffusion_SP_0coh.png'), dpi=300, format='png',bbox_inches="tight")


#Plot boundaries
axs[1,0].plot(np.array([0, maxx]), np.array([-pre_alpha_BT/2, -pre_alpha_BT/2]),linewidth=1,color=precolor,linestyle='-', zorder=13)
axs[1,0].plot(np.array([0, maxx]), np.array([-post_alpha_BT/2, -post_alpha_BT/2]),linewidth=1,color=postcolor,linestyle='-', zorder=14)
axs[1,0].plot(np.array([0, maxx]), np.array([pre_alpha_BT/2, pre_alpha_BT/2]),linewidth=1,color=precolor,linestyle='-', zorder=18)
axs[1,0].plot(np.array([0, maxx]), np.array([post_alpha_BT/2, post_alpha_BT/2]),linewidth=1,color=postcolor,linestyle='-', zorder=19)

#Plot non-decision times
axs[1,0].plot(np.array([0, pre_ndt_BT]), np.array([pre_beta_BT-.5, pre_beta_BT-.5]),linewidth=.5,color=precolorlite)
axs[1,0].plot(np.array([0, post_ndt_BT]), np.array([post_beta_BT-.5, post_beta_BT-.5]),linewidth=.5,color=postcolorlite)


angles_radians_pre = np.empty((pre_drift_BT.shape[0]))
angles_radians_post = np.empty((post_drift_BT.shape[0]))
angletrack = 0
arrowlength = 1000
for c in (0, 3, 6, -4, -1):
    angle_pre = atan2( pre_drift_BT[c] , 1000 ) #Find the angle from start point line, drift rate is the average rate of accumulation per second
    angles_radians_pre[angletrack] = angle_pre
    xlen_pre = arrowlength*cos(angle_pre)
    ylen_pre = arrowlength*sin(angle_pre)
    angle_post = atan2( post_drift_BT[c] , 1000 ) #Find the angle from start point line, drift rate is the average rate of accumulation per second
    angles_radians_post[angletrack] = angle_post
    xlen_post = arrowlength*cos(angle_post)
    ylen_post = arrowlength*sin(angle_post)
    angletrack += 1
    axs[1,0].quiver(pre_ndt_BT, pre_beta_BT-.5, xlen_pre, ylen_pre, color=precolor, alpha=1, angles='xy', scale=8000, width=0.005, zorder=20)
    axs[1,0].quiver(post_ndt_BT, post_beta_BT-.5, xlen_post, ylen_post, color=postcolor, alpha=1, angles='xy', scale=8000, width=0.005, zorder=21)

#Plot hierarchical drift rate (mean across sessions and conditions)
arrowlength = 1000
angle_pre = atan2( pre_Delta_BT , 1000 ) #Find the angle from start point line, drift rate is the average rate of accumulation per second
xlen_pre = arrowlength*cos(angle_pre)
ylen_pre = arrowlength*sin(angle_pre)
angle_post = atan2( post_Delta_BT , 1000 ) #Find the angle from start point line, drift rate is the average rate of accumulation per second
xlen_post = arrowlength*cos(angle_post)
ylen_post = arrowlength*sin(angle_post)
angletrack += 1
axs[1,0].quiver(pre_ndt_BT, pre_beta_BT-.5, xlen_pre, ylen_pre, color=precolor, alpha=1, angles='xy', scale=5000, zorder=20)
axs[1,0].quiver(post_ndt_BT, post_beta_BT-.5, xlen_post, ylen_post, color=postcolor, alpha=1, angles='xy', scale=5000, zorder=20)


#Plot and shade RT distributions
allRT_pre = responsetime[((exactinjectlarge==1) & (monkeylarge==2) & (condition==0))]
allchoice_pre = RF_choice[((exactinjectlarge==1) & (monkeylarge==2) & (condition==0))]
alldata_pre = (allchoice_pre*2 - 3)*allRT_pre
kde_upper_pre = stats.gaussian_kde(alldata_pre[alldata_pre>0]*1000)
axs[1,0].plot(time, kde_upper_pre(time)*400 + pre_alpha_BT/2, linewidth=1,color=precolor,linestyle='-', zorder=9)
kde_lower_pre = stats.gaussian_kde(np.abs(alldata_pre[alldata_pre<0])*1000)
axs[1,0].plot(time, kde_lower_pre(time)*-400 - pre_alpha_BT/2, linewidth=1,color=precolor,linestyle='-', zorder=10)
axs[1,0].fill_between(time, (pre_alpha_BT/2)*np.ones((time.shape[0])), kde_upper_pre(time)*400 + pre_alpha_BT/2, color=precolorlite, alpha=0.5, zorder=5)
axs[1,0].fill_between(time, -(pre_alpha_BT/2)*np.ones((time.shape[0])), -kde_lower_pre(time)*400 - pre_alpha_BT/2, color=precolorlite, alpha=0.5, zorder=6)

allRT_post = responsetime[((exactinjectlarge==2) & (monkeylarge==2) & (condition==0))]
allchoice_post = RF_choice[((exactinjectlarge==2) & (monkeylarge==2) & (condition==0))]
alldata_post = (allchoice_post*2 - 3)*allRT_post
kde_upper_post = stats.gaussian_kde(alldata_post[alldata_post>0]*1000)
axs[1,0].plot(time, kde_upper_post(time)*400 + post_alpha_BT/2, linewidth=1,color=postcolor,linestyle='-', zorder=11)
kde_lower_post = stats.gaussian_kde(np.abs(alldata_post[alldata_post<0])*1000)
axs[1,0].plot(time, -kde_lower_post(time)*400 - post_alpha_BT/2, linewidth=1,color=postcolor,linestyle='-', zorder=12)
axs[1,0].fill_between(time, (post_alpha_BT/2)*np.ones((time.shape[0])), kde_upper_post(time)*400 + post_alpha_BT/2, color=postcolorlite, alpha=0.5, zorder=7)
axs[1,0].fill_between(time, -(post_alpha_BT/2)*np.ones((time.shape[0])), -kde_lower_post(time)*400 - post_alpha_BT/2, color=postcolorlite, alpha=0.5, zorder=8)

#Appropriately label plots
axs[1,0].set_xlabel('Time after Glass Pattern onset (ms)', fontsize=fontsize)
# axs[1,0].set_ylabel('Relative evidence for a choice to the IF', fontsize=fontsize)
axs[1,0].set_xlim([250,maxx]) #X limits in ms
axs[1,0].set_yticks(np.array([-1.0, 0, 1.0]))
# plt.tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
custom_lines = [Line2D([0], [0], color=precolor, linestyle='-', lw=1),
                Line2D([0], [0], color=postcolor, lw=1)]
axs[1,0].legend(custom_lines, ['Monkey B Pre-','Monkey B Post-'],  loc=7, fontsize=legendsize)
axs[1,0].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[1,0].spines['right'].set_visible(False)
axs[1,0].spines['top'].set_visible(False)
axs[1,0].yaxis.set_ticks_position('left')
axs[1,0].xaxis.set_ticks_position('bottom')
axs[1,0].text(0.02, .98, 'u', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[1,0].transAxes)
# plt.savefig((f'../../figures/diffusion_figures/Diffusion_BT_0coh.png'), dpi=300, format='png',bbox_inches="tight")

print(f'Loading {filename} in-sample prediction...')
fileloc = '../modelfits/'
inpredloc = fileloc + filename + '_inpred.mat'
inpred = sio.loadmat(inpredloc)

axs[0,1].plot(np.array([-50, 50]), np.array([50, 50]), color='k', LineStyle='--')
axs[0,1].plot(np.array([0, 0]), np.array([0, 100]), color='k', LineStyle='--')
axs[1,1].plot(np.array([-50, 50]), np.array([50, 50]), color='k', LineStyle='--')
axs[1,1].plot(np.array([0, 0]), np.array([0, 100]), color='k', LineStyle='--')

psychometric_true = np.ones((nses, nconds))*np.nan
psychometric_pred = np.ones((nses, nconds))*np.nan
meanRT_correct_true = np.ones((nses, nconds))*np.nan
meanRT_incorrect_true = np.ones((nses, nconds))*np.nan
meanRT_correct_pred = np.ones((nses, nconds))*np.nan
meanRT_incorrect_pred = np.ones((nses, nconds))*np.nan

meanRT_correct_true_SP_pre = np.ones((nses, nconds))*np.nan
meanRT_incorrect_true_SP_pre = np.ones((nses, nconds))*np.nan
meanRT_correct_pred_SP_pre = np.ones((nses, nconds))*np.nan
meanRT_incorrect_pred_SP_pre = np.ones((nses, nconds))*np.nan
SP_pre_index = 0

meanRT_correct_true_SP_post = np.ones((nses, nconds))*np.nan
meanRT_incorrect_true_SP_post = np.ones((nses, nconds))*np.nan
meanRT_correct_pred_SP_post = np.ones((nses, nconds))*np.nan
meanRT_incorrect_pred_SP_post = np.ones((nses, nconds))*np.nan
SP_post_index = 0

meanRT_correct_true_BT_pre = np.ones((nses, nconds))*np.nan
meanRT_incorrect_true_BT_pre = np.ones((nses, nconds))*np.nan
meanRT_correct_pred_BT_pre = np.ones((nses, nconds))*np.nan
meanRT_incorrect_pred_BT_pre = np.ones((nses, nconds))*np.nan
BT_pre_index = 0

meanRT_correct_true_BT_post = np.ones((nses, nconds))*np.nan
meanRT_incorrect_true_BT_post = np.ones((nses, nconds))*np.nan
meanRT_correct_pred_BT_post = np.ones((nses, nconds))*np.nan
meanRT_incorrect_pred_BT_post = np.ones((nses, nconds))*np.nan
BT_post_index = 0

for s in range(nses):
    for c in range(nconds):
        wherecond = (condition == c)
        wheresession = (session == s)
        whereboth = (wherecond) & (wheresession)
        thisinjectioncond = exactinjectlarge[whereboth][0]
        thismonkey = monkeylarge[whereboth][0]
        thissimpleinject = injectlarge[whereboth][0]
        predcond = (np.squeeze(inpred['condition']) == c)
        predsession = (np.squeeze(inpred['session']) == s)
        predboth = (predcond) & (predsession)
        if thissimpleinject > 2:
            break

        chose_right_when_right_true = (rightwardRT[whereboth] > 0) & (target[whereboth] == 2.)
        chose_right_when_left_true = (rightwardRT[whereboth] > 0) & (target[whereboth] == 1.)
        chose_left_when_right_true = (rightwardRT[whereboth] < 0) & (target[whereboth] == 2.)
        chose_left_when_left_true = (rightwardRT[whereboth] < 0) & (target[whereboth] == 1.)

        if (c >= 7):
            psychometric_true[s, c] = ((np.sum(chose_right_when_left_true) + np.sum(chose_right_when_right_true)) / (chose_right_when_right_true.shape[0]))*100
            psychometric_pred[s, c] = np.mean(inpred['rightwardRT_pred'][predboth,:] > 0,axis=None)*100
            meanRT_correct_true[s, c] = np.mean(np.abs(rightwardRT[whereboth][rightwardRT[whereboth] > 0]))
            meanRT_incorrect_true[s, c] = np.mean(np.abs(rightwardRT[whereboth][rightwardRT[whereboth] < 0]))
            allpredRT = inpred['rightwardRT_pred'][predboth,:].flatten()
            meanRT_correct_pred[s, c] = np.mean(np.abs(allpredRT[allpredRT > 0]))
            meanRT_incorrect_pred[s, c] = np.mean(np.abs(allpredRT[allpredRT < 0]))
        else:
            psychometric_true[s, 6-c] = ((np.sum(chose_right_when_left_true) + np.sum(chose_right_when_right_true)) / (chose_right_when_right_true.shape[0]))*100
            psychometric_pred[s, 6-c] = np.mean(inpred['rightwardRT_pred'][predboth,:] > 0,axis=None)*100
            meanRT_correct_true[s, 6-c] = np.mean(np.abs(rightwardRT[whereboth][rightwardRT[whereboth] < 0]))
            meanRT_incorrect_true[s, 6-c] = np.mean(np.abs(rightwardRT[whereboth][rightwardRT[whereboth] > 0]))
            allpredRT = inpred['rightwardRT_pred'][predboth,:].flatten()
            meanRT_correct_pred[s, 6-c] = np.mean(np.abs(allpredRT[allpredRT < 0]))
            meanRT_incorrect_pred[s, 6-c] = np.mean(np.abs(allpredRT[allpredRT > 0]))
    if (thissimpleinject == 1):
        thiscolor = precolor
        thislite = precolorlite
        if (thismonkey == 1):
            meanRT_correct_true_SP_pre[SP_pre_index,:] = meanRT_correct_true[s, :]
            meanRT_correct_pred_SP_pre[SP_pre_index,:] = meanRT_correct_pred[s, :]
            meanRT_incorrect_true_SP_pre[SP_pre_index,:] = meanRT_incorrect_true[s, :]
            meanRT_incorrect_pred_SP_pre[SP_pre_index,:] = meanRT_incorrect_pred[s, :]
            SP_pre_index += 1
        elif (thismonkey == 2):
            meanRT_correct_true_BT_pre[BT_pre_index,:] = meanRT_correct_true[s, :]
            meanRT_correct_pred_BT_pre[BT_pre_index,:] = meanRT_correct_pred[s, :]
            meanRT_incorrect_true_BT_pre[BT_pre_index,:] = meanRT_incorrect_true[s, :]
            meanRT_incorrect_pred_BT_pre[BT_pre_index,:] = meanRT_incorrect_pred[s, :]
            BT_pre_index += 1
    elif (thissimpleinject == 2):
        thiscolor = postcolor
        thislite = postcolorlite
        if (thismonkey == 1):
            meanRT_correct_true_SP_post[SP_post_index,:] = meanRT_correct_true[s, :]
            meanRT_correct_pred_SP_post[SP_post_index,:] = meanRT_correct_pred[s, :]
            meanRT_incorrect_true_SP_post[SP_post_index,:] = meanRT_incorrect_true[s, :]
            meanRT_incorrect_pred_SP_post[SP_post_index,:] = meanRT_incorrect_pred[s, :]
            SP_post_index += 1
        elif (thismonkey == 2):
            meanRT_correct_true_BT_post[BT_post_index,:] = meanRT_correct_true[s, :]
            meanRT_correct_pred_BT_post[BT_post_index,:] = meanRT_correct_pred[s, :]
            meanRT_incorrect_true_BT_post[BT_post_index,:] = meanRT_incorrect_true[s, :]
            meanRT_incorrect_pred_BT_post[BT_post_index,:] = meanRT_incorrect_pred[s, :]
            BT_post_index += 1
    axs[thismonkey-1,1].plot(evidenceright,psychometric_true[s,:], 'o',markersize=markersize, color=thiscolor)
    # axs[thismonkey-1,2].plot(evidenceright,meanRT_correct_true[s,:]*1000, 'o', color=thiscolor)
    # axs[thismonkey-1,2].plot(evidenceright,meanRT_incorrect_true[s,:]*1000, 's', color=thiscolor)
    # axs[thismonkey-1,2].plot(evidenceright,meanRT_correct_pred[s,:]*1000, linestyle='-', color=thislite)
    # axs[thismonkey-1,2].plot(evidenceright,meanRT_incorrect_pred[s,:]*1000, linestyle='--', color=thislite)
    try:
        popt_pred, pcov = curve_fit(fpsychometric, evidenceright, psychometric_pred[s,:]/100)
        axs[thismonkey-1,1].plot(evidenceright, fpsychometric(evidenceright, *popt_pred)*100, color=thislite)
    except:
        print('PMF curve could not be fit for session %d!' % (s))
        pass
    

axs[0,2].plot(evidenceright,np.nanmean(meanRT_correct_true_SP_pre[0:SP_pre_index,:],axis=0)*1000, 'o', markersize=markersize, color=precolor)
axs[0,3].plot(evidenceright,np.nanmean(meanRT_incorrect_true_SP_pre[0:SP_pre_index,:],axis=0)*1000, 's', markersize=markersize, color=precolor)
axs[0,2].plot(evidenceright,np.nanmean(meanRT_correct_pred_SP_pre[0:SP_pre_index,:],axis=0)*1000, linestyle='-', color=precolorlite)
axs[0,3].plot(evidenceright,np.nanmean(meanRT_incorrect_pred_SP_pre[0:SP_pre_index,:],axis=0)*1000, linestyle='--', color=precolorlite)

jitter = 0
axs[0,2].plot(evidenceright,np.nanmean(meanRT_correct_true_SP_post[0:SP_post_index,:],axis=0)*1000, 'o', markersize=markersize, color=postcolor)
axs[0,3].plot(evidenceright+jitter,np.nanmean(meanRT_incorrect_true_SP_post[0:SP_post_index,:],axis=0)*1000, 's', markersize=markersize, color=postcolor)
axs[0,2].plot(evidenceright,np.nanmean(meanRT_correct_pred_SP_post[0:SP_post_index,:],axis=0)*1000, linestyle='-', color=postcolorlite)
axs[0,3].plot(evidenceright+jitter,np.nanmean(meanRT_incorrect_pred_SP_post[0:SP_post_index,:],axis=0)*1000, linestyle='--', color=postcolorlite)

# CI95upper_SP_pre = np.nanmean(meanRT_correct_true_SP_pre[0:SP_pre_index,:],axis=0)*1000 + 1.96*np.nanstd(meanRT_correct_true_SP_pre[0:SP_pre_index,:],axis=0)*1000
# CI95lower_SP_pre = np.nanmean(meanRT_correct_true_SP_pre[0:SP_pre_index,:],axis=0)*1000 - 1.96*np.nanstd(meanRT_correct_true_SP_pre[0:SP_pre_index,:],axis=0)*1000
# axs[0,2].fill_between(evidenceright, CI95upper_SP_pre, CI95lower_SP_pre, color=precolor, alpha=0.25)

# CI95upper_SP_post = np.nanmean(meanRT_correct_true_SP_post[0:SP_post_index,:],axis=0)*1000 + 1.96*np.nanstd(meanRT_correct_true_SP_post[0:SP_post_index,:],axis=0)*1000
# CI95lower_SP_post = np.nanmean(meanRT_correct_true_SP_post[0:SP_post_index,:],axis=0)*1000 - 1.96*np.nanstd(meanRT_correct_true_SP_post[0:SP_post_index,:],axis=0)*1000
# axs[0,2].fill_between(evidenceright, CI95upper_SP_post, CI95lower_SP_post, color=postcolor, alpha=0.25)

axs[1,2].plot(evidenceright,np.nanmean(meanRT_correct_true_BT_pre[0:BT_pre_index,:],axis=0)*1000, 'o', markersize=markersize, color=precolor)
axs[1,3].plot(evidenceright,np.nanmean(meanRT_incorrect_true_BT_pre[0:BT_pre_index,:],axis=0)*1000, 's', markersize=markersize, color=precolor)
axs[1,2].plot(evidenceright,np.nanmean(meanRT_correct_pred_BT_pre[0:BT_pre_index,:],axis=0)*1000, linestyle='-', color=precolorlite)
axs[1,3].plot(evidenceright,np.nanmean(meanRT_incorrect_pred_BT_pre[0:BT_pre_index,:],axis=0)*1000, linestyle='--', color=precolorlite)

jitter = 0
axs[1,2].plot(evidenceright+jitter,np.nanmean(meanRT_correct_true_BT_post[0:BT_post_index,:],axis=0)*1000, 'o', markersize=markersize, color=postcolor)
axs[1,3].plot(evidenceright+jitter,np.nanmean(meanRT_incorrect_true_BT_post[0:BT_post_index,:],axis=0)*1000, 's', markersize=markersize, color=postcolor)
axs[1,2].plot(evidenceright,np.nanmean(meanRT_correct_pred_BT_post[0:BT_post_index,:],axis=0)*1000, linestyle='-', color=postcolorlite)
axs[1,3].plot(evidenceright,np.nanmean(meanRT_incorrect_pred_BT_post[0:BT_post_index,:],axis=0)*1000, linestyle='--', color=postcolorlite)

# CI95upper_BT_pre = np.nanmean(meanRT_correct_true_BT_pre[0:BT_pre_index,:],axis=0)*1000 + 1.96*np.nanstd(meanRT_correct_true_BT_pre[0:BT_pre_index,:],axis=0)*1000
# CI95lower_BT_pre = np.nanmean(meanRT_correct_true_BT_pre[0:BT_pre_index,:],axis=0)*1000 - 1.96*np.nanstd(meanRT_correct_true_BT_pre[0:BT_pre_index,:],axis=0)*1000
# axs[1,2].fill_between(evidenceright, CI95upper_BT_pre, CI95lower_BT_pre, color=precolor, alpha=0.25)

# CI95upper_BT_post = np.nanmean(meanRT_correct_true_BT_post[0:BT_post_index,:],axis=0)*1000 + 1.96*np.nanstd(meanRT_correct_true_BT_post[0:BT_post_index,:],axis=0)*1000
# CI95lower_BT_post = np.nanmean(meanRT_correct_true_BT_post[0:BT_post_index,:],axis=0)*1000 - 1.96*np.nanstd(meanRT_correct_true_BT_post[0:BT_post_index,:],axis=0)*1000
# axs[1,2].fill_between(evidenceright, CI95upper_BT_post, CI95lower_BT_post, color=postcolor, alpha=0.25)


# axs[0,1].set_xlabel('Coherence (%)', fontsize=fontsize)
axs[0,1].set_xticks(psychometric_xtick)
axs[0,1].set_xticklabels(psychometric_xtick)
axs[0,1].set_xlim((-40,40))
axs[0,1].set_ylim((0, 100))
axs[0,1].set_yticks(np.array([20, 40, 60, 80, 100]))
axs[0,1].set_ylabel('Choices to IF (%)', fontsize=fontsize)
# axs[0,1].get_yaxis().set_label_coords(-0.3,0.5)
custom_lines = [Line2D([0], [0], color=precolor, marker='o', markersize=markersize, linestyle='None'),
        Line2D([0], [0], color=postcolor, marker='o', markersize=markersize, linestyle='None'),
        Line2D([0], [0], color=precolorlite, lw=.5)]
axs[0,1].legend(custom_lines, ['Pre-','Post-', 'HDDM'],  loc=4, fontsize=legendsize)
# axs[0,1].set_ylabel('Choices to IF (%)', fontsize=fontsize)
axs[0,1].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[0,1].spines['right'].set_visible(False)
axs[0,1].spines['top'].set_visible(False)
axs[0,1].get_xaxis().set_visible(False)
axs[0,1].yaxis.set_ticks_position('left')
axs[0,1].xaxis.set_ticks_position('bottom')
axs[0,1].text(0.02, .98, 'r', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[0,1].transAxes)

axs[1,1].set_xlabel('Coherence (%)', fontsize=fontsize)
axs[1,1].set_xticks(psychometric_xtick)
axs[1,1].set_xticklabels(psychometric_xtick)
axs[1,1].set_xlim((-40,40))
axs[1,1].set_ylim((0, 100))
axs[1,1].set_yticks(np.array([20, 40, 60, 80, 100]))
# axs[1,1].set_ylabel('Choices to IF (%)', fontsize=fontsize)
axs[1,1].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[1,1].spines['right'].set_visible(False)
axs[1,1].spines['top'].set_visible(False)
axs[1,1].yaxis.set_ticks_position('left')
axs[1,1].xaxis.set_ticks_position('bottom')
axs[1,1].text(0.02, .98, 'v', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[1,1].transAxes)

# # axs[0,2].set_xlabel('Coherence (%)', fontsize=fontsize)
# axs[0,2].set_xticks(psychometric_xtick)
# axs[0,2].set_xticklabels(psychometric_xtick)
# axs[0,2].set_xlim((-40,40))
# axs[0,2].set_ylim((500,1200))
# axs[0,2].set_yticks(np.array([700, 1000, 1300]))
# axs[0,2].plot(np.array([0, 0]), np.array([500,1400]), color='k', LineStyle='--')
# # axs[0,2].set_ylabel('Mean Reaction Time (ms)', fontsize=fontsize)
# custom_lines = [Line2D([0], [0], color=precolor, marker='o', markersize=markersize, linestyle='None'),
#                 Line2D([0], [0], color=precolor, marker='s', markersize=markersize, linestyle='None'),
#                 Line2D([0], [0], color=precolorlite, lw=.5, linestyle='-'),
#                 Line2D([0], [0], color=precolorlite, lw=.5, linestyle='--')]
# axs[0,2].legend(custom_lines, ['Correct', 'Error', 'Fitted HDDM Correct', 'Fitted HDDM Error'],  loc=3, fontsize=legendsize, ncol=2)
# axs[0,2].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
# axs[0,2].spines['right'].set_visible(False)
# axs[0,2].spines['top'].set_visible(False)
# axs[0,2].get_xaxis().set_visible(False)
# axs[0,2].yaxis.set_ticks_position('left')
# axs[0,2].xaxis.set_ticks_position('bottom')
# axs[0,2].set_ylabel('Mean Reaction Time (ms)', fontsize=fontsize)
# axs[0,2].text(0.02, .98, 'c', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[0,2].transAxes)




# axs[1,2].set_xlabel('Coherence (%)', fontsize=fontsize)
# axs[1,2].set_xticks(psychometric_xtick)
# axs[1,2].set_xticklabels(psychometric_xtick)
# axs[1,2].set_xlim((-40,40))
# axs[1,2].set_ylim((700,1400))
# axs[1,2].set_yticks(np.array([900 , 1200, 1500]))
# axs[1,2].plot(np.array([0, 0]), np.array([700,1600]), color='k', LineStyle='--')
# axs[1,2].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
# axs[1,2].spines['right'].set_visible(False)
# axs[1,2].spines['top'].set_visible(False)
# axs[1,2].yaxis.set_ticks_position('left')
# axs[1,2].xaxis.set_ticks_position('bottom')
# axs[1,2].text(0.02, .98, 'f', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[1,2].transAxes)

# axs[0,2].set_xlabel('Coherence (%)', fontsize=fontsize)
axs[0,2].set_xticks(psychometric_xtick)
axs[0,2].set_xticklabels(psychometric_xtick)
axs[0,2].set_xlim((-40,40))
axs[0,2].set_ylim((500,1200))
axs[0,2].set_yticks(np.array([700, 1000, 1300]))
axs[0,2].plot(np.array([0, 0]), np.array([500,1400]), color='k', LineStyle='--')
# axs[0,2].set_ylabel('Mean Reaction Time (ms)', fontsize=fontsize)
custom_lines = [Line2D([0], [0], color=precolor, marker='o', markersize=markersize, linestyle='None'),
                Line2D([0], [0], color=postcolor, marker='o', markersize=markersize, linestyle='None'),
                Line2D([0], [0], color=precolorlite, lw=.5, linestyle='-')]
axs[0,2].legend(custom_lines, ['Pre- Correct', 'Post- Correct', 'HDDM'],  loc=3, fontsize=legendsize)
axs[0,2].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[0,2].spines['right'].set_visible(False)
axs[0,2].spines['top'].set_visible(False)
axs[0,2].get_xaxis().set_visible(False)
axs[0,2].yaxis.set_ticks_position('left')
axs[0,2].xaxis.set_ticks_position('bottom')
axs[0,2].set_ylabel('Mean Correct Reaction Time (ms)', fontsize=fontsize)
# axs[0,2].get_yaxis().set_label_coords(-0.3,0.5)
axs[0,2].text(0.02, .98, 's', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[0,2].transAxes)

axs[1,2].set_xlabel('Coherence (%)', fontsize=fontsize)
axs[1,2].set_xticks(psychometric_xtick)
axs[1,2].set_xticklabels(psychometric_xtick)
axs[1,2].set_xlim((-40,40))
axs[1,2].set_ylim((700,1400))
axs[1,2].set_yticks(np.array([900 , 1200, 1500]))
axs[1,2].plot(np.array([0, 0]), np.array([700,1600]), color='k', LineStyle='--')
axs[1,2].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[1,2].spines['right'].set_visible(False)
axs[1,2].spines['top'].set_visible(False)
axs[1,2].yaxis.set_ticks_position('left')
axs[1,2].xaxis.set_ticks_position('bottom')
axs[1,2].text(0.02, .98, 'w', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[1,2].transAxes)

axs[0,3].set_xticks(psychometric_xtick)
axs[0,3].set_xticklabels(psychometric_xtick)
axs[0,3].set_xlim((-40,40))
axs[0,3].set_ylim((500,1200))
axs[0,3].set_yticks(np.array([700, 1000, 1300]))
axs[0,3].plot(np.array([0, 0]), np.array([500,1400]), color='k', LineStyle='--')
# axs[0,3].set_ylabel('Mean Reaction Time (ms)', fontsize=fontsize)
custom_lines = [Line2D([0], [0], color=precolor, marker='s', markersize=markersize, linestyle='None'),
                Line2D([0], [0], color=postcolor, marker='s', markersize=markersize, linestyle='None'),
                Line2D([0], [0], color=precolorlite, lw=.5, linestyle='--')]
axs[0,3].legend(custom_lines, ['Pre- Error', 'Post- Error', 'HDDM'],  loc=3, fontsize=legendsize)
axs[0,3].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[0,3].spines['right'].set_visible(False)
axs[0,3].spines['top'].set_visible(False)
axs[0,3].get_xaxis().set_visible(False)
axs[0,3].yaxis.set_ticks_position('left')
axs[0,3].xaxis.set_ticks_position('bottom')
axs[0,3].set_ylabel('Mean Error Reaction Time (ms)', fontsize=fontsize)
# axs[0,3].get_yaxis().set_label_coords(-0.3,0.5)
axs[0,3].text(0.02, .98, 't', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[0,3].transAxes)

axs[1,3].set_xlabel('Coherence (%)', fontsize=fontsize)
axs[1,3].set_xticks(psychometric_xtick)
axs[1,3].set_xticklabels(psychometric_xtick)
axs[1,3].set_xlim((-40,40))
axs[1,3].set_ylim((700,1400))
axs[1,3].set_yticks(np.array([900 , 1200, 1500]))
axs[1,3].plot(np.array([0, 0]), np.array([700,1600]), color='k', LineStyle='--')
axs[1,3].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[1,3].spines['right'].set_visible(False)
axs[1,3].spines['top'].set_visible(False)
axs[1,3].yaxis.set_ticks_position('left')
axs[1,3].xaxis.set_ticks_position('bottom')
axs[1,3].text(0.02, .98, 'x', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[1,3].transAxes)



#For all subplots
# plt.subplots_adjust(hspace=0,left=0, right=0.1, bottom=0, top=0.1)
plt.subplots_adjust(hspace=0, wspace=0.4)
# plt.tight_layout(h_pad=-1,w_pad=0.3)
fig.set_size_inches(mm2inch(218,87),forward=False)

plt.savefig((f'../../figures/diffusion_figures/Figure4.png'), dpi=300, format='png',bbox_inches='tight')
plt.savefig((f'../../figures/diffusion_figures/Figure4.pdf'), dpi=300, format='pdf',bbox_inches='tight')
plt.savefig((f'../../figures/diffusion_figures/Figure4.svg'), dpi=300, format='svg',bbox_inches='tight')