# mky3_Figure5sims.py - Simulates a simple drift-diffusion model process, and generates modeling results plots in subfigures
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
# 03/19/20      Michael Nunez               Converted from mky3_simdiffusion.py
# 03/23/20      Michael Nunez              Place all figures in subplots
# 03/24/20      Michael Nunez            Create graphical representation for monkey SP
# 03/25/20      Michael Nunez             Adjust subplots to match Nature article specifications
# 03/30/20      Michael Nunez              Plot mean of mean RT curves across sessions
# 03/31/20      Michael Nunez              Add confidence interval to RT mean plots of real data
# 04/06/20      Michael Nunez              Remove top and right spines, remove CIs
# 04/13/20      Michael Nunez                  Use parameter estimates from simplified model
# 04/15/20      Michael Nunez                Plot all data (training and test) in the last two rows, plotting layout fixed
# 04/16/20      Michael Nunez             Bug fixes and plotting adjustments use tight_layout to enforce true figure size
# 04/20/20      Michael Nunez                  Shrink legend font size
# 04/21/20      Michael Nunez             Change Panel c legend location, save pdf and svg
# 04/24/20      Michael Nunez            Extract lapse proportion from model fitting data
# 04/28/20      Michael Nunez   Increase final mm2inch in set_size_inches so that the ultimate image is 183mm * 183mm after white space removal (checked with GIMP)
# 05/28/20      Michael Nunez                 Reduce figure to 3 rows, change transparency of RT curves
# 06/24/20      Michael Nunez                         pre_alpha == post_alpha
# 07/20/20      Michael Nunez          Change colors and drift-rate vector line thickness, 4 panels
# 07/21/20      Michael Nunez                Reordered plots

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


#alpha - amount of evidence units necessary to make a decision (evidence units)
#zeta - initial evidence bias parameter as an absolute value of evidence
#drift - average rate of evidence accumulation within a trial (evidence units per second)
#dcoef - diffusion coefficient, standard deviation of each sample of the evidence accumulation process  (evidence units per second)
#eta - standard deviation of the drift rate distribution across trials
#dt - length of time per sample (seconds)
#n - number of obervations of RT and choice
#maxsample - maximum number of samples allowed in the simulation


def fixedstartDDM(zeta,drift,aU,aL,dcoef,dt,n,maxsample,eta=0) :
    resp = []
    rt = []
    data = []
    allevidence = np.empty((n,maxsample))
    for i in range(n):
        evidence = []
        x=zeta
        sample=0
        thisdrift = drift + eta*np.random.normal()
        while (sample<maxsample):
            sample=sample+1
            evidence.append(x)
            x = x + dt*thisdrift + dcoef*np.sqrt(dt)*np.random.normal()
            if (x>=aU):
                resp.append(float(1.0))
                break
            if (x<=aL):
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


#Number of conditions
nconditions = 13
nmonkeys = 2

conditionspan_offset = np.linspace(-3,3,num=nconditions)
multiplier = .5

showdrifts_offset = np.array([conditionspan_offset[0], conditionspan_offset[3], conditionspan_offset[6], 
    conditionspan_offset[-4], conditionspan_offset[-1]])

#Simulated parameters
pre_alpha = 1.4
post_alpha = 1.4
pre_drift = 0
post_drift = -1
pre_ndt_ms = 0.45*1000
post_ndt_ms = 0.45*1000
pre_beta = 0.5
post_beta = pre_beta
alt_beta = 0.25

eta = 2 #The standard deviation across trials of the drift-rate

maxx = 2250 #Maximum number of ms to plot


dt = .001 #Seconds
lastsec = 2
time = np.arange(np.round(lastsec/float(dt)))*dt*1000


fig, axs = plt.subplots(4, 4,figsize=mm2inch(244,110), sharex='col')


#Plot boundaries
axs[3, 0].plot(np.array([0, maxx]), np.array([-pre_alpha/2, -pre_alpha/2]),linewidth=1,color=precolor,linestyle='-', zorder=13)
axs[3, 0].plot(np.array([0, maxx]), np.array([-post_alpha/2, -post_alpha/2]),linewidth=1,color=postcolor,linestyle='-', zorder=14)
axs[3, 0].plot(np.array([0, maxx]), np.array([pre_alpha/2, pre_alpha/2]),linewidth=1,color=precolor,linestyle='-', zorder=18)
axs[3, 0].plot(np.array([0, maxx]), np.array([post_alpha/2, post_alpha/2]),linewidth=1,color=postcolor,linestyle='-', zorder=19)

#Plot non-decision times
axs[3, 0].plot(np.array([0, pre_ndt_ms]), np.array([pre_beta-.5, pre_beta-.5]),linewidth=.5,color=precolorlite)
axs[3, 0].plot(np.array([0, post_ndt_ms]), np.array([post_beta-.5, post_beta-.5]),linewidth=.5,color=postcolorlite)


angles_radians_pre = np.empty((showdrifts_offset.shape[0]))
angles_radians_post = np.empty((showdrifts_offset.shape[0]))
angletrack = 0
for n in showdrifts_offset:
    if n == 0:
        scaleit = 5000
    else:
        scaleit = 8000
    arrowlength = 1000
    angle_pre = atan2( (n + pre_drift) , 1000 ) #Find the angle from start point line, drift rate is the average rate of accumulation per second
    angles_radians_pre[angletrack] = angle_pre
    xlen_pre = arrowlength*cos(angle_pre)
    ylen_pre = arrowlength*sin(angle_pre)
    angle_post = atan2( (n + post_drift) , 1000 ) #Find the angle from start point line, drift rate is the average rate of accumulation per second
    angles_radians_post[angletrack] = angle_post
    xlen_post = arrowlength*cos(angle_post)
    ylen_post = arrowlength*sin(angle_post)
    angletrack += 1
    if n == 0:
        axs[3, 0].quiver(pre_ndt_ms, pre_beta-.5, xlen_pre, ylen_pre, color=precolor, alpha=1, angles='xy', scale=scaleit, zorder=20)
        axs[3, 0].quiver(post_ndt_ms, post_beta-.5, xlen_post, ylen_post, color=postcolor, alpha=1, angles='xy', scale=scaleit, zorder=21)
    else:
        axs[3, 0].quiver(pre_ndt_ms, pre_beta-.5, xlen_pre, ylen_pre, color=precolor, alpha=1, angles='xy', scale=scaleit, width=0.005, zorder=20)
        axs[3, 0].quiver(post_ndt_ms, post_beta-.5, xlen_post, ylen_post, color=postcolor, alpha=1, angles='xy', scale=scaleit, width=0.005, zorder=21)

#Generate and plot evidence traces and RT data
(alldata_pre, allevidence_pre) = ratcliffDDM(beta=pre_beta, drift=pre_drift, alpha=pre_alpha, dcoef=1, eta=eta, dt=dt, n=500, maxsample=5000)
(alldata_post, allevidence_post) = ratcliffDDM(beta=post_beta, drift=post_drift, alpha=post_alpha, dcoef=1, eta=eta, dt=dt, n=500, maxsample=5000)

#Plot and shade RT distributions
kde_upper_pre = stats.gaussian_kde(alldata_pre[alldata_pre>0]*1000+pre_ndt_ms)
axs[3, 0].plot(time, kde_upper_pre(time)*400 + pre_alpha/2, linewidth=1,color=precolor,linestyle='-', zorder=9)
kde_lower_pre = stats.gaussian_kde(np.abs(alldata_pre[alldata_pre<0])*1000+pre_ndt_ms)
axs[3, 0].plot(time, kde_lower_pre(time)*-400 - pre_alpha/2, linewidth=1,color=precolor,linestyle='-', zorder=10)
axs[3, 0].fill_between(time, (pre_alpha/2)*np.ones((time.shape[0])), kde_upper_pre(time)*400 + pre_alpha/2, color=precolorlite, alpha=.5, zorder=5)
axs[3, 0].fill_between(time, -(pre_alpha/2)*np.ones((time.shape[0])), -kde_lower_pre(time)*400 - pre_alpha/2, color=precolorlite, alpha=.5, zorder=6)


kde_upper_post = stats.gaussian_kde(alldata_post[alldata_post>0]*1000+post_ndt_ms)
axs[3, 0].plot(time, kde_upper_post(time)*400 + post_alpha/2, linewidth=1,color=postcolor,linestyle='-', zorder=11)
kde_lower_post = stats.gaussian_kde(np.abs(alldata_post[alldata_post<0])*1000+post_ndt_ms)
axs[3, 0].plot(time, -kde_lower_post(time)*400 - post_alpha/2, linewidth=1,color=postcolor,linestyle='-', zorder=12)
axs[3, 0].fill_between(time, (post_alpha/2)*np.ones((time.shape[0])), kde_upper_post(time)*400 + post_alpha/2, color=postcolorlite, alpha=.5, zorder=7)
axs[3, 0].fill_between(time, -(post_alpha/2)*np.ones((time.shape[0])), -kde_lower_post(time)*400 - post_alpha/2, color=postcolorlite, alpha=.5, zorder=8)

# #Plot the evidence traces
allevidence_pre[allevidence_pre > pre_alpha] = pre_alpha #Correction for plotting first time pass distribution
allevidence_post[allevidence_post > post_alpha] = post_alpha #Correction for plotting first time pass distribution
allevidence_pre[allevidence_pre < 0] = 0 #Correction for plotting first time pass distribution
allevidence_post[allevidence_post < 0] = 0 #Correction for plotting first time pass distribution
axs[3, 0].plot(time + pre_ndt_ms, allevidence_pre[0:20,0:(time.shape[0])].T -pre_alpha/2, linewidth=.5,color=precolorlite, alpha=.25, zorder=1)
axs[3, 0].plot(time + post_ndt_ms, allevidence_post[0:20,0:(time.shape[0])].T -post_alpha/2, linewidth=.5,color=postcolorlite, alpha=.25, zorder=2)


#Appropriately label plots
# axs[3,0].set_xlabel('Time after Glass Pattern onset (ms)', fontsize=fontsize)
# axs[3,0].set_ylabel('Relative evidence for a choice to the IF', fontsize=fontsize)
axs[3,0].set_xlim([250,maxx]) #X limits in ms
axs[3,0].set_yticks(np.array([-1.0, 0, 1.0]))
# plt.tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[3,0].text(0.02, .98, 'm', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[3,0].transAxes)
axs[3,0].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[3,0].spines['right'].set_visible(False)
axs[3,0].spines['top'].set_visible(False)
axs[3,0].get_xaxis().set_visible(False)
axs[3,0].yaxis.set_ticks_position('left')
axs[3,0].xaxis.set_ticks_position('bottom')
# plt.savefig((f'../../figures/diffusion_figures/Diffusion_DriftOffset_example.png'), dpi=300, format='png',bbox_inches="tight")


#Plot the simulated psychometric curves
evidenceright = np.array([-36, -24, -17, -10, -5, -3, 0, 3, 5, 10, 17, 24, 36])
psychometric_xtick = np.array([-36, -17, -5, 0, 5, 17, 36])

psychometric_pre = np.ones((nconditions))*np.nan
psychometric_post = np.ones((nconditions))*np.nan
meanRT_correct_pre = np.ones((nconditions))*np.nan
meanRT_correct_post = np.ones((nconditions))*np.nan
meanRT_incorrect_pre = np.ones((nconditions))*np.nan
meanRT_incorrect_post = np.ones((nconditions))*np.nan
psychoindex = 0 
for n in conditionspan_offset:
    (alldata_pre, allevidence_pre) = ratcliffDDM(beta=pre_beta, drift=(n + pre_drift), alpha=pre_alpha, dcoef=1, eta=eta, dt=dt, n=5000, maxsample=5000)
    psychometric_pre[psychoindex] = np.sum(alldata_pre > 0)/alldata_pre.shape[0]
    if (psychoindex < np.round(nconditions/2)):
        meanRT_correct_pre[psychoindex] = np.mean(np.abs(alldata_pre[alldata_pre < 0])) + pre_ndt_ms/1000
        meanRT_incorrect_pre[psychoindex] = np.mean(alldata_pre[alldata_pre > 0]) + pre_ndt_ms/1000
    elif (psychoindex > np.round(nconditions/2)):
        meanRT_correct_pre[psychoindex] = np.mean(alldata_pre[alldata_pre > 0]) + pre_ndt_ms/1000
        meanRT_incorrect_pre[psychoindex] = np.mean(np.abs(alldata_pre[alldata_pre < 0])) + pre_ndt_ms/1000
    else:
        meanRT_correct_pre[psychoindex] = np.mean(np.abs(alldata_pre)) + pre_ndt_ms/1000
        meanRT_incorrect_pre[psychoindex] = np.mean(np.abs(alldata_pre)) + pre_ndt_ms/1000        
    (alldata_post, allevidence_post) = ratcliffDDM(beta=post_beta, drift=(n + post_drift), alpha=post_alpha, dcoef=1, eta=eta, dt=dt, n=5000, maxsample=5000)
    psychometric_post[psychoindex] = np.sum(alldata_post > 0)/alldata_post.shape[0]
    if (psychoindex < np.round(nconditions/2)):
        meanRT_correct_post[psychoindex] = np.mean(np.abs(alldata_post[alldata_post < 0])) + post_ndt_ms/1000
        meanRT_incorrect_post[psychoindex] = np.mean(alldata_post[alldata_post > 0]) + post_ndt_ms/1000
    elif (psychoindex > np.round(nconditions/2)):
        meanRT_correct_post[psychoindex] = np.mean(alldata_post[alldata_post > 0]) + post_ndt_ms/1000
        meanRT_incorrect_post[psychoindex] = np.mean(np.abs(alldata_post[alldata_post < 0])) + post_ndt_ms/1000
    else:
        meanRT_correct_post[psychoindex] = np.mean(np.abs(alldata_post)) + post_ndt_ms/1000
        meanRT_incorrect_post[psychoindex] = np.mean(np.abs(alldata_post)) + post_ndt_ms/1000

    psychoindex += 1




axs[3,1].plot(np.array([-50, 50]), np.array([50, 50]), color='k', LineStyle='--')
axs[3,1].plot(np.array([0, 0]), np.array([0, 100]), color='k', LineStyle='--')

popt_pre, pcov_pre = curve_fit(fpsychometric, evidenceright, psychometric_pre)
axs[3,1].plot(evidenceright, psychometric_pre*100, 'o',color=precolor, markersize=markersize)
axs[3,1].plot(evidenceright, fpsychometric(evidenceright, *popt_pre)*100, color=precolor, linewidth = 1)
popt_post, pcov_post = curve_fit(fpsychometric, evidenceright, psychometric_post)
axs[3,1].plot(evidenceright, psychometric_post*100, 'o',color=postcolor, markersize=markersize)
axs[3,1].plot(evidenceright, fpsychometric(evidenceright, *popt_post)*100, color=postcolor, linewidth = 1)

# axs[3,1].set_xlabel('Coherence (%)', fontsize=fontsize)
axs[3,1].set_xticks(psychometric_xtick)
axs[3,1].set_xticklabels(psychometric_xtick)
axs[3,1].set_ylim((0, 100))
axs[3,1].set_xlim((-40,40))
axs[3,1].set_yticks(np.array([20, 40, 60, 80, 100]))
# axs[3,1].set_ylabel('Choices to IF (%)', fontsize=fontsize)
axs[3,1].text(0.02, .98, 'n', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[3,1].transAxes)
axs[3,1].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[3,1].spines['right'].set_visible(False)
axs[3,1].spines['top'].set_visible(False)
axs[3,1].get_xaxis().set_visible(False)
axs[3,1].yaxis.set_ticks_position('left')
axs[3,1].xaxis.set_ticks_position('bottom')
# plt.savefig((f'../../figures/diffusion_figures/PSF_DriftOffset_example.png'), dpi=300, format='png',bbox_inches="tight")


axs[3,2].plot(evidenceright,meanRT_correct_pre*1000, 'o', color=precolor, markersize=markersize)
axs[3,2].plot(evidenceright,meanRT_correct_pre*1000, linewidth=1,color=precolor, zorder=2)
axs[3,3].plot(evidenceright,meanRT_incorrect_pre*1000, 's', color=precolor, markersize=markersize)
axs[3,3].plot(evidenceright,meanRT_incorrect_pre*1000, color=precolor, linewidth=1, linestyle='--', zorder=1)
axs[3,2].plot(evidenceright,meanRT_correct_post*1000, 'o', color=postcolor, markersize=markersize)
axs[3,2].plot(evidenceright,meanRT_correct_post*1000, linewidth=1,color=postcolor, zorder= 4)
axs[3,3].plot(evidenceright,meanRT_incorrect_post*1000, 's', color=postcolor, markersize=markersize)
axs[3,3].plot(evidenceright,meanRT_incorrect_post*1000, color=postcolor, linestyle='--', linewidth=1, zorder=3)

# axs[3,2].set_xlabel('Coherence (%)', fontsize=fontsize)
# axs[3,3].set_xlabel('Coherence (%)', fontsize=fontsize)

axs[3,2].set_xticks(psychometric_xtick)
axs[3,2].set_xticklabels(psychometric_xtick)
axs[3,2].set_ylim((500,1200))
axs[3,2].set_yticks(np.array([700, 900, 1100]))
# axes = plt.gca()
# ymin, ymax = axes.get_ylim()
axs[3,2].plot(np.array([0, 0]), np.array([500,1200]), color='k', LineStyle='--')
# axs[3,2].set_ylabel('Mean Reaction Time (ms)', fontsize=fontsize)
axs[3,2].text(0.02, .98, 'o', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[3,2].transAxes)
axs[3,2].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[3,2].spines['right'].set_visible(False)
axs[3,2].spines['top'].set_visible(False)
axs[3,2].get_xaxis().set_visible(False)
axs[3,2].yaxis.set_ticks_position('left')
axs[3,2].xaxis.set_ticks_position('bottom')

axs[3,3].set_xticks(psychometric_xtick)
axs[3,3].set_xticklabels(psychometric_xtick)
axs[3,3].set_ylim((500,1200))
axs[3,3].set_yticks(np.array([700, 900, 1100]))
# axes = plt.gca()
# ymin, ymax = axes.get_ylim()
axs[3,3].plot(np.array([0, 0]), np.array([500,1200]), color='k', LineStyle='--')
# axs[3,3].set_ylabel('Mean Reaction Time (ms)', fontsize=fontsize)
axs[3,3].text(0.02, .98, 'p', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[3,3].transAxes)
axs[3,3].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[3,3].spines['right'].set_visible(False)
axs[3,3].spines['top'].set_visible(False)
axs[3,3].get_xaxis().set_visible(False)
axs[3,3].yaxis.set_ticks_position('left')
axs[3,3].xaxis.set_ticks_position('bottom')

#### Show the effect of a multiplier on drift rate ####


#Plot boundaries
axs[0,0].plot(np.array([0, maxx]), np.array([-pre_alpha/2, -pre_alpha/2]),linewidth=1,color=precolor,linestyle='-', zorder=13)
axs[0,0].plot(np.array([0, maxx]), np.array([-post_alpha/2, -post_alpha/2]),linewidth=1,color=postcolor,linestyle='-', zorder=14)
axs[0,0].plot(np.array([0, maxx]), np.array([pre_alpha/2, pre_alpha/2]),linewidth=1,color=precolor,linestyle='-', zorder=18)
axs[0,0].plot(np.array([0, maxx]), np.array([post_alpha/2, post_alpha/2]),linewidth=1,color=postcolor,linestyle='-', zorder=19)

#Plot non-decision times
axs[0,0].plot(np.array([0, pre_ndt_ms]), np.array([pre_beta-.5, pre_beta-.5]),linewidth=.5,color=precolorlite)
axs[0,0].plot(np.array([0, post_ndt_ms]), np.array([post_beta-.5, post_beta-.5]),linewidth=.5,color=postcolorlite)


angles_radians_pre = np.empty((showdrifts_offset.shape[0]))
angles_radians_post = np.empty((showdrifts_offset.shape[0]))
angletrack = 0
for n in showdrifts_offset:
    if n == 0:
        scaleit = 5000
    else:
        scaleit = 8000
    arrowlength = 1000
    angle_pre = atan2( (n + pre_drift) , 1000 ) #Find the angle from start point line, drift rate is the average rate of accumulation per second
    angles_radians_pre[angletrack] = angle_pre
    xlen_pre = arrowlength*cos(angle_pre)
    ylen_pre = arrowlength*sin(angle_pre)
    angle_post = atan2( multiplier*(n + pre_drift) , 1000 ) #Find the angle from start point line, drift rate is the average rate of accumulation per second
    angles_radians_post[angletrack] = angle_post
    xlen_post = arrowlength*cos(angle_post)
    ylen_post = arrowlength*sin(angle_post)
    angletrack += 1
    if n == 0:
        axs[0,0].quiver(pre_ndt_ms, pre_beta-.5, xlen_pre, ylen_pre, color=precolor, alpha=1, angles='xy', scale=scaleit, zorder=20)
        axs[0,0].quiver(post_ndt_ms, post_beta-.5, xlen_post, ylen_post, color=postcolor, alpha=1, angles='xy', scale=scaleit, zorder=21)
    else:
        axs[0,0].quiver(pre_ndt_ms, pre_beta-.5, xlen_pre, ylen_pre, color=precolor, alpha=1, angles='xy', scale=scaleit, width=0.005, zorder=20)
        axs[0,0].quiver(post_ndt_ms, post_beta-.5, xlen_post, ylen_post, color=postcolor, alpha=1, angles='xy', scale=scaleit, width=0.005, zorder=21)

#Generate and plot evidence traces and RT data
(alldata_pre, allevidence_pre) = ratcliffDDM(beta=pre_beta, drift=pre_drift, alpha=pre_alpha, dcoef=1, eta=eta, dt=dt, n=500, maxsample=5000)
(alldata_post, allevidence_post) = ratcliffDDM(beta=post_beta, drift=pre_drift, alpha=post_alpha, dcoef=1, eta=eta, dt=dt, n=500, maxsample=5000)

#Plot and shade RT distributions
kde_upper_pre = stats.gaussian_kde(alldata_pre[alldata_pre>0]*1000+pre_ndt_ms)
axs[0,0].plot(time, kde_upper_pre(time)*400 + pre_alpha/2, linewidth=1,color=precolor,linestyle='-', zorder=9)
kde_lower_pre = stats.gaussian_kde(np.abs(alldata_pre[alldata_pre<0])*1000+pre_ndt_ms)
axs[0,0].plot(time, kde_lower_pre(time)*-400 - pre_alpha/2, linewidth=1,color=precolor,linestyle='-', zorder=10)
axs[0,0].fill_between(time, (pre_alpha/2)*np.ones((time.shape[0])), kde_upper_pre(time)*400 + pre_alpha/2, color=precolorlite, alpha=.5, zorder=5)
axs[0,0].fill_between(time, -(pre_alpha/2)*np.ones((time.shape[0])), -kde_lower_pre(time)*400 - pre_alpha/2, color=precolorlite, alpha=.5, zorder=6)


kde_upper_post = stats.gaussian_kde(alldata_post[alldata_post>0]*1000+post_ndt_ms)
axs[0,0].plot(time, kde_upper_post(time)*400 + post_alpha/2, linewidth=1,color=postcolor,linestyle='-', zorder=11)
kde_lower_post = stats.gaussian_kde(np.abs(alldata_post[alldata_post<0])*1000+post_ndt_ms)
axs[0,0].plot(time, -kde_lower_post(time)*400 - post_alpha/2, linewidth=1,color=postcolor,linestyle='-', zorder=12)
axs[0,0].fill_between(time, (post_alpha/2)*np.ones((time.shape[0])), kde_upper_post(time)*400 + post_alpha/2, color=postcolorlite, alpha=.5, zorder=7)
axs[0,0].fill_between(time, -(post_alpha/2)*np.ones((time.shape[0])), -kde_lower_post(time)*400 - post_alpha/2, color=postcolorlite, alpha=.5, zorder=8)

# #Plot the evidence traces
allevidence_pre[allevidence_pre > pre_alpha] = pre_alpha #Correction for plotting first time pass distribution
allevidence_post[allevidence_post > post_alpha] = post_alpha #Correction for plotting first time pass distribution
allevidence_pre[allevidence_pre < 0] = 0 #Correction for plotting first time pass distribution
allevidence_post[allevidence_post < 0] = 0 #Correction for plotting first time pass distribution
axs[0,0].plot(time + pre_ndt_ms, allevidence_pre[0:20,0:(time.shape[0])].T -pre_alpha/2, linewidth=.5,color=precolorlite, alpha=.25, zorder=1)
axs[0,0].plot(time + post_ndt_ms, allevidence_post[0:20,0:(time.shape[0])].T -post_alpha/2, linewidth=.5,color=postcolorlite, alpha=.25, zorder=2)


#Appropriately label plots
# axs[0,0].set_xlabel('Time after Glass Pattern onset (ms)', fontsize=fontsize)
# axs[0,0].set_ylabel('Relative evidence for a choice to the IF', fontsize=fontsize)
axs[0,0].set_xlim([250,maxx]) #X limits in ms
axs[0,0].set_yticks(np.array([-1.0, 0, 1.0]))
axs[0,0].text(0.02, .98, 'a', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[0,0].transAxes)
axs[0,0].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[0,0].spines['right'].set_visible(False)
axs[0,0].spines['top'].set_visible(False)
axs[0,0].get_xaxis().set_visible(False)
axs[0,0].yaxis.set_ticks_position('left')
axs[0,0].xaxis.set_ticks_position('bottom')
custom_lines = [Line2D([0], [0], color=precolor, linestyle='-', lw=1),
                Line2D([0], [0], color=postcolor, lw=1)]
axs[0,0].legend(custom_lines, ['Pre-','Post-'],  loc=4, fontsize=legendsize)
# plt.savefig((f'../../figures/diffusion_figures/Diffusion_DriftMulti_example.png'), dpi=300, format='png',bbox_inches="tight")


#Plot the simulated psychometric curves
evidenceright = np.array([-36, -24, -17, -10, -5, -3, 0, 3, 5, 10, 17, 24, 36])

psychometric_pre = np.ones((nconditions))*np.nan
psychometric_post = np.ones((nconditions))*np.nan
meanRT_correct_pre = np.ones((nconditions))*np.nan
meanRT_correct_post = np.ones((nconditions))*np.nan
meanRT_incorrect_pre = np.ones((nconditions))*np.nan
meanRT_incorrect_post = np.ones((nconditions))*np.nan
psychoindex = 0 
for n in conditionspan_offset:
    (alldata_pre, allevidence_pre) = ratcliffDDM(beta=pre_beta, drift=(n + pre_drift), alpha=pre_alpha, dcoef=1, eta=eta, dt=dt, n=5000, maxsample=5000)
    psychometric_pre[psychoindex] = np.sum(alldata_pre > 0)/alldata_pre.shape[0]
    if (psychoindex < np.round(nconditions/2)):
        meanRT_correct_pre[psychoindex] = np.mean(np.abs(alldata_pre[alldata_pre < 0])) + pre_ndt_ms/1000
        meanRT_incorrect_pre[psychoindex] = np.mean(alldata_pre[alldata_pre > 0]) + pre_ndt_ms/1000
    elif (psychoindex > np.round(nconditions/2)):
        meanRT_correct_pre[psychoindex] = np.mean(alldata_pre[alldata_pre > 0]) + pre_ndt_ms/1000
        meanRT_incorrect_pre[psychoindex] = np.mean(np.abs(alldata_pre[alldata_pre < 0])) + pre_ndt_ms/1000
    else:
        meanRT_correct_pre[psychoindex] = np.mean(np.abs(alldata_pre)) + pre_ndt_ms/1000
        meanRT_incorrect_pre[psychoindex] = np.mean(np.abs(alldata_pre)) + pre_ndt_ms/1000        
    (alldata_post, allevidence_post) = ratcliffDDM(beta=post_beta, drift=multiplier*(n + pre_drift), alpha=post_alpha, dcoef=1, eta=eta, dt=dt, n=5000, maxsample=5000)
    psychometric_post[psychoindex] = np.sum(alldata_post > 0)/alldata_post.shape[0]
    if (psychoindex < np.round(nconditions/2)):
        meanRT_correct_post[psychoindex] = np.mean(np.abs(alldata_post[alldata_post < 0])) + post_ndt_ms/1000
        meanRT_incorrect_post[psychoindex] = np.mean(alldata_post[alldata_post > 0]) + post_ndt_ms/1000
    elif (psychoindex > np.round(nconditions/2)):
        meanRT_correct_post[psychoindex] = np.mean(alldata_post[alldata_post > 0]) + post_ndt_ms/1000
        meanRT_incorrect_post[psychoindex] = np.mean(np.abs(alldata_post[alldata_post < 0])) + post_ndt_ms/1000
    else:
        meanRT_correct_post[psychoindex] = np.mean(np.abs(alldata_post)) + post_ndt_ms/1000
        meanRT_incorrect_post[psychoindex] = np.mean(np.abs(alldata_post)) + post_ndt_ms/1000

    psychoindex += 1




axs[0,1].plot(np.array([-50, 50]), np.array([50, 50]), color='k', LineStyle='--')
axs[0,1].plot(np.array([0, 0]), np.array([0, 100]), color='k', LineStyle='--')

popt_pre, pcov_pre = curve_fit(fpsychometric, evidenceright, psychometric_pre)
axs[0,1].plot(evidenceright, psychometric_pre*100, 'o',color=precolor, markersize=markersize)
axs[0,1].plot(evidenceright, fpsychometric(evidenceright, *popt_pre)*100, color=precolor, linewidth = 1)
popt_post, pcov_post = curve_fit(fpsychometric, evidenceright, psychometric_post)
axs[0,1].plot(evidenceright, psychometric_post*100, 'o',color=postcolor, markersize=markersize)
axs[0,1].plot(evidenceright, fpsychometric(evidenceright, *popt_post)*100, color=postcolor, linewidth = 1)

# axs[0,1].set_xlabel('Coherence (%)', fontsize=fontsize)
axs[0,1].set_xticks(psychometric_xtick)
axs[0,1].set_xticklabels(psychometric_xtick)
axs[0,1].set_ylim((0, 100))
axs[0,1].set_xlim((-40,40))
axs[0,1].set_yticks(np.array([20, 40, 60, 80, 100]))
axs[0,1].text(0.02, .98, 'b', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[0,1].transAxes)
axs[0,1].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[0,1].spines['right'].set_visible(False)
axs[0,1].spines['top'].set_visible(False)
axs[0,1].get_xaxis().set_visible(False)
axs[0,1].yaxis.set_ticks_position('left')
axs[0,1].xaxis.set_ticks_position('bottom')
custom_lines = [Line2D([0], [0], color=precolor, linestyle='-', lw=1),
                Line2D([0], [0], color=postcolor, lw=1)]
axs[0,1].legend(custom_lines, ['Pre-','Post-'],  loc=4, fontsize=legendsize)

axs[0,2].plot(evidenceright,meanRT_correct_pre*1000, 'o', color=precolor, markersize=markersize)
axs[0,2].plot(evidenceright,meanRT_correct_pre*1000, linewidth=1,color=precolor, zorder=2)
axs[0,3].plot(evidenceright,meanRT_incorrect_pre*1000, 's', color=precolor, markersize=markersize)
axs[0,3].plot(evidenceright,meanRT_incorrect_pre*1000, color=precolor, linewidth=1, linestyle='--', zorder=1)
axs[0,2].plot(evidenceright,meanRT_correct_post*1000, 'o', color=postcolor, markersize=markersize)
axs[0,2].plot(evidenceright,meanRT_correct_post*1000, linewidth=1,color=postcolor, zorder= 4)
axs[0,3].plot(evidenceright,meanRT_incorrect_post*1000, 's', color=postcolor, markersize=markersize)
axs[0,3].plot(evidenceright,meanRT_incorrect_post*1000, color=postcolor, linestyle='--', linewidth=1, zorder=3)

axs[0,2].set_xticks(psychometric_xtick)
axs[0,2].set_xticklabels(psychometric_xtick)
axs[0,2].set_ylim((500,1200))
axs[0,2].set_yticks(np.array([700, 900, 1100]))
# axes = plt.gca()
# ymin, ymax = axes.get_ylim()
axs[0,2].plot(np.array([0, 0]), np.array([500,1200]), color='k', LineStyle='--')
# axs[0,2].set_ylabel('Mean Reaction Time (ms)', fontsize=fontsize)
custom_lines = [Line2D([0], [0], marker='o', color=precolor, lw=.5, markersize=markersize),
                Line2D([0], [0], marker='o', color=postcolor, lw=.5, markersize=markersize)]
axs[0,2].legend(custom_lines, ['Pre- Correct','Post- Correct'],  loc=3, fontsize=legendsize)
axs[0,2].text(0.02, .98, 'c', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[0,2].transAxes)
axs[0,2].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[0,2].spines['right'].set_visible(False)
axs[0,2].spines['top'].set_visible(False)
axs[0,2].get_xaxis().set_visible(False)
axs[0,2].yaxis.set_ticks_position('left')
axs[0,2].xaxis.set_ticks_position('bottom')

axs[0,3].set_xticks(psychometric_xtick)
axs[0,3].set_xticklabels(psychometric_xtick)
axs[0,3].set_ylim((500,1200))
axs[0,3].set_yticks(np.array([700, 900, 1100]))
# axes = plt.gca()
# ymin, ymax = axes.get_ylim()
axs[0,3].plot(np.array([0, 0]), np.array([500,1200]), color='k', LineStyle='--')
# axs[0,3].set_ylabel('Mean Reaction Time (ms)', fontsize=fontsize)
custom_lines = [Line2D([0], [0], marker='s', linestyle='--', color=precolor, lw=.5, markersize=markersize),
                Line2D([0], [0], marker='s', linestyle='--', color=postcolor, lw=.5, markersize=markersize)]
axs[0,3].legend(custom_lines, ['Pre- Error', 'Post- Error'],  loc=3, fontsize=legendsize)
axs[0,3].text(0.02, .98, 'd', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[0,3].transAxes)
axs[0,3].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[0,3].spines['right'].set_visible(False)
axs[0,3].spines['top'].set_visible(False)
axs[0,3].get_xaxis().set_visible(False)
axs[0,3].yaxis.set_ticks_position('left')
axs[0,3].xaxis.set_ticks_position('bottom')

#### Show the effect of a start point change ####

#Plot boundaries
axs[1,0].plot(np.array([0, maxx]), np.array([-pre_alpha/2, -pre_alpha/2]),linewidth=1,color=precolor,linestyle='-', zorder=13)
axs[1,0].plot(np.array([0, maxx]), np.array([-post_alpha/2, -post_alpha/2]),linewidth=1,color=postcolor,linestyle='-', zorder=14)
axs[1,0].plot(np.array([0, maxx]), np.array([pre_alpha/2, pre_alpha/2]),linewidth=1,color=precolor,linestyle='-', zorder=18)
axs[1,0].plot(np.array([0, maxx]), np.array([post_alpha/2, post_alpha/2]),linewidth=1,color=postcolor,linestyle='-', zorder=19)
# plt.plot(np.array([0, maxx]), np.array([0, 0]),linewidth=1,color='k',linestyle='-',zorder=15)

#Plot non-decision times
axs[1,0].plot(np.array([0, pre_ndt_ms]), np.array([pre_beta-.5, pre_beta-.5]),linewidth=.5,color=precolorlite)
axs[1,0].plot(np.array([0, post_ndt_ms]), np.array([alt_beta-.5, alt_beta-.5]),linewidth=.5,color=postcolorlite)


angles_radians_pre = np.empty((showdrifts_offset.shape[0]))
angles_radians_post = np.empty((showdrifts_offset.shape[0]))
angletrack = 0
for n in showdrifts_offset:
    if n == 0:
        scaleit = 5000
    else:
        scaleit = 8000
    arrowlength = 1000
    angle_pre = atan2( (n + pre_drift) , 1000 ) #Find the angle from start point line, drift rate is the average rate of accumulation per second
    angles_radians_pre[angletrack] = angle_pre
    xlen_pre = arrowlength*cos(angle_pre)
    ylen_pre = arrowlength*sin(angle_pre)
    angle_post = atan2( (n + pre_drift) , 1000 ) #Find the angle from start point line, drift rate is the average rate of accumulation per second
    angles_radians_post[angletrack] = angle_post
    xlen_post = arrowlength*cos(angle_post)
    ylen_post = arrowlength*sin(angle_post)
    angletrack += 1
    if n == 0:
        axs[1,0].quiver(pre_ndt_ms, pre_beta-.5, xlen_pre, ylen_pre, color=precolor, alpha=1, angles='xy', scale=scaleit, zorder=20)
        axs[1,0].quiver(post_ndt_ms, alt_beta-.5, xlen_post, ylen_post, color=postcolor, alpha=1, angles='xy', scale=scaleit, zorder=21)
    else:
        axs[1,0].quiver(pre_ndt_ms, pre_beta-.5, xlen_pre, ylen_pre, color=precolor, alpha=1, angles='xy', scale=scaleit, width=0.005, zorder=20)
        axs[1,0].quiver(post_ndt_ms, alt_beta-.5, xlen_post, ylen_post, color=postcolor, alpha=1, angles='xy', scale=scaleit, width=0.005, zorder=21)


#Generate and plot evidence traces and RT data
(alldata_pre, allevidence_pre) = ratcliffDDM(beta=pre_beta, drift=pre_drift, alpha=pre_alpha, dcoef=1, eta=eta, dt=dt, n=500, maxsample=5000)
(alldata_post, allevidence_post) = ratcliffDDM(beta=alt_beta, drift=pre_drift, alpha=post_alpha, dcoef=1, eta=eta, dt=dt, n=500, maxsample=5000)

#Plot and shade RT distributions
kde_upper_pre = stats.gaussian_kde(alldata_pre[alldata_pre>0]*1000+pre_ndt_ms)
axs[1,0].plot(time, kde_upper_pre(time)*400 + pre_alpha/2, linewidth=1,color=precolor,linestyle='-', zorder=9)
kde_lower_pre = stats.gaussian_kde(np.abs(alldata_pre[alldata_pre<0])*1000+pre_ndt_ms)
axs[1,0].plot(time, kde_lower_pre(time)*-400 - pre_alpha/2, linewidth=1,color=precolor,linestyle='-', zorder=10)
axs[1,0].fill_between(time, (pre_alpha/2)*np.ones((time.shape[0])), kde_upper_pre(time)*400 + pre_alpha/2, color=precolorlite, alpha=.5, zorder=5)
axs[1,0].fill_between(time, -(pre_alpha/2)*np.ones((time.shape[0])), -kde_lower_pre(time)*400 - pre_alpha/2, color=precolorlite, alpha=.5, zorder=6)


kde_upper_post = stats.gaussian_kde(alldata_post[alldata_post>0]*1000+post_ndt_ms)
axs[1,0].plot(time, kde_upper_post(time)*400 + post_alpha/2, linewidth=1,color=postcolor,linestyle='-', zorder=11)
kde_lower_post = stats.gaussian_kde(np.abs(alldata_post[alldata_post<0])*1000+post_ndt_ms)
axs[1,0].plot(time, -kde_lower_post(time)*400 - post_alpha/2, linewidth=1,color=postcolor,linestyle='-', zorder=12)
axs[1,0].fill_between(time, (post_alpha/2)*np.ones((time.shape[0])), kde_upper_post(time)*400 + post_alpha/2, color=postcolorlite, alpha=.5, zorder=7)
axs[1,0].fill_between(time, -(post_alpha/2)*np.ones((time.shape[0])), -kde_lower_post(time)*400 - post_alpha/2, color=postcolorlite, alpha=.5, zorder=8)

#Plot the evidence traces
allevidence_pre[allevidence_pre > pre_alpha] = pre_alpha #Correction for plotting first time pass distribution
allevidence_post[allevidence_post > post_alpha] = post_alpha #Correction for plotting first time pass distribution
allevidence_pre[allevidence_pre < 0] = 0 #Correction for plotting first time pass distribution
allevidence_post[allevidence_post < 0] = 0 #Correction for plotting first time pass distribution
axs[1,0].plot(time + pre_ndt_ms, allevidence_pre[0:20,0:(time.shape[0])].T -pre_alpha/2, linewidth=.5,color=precolorlite, alpha=.25, zorder=1)
axs[1,0].plot(time + post_ndt_ms, allevidence_post[0:20,0:(time.shape[0])].T -post_alpha/2, linewidth=.5,color=postcolorlite, alpha=.25, zorder=2)


#Appropriately label plots
# axs[1,0].set_xlabel('Time after Glass Pattern onset (ms)', fontsize=fontsize)
axs[1,0].set_ylabel('Relative evidence for a choice to the IF', fontsize=fontsize)
axs[1,0].set_xlim([250,maxx]) #X limits in ms
axs[1,0].set_yticks(np.array([-1.0, 0, 1.0]))
axs[1,0].text(0.02, .98, 'e', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[1,0].transAxes)
# custom_lines = [Line2D([0], [0], color=precolor, lw=1),
                # Line2D([0], [0], color=postcolor, lw=1)]
# axs[1,0].legend(custom_lines, ['Pre-','Post-'],  loc=2, fontsize=legendsize)
# plt.savefig((f'../../figures/diffusion_figures/Diffusion_SP_example.png'), dpi=300, format='png',bbox_inches="tight")
axs[1,0].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[1,0].spines['right'].set_visible(False)
axs[1,0].spines['top'].set_visible(False)
axs[1,0].get_xaxis().set_visible(False)
axs[1,0].yaxis.set_ticks_position('left')
axs[1,0].xaxis.set_ticks_position('bottom')


psychometric_pre = np.ones((nconditions))*np.nan
psychometric_post = np.ones((nconditions))*np.nan
meanRT_correct_pre = np.ones((nconditions))*np.nan
meanRT_correct_post = np.ones((nconditions))*np.nan
meanRT_incorrect_pre = np.ones((nconditions))*np.nan
meanRT_incorrect_post = np.ones((nconditions))*np.nan
psychoindex = 0 
for n in conditionspan_offset:
    (alldata_pre, allevidence_pre) = ratcliffDDM(beta=pre_beta, drift=(n + pre_drift), alpha=pre_alpha, dcoef=1, eta=eta, dt=dt, n=5000, maxsample=5000)
    psychometric_pre[psychoindex] = np.sum(alldata_pre > 0)/alldata_pre.shape[0]
    if (psychoindex < np.round(nconditions/2)):
        meanRT_correct_pre[psychoindex] = np.mean(np.abs(alldata_pre[alldata_pre < 0])) + pre_ndt_ms/1000
        meanRT_incorrect_pre[psychoindex] = np.mean(alldata_pre[alldata_pre > 0]) + pre_ndt_ms/1000
    elif (psychoindex > np.round(nconditions/2)):
        meanRT_correct_pre[psychoindex] = np.mean(alldata_pre[alldata_pre > 0]) + pre_ndt_ms/1000
        meanRT_incorrect_pre[psychoindex] = np.mean(np.abs(alldata_pre[alldata_pre < 0])) + pre_ndt_ms/1000
    else:
        meanRT_correct_pre[psychoindex] = np.mean(np.abs(alldata_pre)) + pre_ndt_ms/1000
        meanRT_incorrect_pre[psychoindex] = np.mean(np.abs(alldata_pre)) + pre_ndt_ms/1000        
    (alldata_post, allevidence_post) = ratcliffDDM(beta=alt_beta, drift=(n + pre_drift), alpha=post_alpha, dcoef=1, eta=eta, dt=dt, n=5000, maxsample=5000)
    psychometric_post[psychoindex] = np.sum(alldata_post > 0)/alldata_post.shape[0]
    if (psychoindex < np.round(nconditions/2)):
        meanRT_correct_post[psychoindex] = np.mean(np.abs(alldata_post[alldata_post < 0])) + post_ndt_ms/1000
        meanRT_incorrect_post[psychoindex] = np.mean(alldata_post[alldata_post > 0]) + post_ndt_ms/1000
    elif (psychoindex > np.round(nconditions/2)):
        meanRT_correct_post[psychoindex] = np.mean(alldata_post[alldata_post > 0]) + post_ndt_ms/1000
        meanRT_incorrect_post[psychoindex] = np.mean(np.abs(alldata_post[alldata_post < 0])) + post_ndt_ms/1000
    else:
        meanRT_correct_post[psychoindex] = np.mean(np.abs(alldata_post)) + post_ndt_ms/1000
        meanRT_incorrect_post[psychoindex] = np.mean(np.abs(alldata_post)) + post_ndt_ms/1000

    psychoindex += 1




axs[1,1].plot(np.array([-50, 50]), np.array([50, 50]), color='k', LineStyle='--')
axs[1,1].plot(np.array([0, 0]), np.array([0, 100]), color='k', LineStyle='--')

popt_pre, pcov_pre = curve_fit(fpsychometric, evidenceright, psychometric_pre)
axs[1,1].plot(evidenceright, psychometric_pre*100, 'o',color=precolor, markersize=markersize)
axs[1,1].plot(evidenceright, fpsychometric(evidenceright, *popt_pre)*100, color=precolor, linewidth = 1)
popt_post, pcov_post = curve_fit(fpsychometric, evidenceright, psychometric_post)
axs[1,1].plot(evidenceright, psychometric_post*100, 'o',color=postcolor, markersize=markersize)
axs[1,1].plot(evidenceright, fpsychometric(evidenceright, *popt_post)*100, color=postcolor, linewidth = 1)

# axs[1,1].set_xlabel('Coherence (%)', fontsize=fontsize)
axs[1,1].set_xticks(psychometric_xtick)
axs[1,1].set_xticklabels(psychometric_xtick)
axs[1,1].set_ylim((0, 100))
axs[1,1].set_xlim((-40,40))
axs[1,1].set_yticks(np.array([20, 40, 60, 80, 100]))
axs[1,1].set_ylabel('Choices to IF (%)', fontsize=fontsize)
axs[1,1].text(0.02, .98, 'f', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[1,1].transAxes)
axs[1,1].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[1,1].spines['right'].set_visible(False)
axs[1,1].spines['top'].set_visible(False)
axs[1,1].get_xaxis().set_visible(False)
axs[1,1].yaxis.set_ticks_position('left')
axs[1,1].xaxis.set_ticks_position('bottom')
# axs[1,1].legend(custom_lines, ['Pre-','Post-'],  loc=2, fontsize=legendsize)
# plt.savefig((f'../../figures/diffusion_figures/PSF_SP_example.png'), dpi=300, format='png',bbox_inches="tight")

axs[1,2].plot(evidenceright,meanRT_correct_pre*1000, 'o', color=precolor, markersize=markersize)
axs[1,2].plot(evidenceright,meanRT_correct_pre*1000, linewidth=1,color=precolor, zorder=2)
axs[1,3].plot(evidenceright,meanRT_incorrect_pre*1000, 's', color=precolor, markersize=markersize)
axs[1,3].plot(evidenceright,meanRT_incorrect_pre*1000, color=precolor, linewidth=1, linestyle='--', zorder=1)
axs[1,2].plot(evidenceright,meanRT_correct_post*1000, 'o', color=postcolor, markersize=markersize)
axs[1,2].plot(evidenceright,meanRT_correct_post*1000, linewidth=1,color=postcolor, zorder= 4)
axs[1,3].plot(evidenceright,meanRT_incorrect_post*1000, 's', color=postcolor, markersize=markersize)
axs[1,3].plot(evidenceright,meanRT_incorrect_post*1000, color=postcolor, linestyle='--', linewidth=1, zorder=3)

# axs[1,2].set_xlabel('Coherence (%)', fontsize=fontsize)
axs[1,2].set_xticks(psychometric_xtick)
axs[1,2].set_xticklabels(psychometric_xtick)
axs[1,2].set_ylim((500,1200))
axs[1,2].set_yticks(np.array([700, 900, 1100]))
# axes = plt.gca()
# ymin, ymax = axes.get_ylim()
axs[1,2].plot(np.array([0, 0]), np.array([500,1200]), color='k', LineStyle='--')
axs[1,2].set_ylabel('Mean Correct Reaction Time (ms)', fontsize=fontsize)
axs[1,2].text(0.02, .98, 'g', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[1,2].transAxes)
axs[1,2].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[1,2].spines['right'].set_visible(False)
axs[1,2].spines['top'].set_visible(False)
axs[1,2].get_xaxis().set_visible(False)
axs[1,2].yaxis.set_ticks_position('left')
axs[1,2].xaxis.set_ticks_position('bottom')

# axs[1,2].set_xlabel('Coherence (%)', fontsize=fontsize)
axs[1,3].set_xticks(psychometric_xtick)
axs[1,3].set_xticklabels(psychometric_xtick)
axs[1,3].set_ylim((500,1200))
axs[1,3].set_yticks(np.array([700, 900, 1100]))
# axes = plt.gca()
# ymin, ymax = axes.get_ylim()
axs[1,3].plot(np.array([0, 0]), np.array([500,1200]), color='k', LineStyle='--')
axs[1,3].set_ylabel('Mean Error Reaction Time (ms)', fontsize=fontsize)
axs[1,3].text(0.02, .98, 'h', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[1,3].transAxes)
axs[1,3].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[1,3].spines['right'].set_visible(False)
axs[1,3].spines['top'].set_visible(False)
axs[1,3].get_xaxis().set_visible(False)
axs[1,3].yaxis.set_ticks_position('left')
axs[1,3].xaxis.set_ticks_position('bottom')


### Simulate decrease in one boundary ###

#Simulated parameters
pre_aU = .7
post_aU = 1
pre_aL = -.7
post_aL = -.7
pre_drift = 0
post_drift = 0
pre_ndt_ms = 0.45*1000
post_ndt_ms = 0.45*1000
pre_zeta = 0
post_zeta = pre_zeta

#Plot boundaries
axs[2,0].plot(np.array([0, maxx]), np.array([pre_aL, pre_aL]),linewidth=1,color=precolor,linestyle='-', zorder=13)
axs[2,0].plot(np.array([0, maxx]), np.array([post_aL, post_aL]),linewidth=1,color=postcolor,linestyle='-', zorder=14)
axs[2,0].plot(np.array([0, maxx]), np.array([pre_aU, pre_aU]),linewidth=1,color=precolor,linestyle='-', zorder=18)
axs[2,0].plot(np.array([0, maxx]), np.array([post_aU, post_aU]),linewidth=1,color=postcolor,linestyle='-', zorder=19)


#Plot non-decision times
axs[2,0].plot(np.array([0, pre_ndt_ms]), np.array([pre_zeta, pre_zeta]),linewidth=.5,color=precolorlite)
axs[2,0].plot(np.array([0, pre_ndt_ms]), np.array([post_zeta, post_zeta]),linewidth=.5,color=postcolorlite)


angles_radians_pre = np.empty((showdrifts_offset.shape[0]))
angles_radians_post = np.empty((showdrifts_offset.shape[0]))
angletrack = 0
for n in showdrifts_offset:
    if n == 0:
        scaleit = 5000
    else:
        scaleit = 8000
    arrowlength = 1000
    angle_pre = atan2( (n + pre_drift) , 1000 ) #Find the angle from start point line, drift rate is the average rate of accumulation per second
    angles_radians_pre[angletrack] = angle_pre
    xlen_pre = arrowlength*cos(angle_pre)
    ylen_pre = arrowlength*sin(angle_pre)
    angle_post = atan2( (n + post_drift) , 1000 ) #Find the angle from start point line, drift rate is the average rate of accumulation per second
    angles_radians_post[angletrack] = angle_post
    xlen_post = arrowlength*cos(angle_post)
    ylen_post = arrowlength*sin(angle_post)
    angletrack += 1
    if n == 0:
        axs[2,0].quiver(pre_ndt_ms, pre_zeta, xlen_pre, ylen_pre, color=precolor, alpha=1, angles='xy', scale=scaleit, zorder=20)
        axs[2,0].quiver(pre_ndt_ms, post_zeta, xlen_post, ylen_post, color=postcolor, alpha=1, angles='xy', scale=scaleit, zorder=21)
    else:
        axs[2,0].quiver(pre_ndt_ms, pre_zeta, xlen_pre, ylen_pre, color=precolor, alpha=1, angles='xy', scale=scaleit, width=0.005, zorder=20)
        axs[2,0].quiver(pre_ndt_ms, post_zeta, xlen_post, ylen_post, color=postcolor, alpha=1, angles='xy', scale=scaleit, width=0.005, zorder=21)


#Generate and plot evidence traces and RT data
(alldata_pre, allevidence_pre) = fixedstartDDM(zeta=pre_zeta, drift=pre_drift, aU=pre_aU, aL=pre_aL, dcoef=1, eta=eta, dt=dt, n=500, maxsample=5000)
(alldata_post, allevidence_post) = fixedstartDDM(zeta=post_zeta, drift=post_drift, aU=post_aU, aL=post_aL, dcoef=1, eta=eta, dt=dt, n=500, maxsample=5000)


#Plot and shade RT distributions
kde_upper_pre = stats.gaussian_kde(alldata_pre[alldata_pre>0]*1000+pre_ndt_ms)
axs[2,0].plot(time, kde_upper_pre(time)*400 + pre_aU, linewidth=1,color=precolor,linestyle='-', zorder=9)
kde_lower_pre = stats.gaussian_kde(np.abs(alldata_pre[alldata_pre<0])*1000+pre_ndt_ms)
axs[2,0].plot(time, kde_lower_pre(time)*-400 + pre_aL, linewidth=1,color=precolor,linestyle='-', zorder=10)
axs[2,0].fill_between(time, pre_aU*np.ones((time.shape[0])), kde_upper_pre(time)*400 + pre_aU, color=precolorlite, alpha=.5, zorder=5)
axs[2,0].fill_between(time, pre_aL*np.ones((time.shape[0])), -kde_lower_pre(time)*400 + pre_aL, color=precolorlite, alpha=.5,  zorder=6)


kde_upper_post = stats.gaussian_kde(alldata_post[alldata_post>0]*1000+pre_ndt_ms)
axs[2,0].plot(time, kde_upper_post(time)*400 + post_aU, linewidth=1,color=postcolor,linestyle='-', zorder=11)
kde_lower_post = stats.gaussian_kde(np.abs(alldata_post[alldata_post<0])*1000+pre_ndt_ms)
axs[2,0].plot(time, -kde_lower_post(time)*400 + post_aL, linewidth=1,color=postcolor,linestyle='-', zorder=12)
axs[2,0].fill_between(time, post_aU*np.ones((time.shape[0])), kde_upper_post(time)*400 + post_aU, color=postcolorlite, alpha=.5, zorder=7)
axs[2,0].fill_between(time, post_aL*np.ones((time.shape[0])), -kde_lower_post(time)*400 + post_aL, color=postcolorlite, alpha=.5, zorder=8)


#Plot the evidence traces
allevidence_pre[allevidence_pre > pre_aU] = pre_aU #Correction for plotting first time pass distribution
allevidence_post[allevidence_post > post_aU] = post_aU #Correction for plotting first time pass distribution
allevidence_pre[allevidence_pre < pre_aL] = pre_aL #Correction for plotting first time pass distribution
allevidence_post[allevidence_post < post_aL] = post_aL #Correction for plotting first time pass distribution
axs[2,0].plot(time + pre_ndt_ms, allevidence_pre[0:20,0:(time.shape[0])].T, linewidth=.5,color=precolorlite, alpha=0.25, zorder=1)
axs[2,0].plot(time + post_ndt_ms, allevidence_post[0:20,0:(time.shape[0])].T, linewidth=.5,color=postcolorlite, alpha=0.25, zorder=2)


#Appropriately label plots
# axs[2,0].set_ylabel('Relative evidence for a choice to the IF', fontsize=fontsize)
axs[2,0].set_xlim([250,maxx]) #X limits in ms
axs[2,0].set_yticks(np.array([-1.0, 0, 1.0]))
# plt.tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[2,0].text(0.02, .98, 'i', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[2,0].transAxes)
axs[2,0].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[2,0].spines['right'].set_visible(False)
axs[2,0].spines['top'].set_visible(False)
axs[2,0].get_xaxis().set_visible(False)
axs[2,0].yaxis.set_ticks_position('left')
axs[2,0].xaxis.set_ticks_position('bottom')

#Plot the simulated psychometric curves
evidenceright = np.array([-36, -24, -17, -10, -5, -3, 0, 3, 5, 10, 17, 24, 36])
psychometric_xtick = np.array([-36, -17, -5, 0, 5, 17, 36])

psychometric_pre = np.ones((nconditions))*np.nan
psychometric_post = np.ones((nconditions))*np.nan
meanRT_correct_pre = np.ones((nconditions))*np.nan
meanRT_correct_post = np.ones((nconditions))*np.nan
meanRT_incorrect_pre = np.ones((nconditions))*np.nan
meanRT_incorrect_post = np.ones((nconditions))*np.nan
psychoindex = 0 
for n in conditionspan_offset:
    (alldata_pre, allevidence_pre) = fixedstartDDM(zeta=pre_zeta, drift=(n+pre_drift), aU=pre_aU, aL=pre_aL, dcoef=1, eta=eta, dt=dt, n=5000, maxsample=5000)
    psychometric_pre[psychoindex] = np.sum(alldata_pre > 0)/alldata_pre.shape[0]
    if (psychoindex < np.round(nconditions/2)):
        meanRT_correct_pre[psychoindex] = np.mean(np.abs(alldata_pre[alldata_pre < 0])) + pre_ndt_ms/1000
        meanRT_incorrect_pre[psychoindex] = np.mean(alldata_pre[alldata_pre > 0]) + pre_ndt_ms/1000
    elif (psychoindex > np.round(nconditions/2)):
        meanRT_correct_pre[psychoindex] = np.mean(alldata_pre[alldata_pre > 0]) + pre_ndt_ms/1000
        meanRT_incorrect_pre[psychoindex] = np.mean(np.abs(alldata_pre[alldata_pre < 0])) + pre_ndt_ms/1000
    else:
        meanRT_correct_pre[psychoindex] = np.mean(np.abs(alldata_pre)) + pre_ndt_ms/1000
        meanRT_incorrect_pre[psychoindex] = np.mean(np.abs(alldata_pre)) + pre_ndt_ms/1000        
    (alldata_post, allevidence_post) = fixedstartDDM(zeta=post_zeta, drift=(n+post_drift), aU=post_aU, aL=post_aL, dcoef=1, eta=eta, dt=dt, n=5000, maxsample=5000)
    psychometric_post[psychoindex] = np.sum(alldata_post > 0)/alldata_post.shape[0]
    if (psychoindex < np.round(nconditions/2)):
        meanRT_correct_post[psychoindex] = np.mean(np.abs(alldata_post[alldata_post < 0])) + pre_ndt_ms/1000
        meanRT_incorrect_post[psychoindex] = np.mean(alldata_post[alldata_post > 0]) + pre_ndt_ms/1000
    elif (psychoindex > np.round(nconditions/2)):
        meanRT_correct_post[psychoindex] = np.mean(alldata_post[alldata_post > 0]) + pre_ndt_ms/1000
        meanRT_incorrect_post[psychoindex] = np.mean(np.abs(alldata_post[alldata_post < 0])) + pre_ndt_ms/1000
    else:
        meanRT_correct_post[psychoindex] = np.mean(np.abs(alldata_post)) + pre_ndt_ms/1000
        meanRT_incorrect_post[psychoindex] = np.mean(np.abs(alldata_post)) + pre_ndt_ms/1000

    psychoindex += 1




axs[2,1].plot(np.array([-50, 50]), np.array([50, 50]), color='k', LineStyle='--')
axs[2,1].plot(np.array([0, 0]), np.array([0, 100]), color='k', LineStyle='--')

popt_pre, pcov_pre = curve_fit(fpsychometric, evidenceright, psychometric_pre)
axs[2,1].plot(evidenceright, psychometric_pre*100, 'o',color=precolor, markersize=markersize)
axs[2,1].plot(evidenceright, fpsychometric(evidenceright, *popt_pre)*100, color=precolor, linewidth = 1)
popt_post, pcov_post = curve_fit(fpsychometric, evidenceright, psychometric_post)
axs[2,1].plot(evidenceright, psychometric_post*100, 'o',color=postcolor, markersize=markersize)
axs[2,1].plot(evidenceright, fpsychometric(evidenceright, *popt_post)*100, color=postcolor, linewidth = 1)

axs[2,1].set_xlabel('Coherence (%)', fontsize=fontsize)
axs[2,1].set_xticks(psychometric_xtick)
axs[2,1].set_xticklabels(psychometric_xtick)
axs[2,1].set_ylim((0, 100))
axs[2,1].set_xlim((-40,40))
axs[2,1].set_yticks(np.array([20, 40, 60, 80, 100]))
# axs[2,1].set_ylabel('Choices to IF (%)', fontsize=fontsize)
axs[2,1].text(0.02, .98, 'j', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[2,1].transAxes)
axs[2,1].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[2,1].spines['right'].set_visible(False)
axs[2,1].spines['top'].set_visible(False)
axs[2,1].get_xaxis().set_visible(False)
axs[2,1].yaxis.set_ticks_position('left')
axs[2,1].xaxis.set_ticks_position('bottom')


axs[2,2].plot(evidenceright,meanRT_correct_pre*1000, 'o', color=precolor, markersize=markersize)
axs[2,2].plot(evidenceright,meanRT_correct_pre*1000, linewidth=1,color=precolor, zorder=2)
axs[2,3].plot(evidenceright,meanRT_incorrect_pre*1000, 's', color=precolor, markersize=markersize)
axs[2,3].plot(evidenceright,meanRT_incorrect_pre*1000, color=precolor, linewidth=1, linestyle='--', zorder=1)
axs[2,2].plot(evidenceright,meanRT_correct_post*1000, 'o', color=postcolor, markersize=markersize)
axs[2,2].plot(evidenceright,meanRT_correct_post*1000, linewidth=1,color=postcolor, zorder= 4)
axs[2,3].plot(evidenceright,meanRT_incorrect_post*1000, 's', color=postcolor, markersize=markersize)
axs[2,3].plot(evidenceright,meanRT_incorrect_post*1000, color=postcolor, linestyle='--', linewidth=1, zorder=3)

axs[2,2].set_xlabel('Coherence (%)', fontsize=fontsize)
axs[2,2].set_xticks(psychometric_xtick)
axs[2,2].set_xticklabels(psychometric_xtick)
# axs[2,2].set_ylim((300,1000))
# axs[2,2].set_yticks(np.array([500, 700, 900]))
axs[2,2].set_ylim((500,1200))
axs[2,2].set_yticks(np.array([700, 900, 1100]))
# axes = plt.gca()
# ymin, ymax = axes.get_ylim()
# axs[2,2].plot(np.array([0, 0]), np.array([300,1000]), color='k', LineStyle='--')
axs[2,2].plot(np.array([0, 0]), np.array([500,1200]), color='k', LineStyle='--')
axs[2,2].text(0.02, .98, 'k', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[2,2].transAxes)
axs[2,2].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[2,2].spines['right'].set_visible(False)
axs[2,2].spines['top'].set_visible(False)
axs[2,2].get_xaxis().set_visible(False)
axs[2,2].yaxis.set_ticks_position('left')
axs[2,2].xaxis.set_ticks_position('bottom')

# axs[2,3].set_xlabel('Coherence (%)', fontsize=fontsize)
axs[2,3].set_xticks(psychometric_xtick)
axs[2,3].set_xticklabels(psychometric_xtick)
# axs[2,3].set_ylim((300,1000))
# axs[2,3].set_yticks(np.array([500, 700, 900]))
axs[2,3].set_ylim((500,1200))
axs[2,3].set_yticks(np.array([700, 900, 1100]))
# axes = plt.gca()
# ymin, ymax = axes.get_ylim()
# axs[2,3].plot(np.array([0, 0]), np.array([300,1000]), color='k', LineStyle='--')
axs[2,3].plot(np.array([0, 0]), np.array([500,1200]), color='k', LineStyle='--')
axs[2,3].text(0.02, .98, 'l', horizontalalignment='left', verticalalignment='top', weight='bold',fontsize=panelsize,transform=axs[2,3].transAxes)
axs[2,3].tick_params(axis='both', which='major', labelsize=fontsize, direction='in')
axs[2,3].spines['right'].set_visible(False)
axs[2,3].spines['top'].set_visible(False)
axs[2,3].get_xaxis().set_visible(False)
axs[2,3].yaxis.set_ticks_position('left')
axs[2,3].xaxis.set_ticks_position('bottom')



#For all subplots
# plt.subplots_adjust(hspace=0,left=0, right=0.1, bottom=0, top=0.1)
plt.subplots_adjust(hspace=0, wspace=0.4)
# plt.tight_layout(h_pad=-1,w_pad=0.3)
fig.set_size_inches(mm2inch(218,147),forward=False)

plt.savefig((f'../../figures/diffusion_figures/Figure4sims.png'), dpi=300, format='png',bbox_inches='tight')
plt.savefig((f'../../figures/diffusion_figures/Figure4sims.pdf'), dpi=300, format='pdf',bbox_inches='tight')
plt.savefig((f'../../figures/diffusion_figures/Figure4sims.svg'), dpi=300, format='svg',bbox_inches='tight')