# SC Inactivation Project
#### (Repository version 0.1.1)

### Citation

Jun, E. J., Bautista, A., Nunez, M. D., Allen, D. C., Tak, J. H., Alvarez, E., Basso, M. A. (2021). [Causal role for the primate superior colliculus in the computation of evidence for perceptual decisions](https://www.nature.com/articles/s41593-021-00878-6.epdf?sharing_token=alTJNGt7dLdwsp-jZk6OWNRgN0jAjWel9jnR3ZoTv0PCntPI-NBEHPcVRFJdN8ENsUE-yHy7z3ccBLF113IPsFSXPgWX8PXvKIPeExS8kFVrEggTxLN5s5xgpQEBOto7oT01GkMMdTMh5tVpNNkYxk2iUAu7AXcI-X8k87E23vc%3D). Nature Neuroscience. 24, 1121–1131.

## Getting Started

### Prerequisites

[Python 3 and Scientific Python libraries](https://www.anaconda.com/products/individual)

[MCMC Sampling Program: JAGS](http://mcmc-jags.sourceforge.net/)

[Program: JAGS Wiener module](https://sourceforge.net/projects/jags-wiener/)

For JAGS and JAGS Wiener module install steps in Ubuntu, see [this document](https://github.com/mdnunez/pyhddmjags/blob/master/jags_wiener_ubuntu.md).

[Python Repository: pyjags](https://github.com/michaelnowotny/pyjags), can use pip:
```bash
pip install pyjags
```

[R](https://www.r-project.org/)

[ChartR](https://github.com/mailchand/CHaRTr)

[Palamedes](http://www.palamedestoolbox.org/)

### Processing Data in MATLAB
For processing/analyzing decision task data -- 

Process raw data for each session and plot psychometric functions, mean RTs, RT distributions: run decisiontask_masterscript_plotter.m

Processes aggregate data across injection sessions and plot aggregate psychometric functions, mean RTs, accuracies for toIF and awayIF, and peak velocities:
-	for muscimol injections --> run decision_task_allmuscimolinjections.m
-	for saline injections --> run decisiontask_allsalineinjections.m


For processing/analyzing selection task data -- 

Process raw data for each session and plot accuracies for toIF and awayIF choices: run selectiontask_masterscript_plotter.m

Process aggreggate data across injection sessions and plot aggreggate accuracies and peak velocities for toIF and awayIF choices: run selectiontask_allinjections.m


For procesing/analyzing visually guided saccade task data --

To process the raw data into LMData: run VGST_rawdataprocessor.m 

To get the heat map of difference in peak velocity (post - pre) in the format of visual field in cartesian coordinates, in polar coordinates, and on the SC map: run pv_maps_and_analysis.m

Analysis of aggregate data, collapsed across injection sessions: run VGST_task_all_injections.m



### Downloading

The repository can be cloned with `git clone https://gitlab.com/fuster-lab-cognitive-neuroscience/sc-inactivation-project.git`

The repository can also be may download via the _Clone_ button above.


## License

SCInject is licensed under the GNU General Public License v3.0 and written by Michael D. Nunez and Elizabeth J. Jun from the Department of Psychiatry and Biobehavioral Sciences at the University of California, Los Angeles (UCLA).

