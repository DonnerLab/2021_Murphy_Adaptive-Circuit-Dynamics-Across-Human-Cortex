## 2021-Adaptive-Circuit-Dynamics-Across-Human-Cortex
Task and analysis code for

[**Murphy PR, Wilming N, Hernandez Bocanegra DC, Prat Ortega G & Donner TH (2021). Adaptive circuit dynamics across human cortex during evidence accumulation in changing environments. bioRxiv.**](https://www.biorxiv.org/content/10.1101/2020.01.29.924795v4)

Raw behavioural and eye-tracking data are available at [WWW](). Raw MEG data are available at [XXX](). Source reconstructed MEG data are available at [YYY](). Source data for all main text figures and code to reproduce them are available at [ZZZ](). All are shared under a CC-BY 4.0 license.
Code shared here was developed and tested using Matlab R2015a, [FieldTrip Toolbox](https://www.fieldtriptoolbox.org/) version 20160221, Python 3.6 and [mne](https://mne.tools/stable/index.html) version 0.16.2.

Further detail on analysis is provided below. For questions, contact murphyp7@tcd.ie.

#### behav_modelling:
Scripts for fitting a variety of illustrative model variants to choice data, some basic (e.g. `Glaze_basic` for normative model with three free parameters) and some more complex (e.g. `FitPsi_data_npLLR_InconUp` for a model in which both stimulus->log-likelihood ratio and posterior belief->next-sample prior are estimated as interpolated functions, and which includes a term reflecting over-/under-weighting of new evidence that is inconsistent with the existing belief). Code for each model variant is contained in individual subdirectories.
-	Models are fit via particle swarm optimization ([Birge, 2003, *IEEE Swarm Intelligence Symposium*](https://ieeexplore.ieee.org/document/1202265); [code here](https://www.mathworks.com/matlabcentral/fileexchange/7506-particle-swarm-optimization-toolbox)), with very slightly adapted version of code included in the `particle_swarm_Glaze` subdirectory (adaptation made code work with global variables passed to it via higher-level scripts, and is commented with `% PM`).
-	`pull_behav_fits.m`: example script that loads participant behavioural data and associated model fits and computes variety of measures (accuracy, psychophysical kernels, etc.)
#### circuit_modelling:
-	`run_trials.m`: top-level script that specifies circuit and task parameters of interest, loads stimulus inputs provided to participants (and provided here in `norm_vars.mat`), and simulates 50 trials of the circuit model provided with these inputs.
#### motor_localizer:
Scripts for preprocessing sensor-level MEG data from both decision-making and motor localizer tasks, estimating spatio-spectral motor preparatory filter weights from the localizer task, and applying these filters in analyses of data from the decision-making task.
-	`MEGpreprocessing.m` and `MEGpreprocessingMotor.m`: preprocess raw MEG data and reject trials with artifacts.
-	`runFreqPrep*.m`: extract behavioural/model-derived/pupillometric variables for each clean MEG trial and compute time-frequency representations of MEG data.
-	`runFreqAnalysisMotor*.m`: computes hemisphere-lateralized MEG power from the motor localizer data, over time and frequency.
-	`plotTF_clustercorr_LIhalf_AllFreqs.m`: plots lateralization signal during delay period of motor localizer task and constructs spatio-spectral filters to be applied to decision-making task.
-	`runFreqAnalysis_appMotorLocalizer*.m`: applies motor filters to decision-making task data in a variety of ways (single-trial regressions using model-derived variables; computing ‘belief encoding functions’; binning of trials by timing of change-points).
-	End product: `superBIC_appMotorLoc.m` to compute ‘super BIC scores’ comparing different linear models of the motor signal; `plot_CPcentered_tonly.m` to plot motor signal for trials binned by change-point position; `cluster_corr*.m` to run cluster-corrected statistical tests on effects of interest.
#### source_reconstruct:
Scripts for preprocessing MEG data for the decision-making task specifically in preparation for mne- python-based source reconstruction routine; running source reconstruction; and conducting various analyses of source-reconstructed data (see below). For more detailed step-by-step instructions for running source reconstruction routine itself, see included `source_recon_steps.docx`.
-	`BehavPupilPreproc4mne.m`: creates structures containing behavioural, model-derived and pupillometric measures that align with trials used in source reconstruction.
-	`pymeg\lcmv_decoding_PM.py` and `pymeg\lcmv_decoding_ERF_PM.py`: decode variables of interest from ROI-specific source reconstructed data, see.
-	`pymeg\fit_lin_fooof.py`: estimate power spectra exponents (via the [FOOOF toolbox]( https://fooof-tools.github.io/fooof/index.html); [Donoghue et al., 2020, *Nature Neuroscience*](https://www.nature.com/articles/s41593-020-00744-x)).
-	`runSRanalysis_LIhalf.m`: example script for loading ROI-specific source reconstructed data (in Matlab) and regressing single-trial estimates onto computational variables of interest. See also `runSRanalysis_LIhalf_PPI.m` for analysis of effects of signal residuals on behavioural choice; and `clust_stat_F_LIhalf.m` for cluster-corrected statistical testing of the kind of output generated by these analyses. 
-	`plot_ERFdec_wperm_subclust.m` and `plot_LLRenc_wperm_subclust.m`: scripts for computing latency to half maxima and exponential decay timescales of i) decoded *LLR* from event-related field responses and ii) *LLR* encoding in alpha-band power, respectively. Compares estimates between ROI clusters via weighted permutation tests. ERF decoding runs on output from `lcmv_decoding_ERF_PM.py`; *LLR* encoding runs on output from `runSRanalysis_LIhalf.m`.
#### pupil:
Scripts for preprocessing and analysing pupillometric and eye-tracking data.
-	`loop_proc_asc.m`: read eye-tracker data into Matlab and pull timing of blinks and saccades. Assumes raw .edf data files from the Eyelink 1000 (SR Research) have been converted into .asc text files via the manufacturer’s ‘edf2asc’ software utility.
-	`loop_interp.m`: interpolates pupil and gaze position time-series.
-	`analyse_pupil.m`: example script for analysing pupil data. In this case, script pulls and plots trial-averaged pupil responses, regresses pupil diameter onto computational variables of interest in a time-resolved fashion, and runs ‘psychophysiological interaction’ style choice regressions to quantify the impact that residual pupil fluctuations have on the weighting of evidence samples in choice.
#### task_decision and task_motor:
Scripts that use [Psychotoolbox-3](http://psychtoolbox.org/) to run the decision-making and motor-localizer tasks.
-	Decision-making task: run training routine with `Surprise_radial_checkers_intro.m` and main task with `Surprise_radial_checkers_block.m`.
-	Motor localizer task: run with `Motor_localizer_block.m`.
