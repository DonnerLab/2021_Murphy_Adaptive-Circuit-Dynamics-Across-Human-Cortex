function [] = runFreqPrepMotorHiFreq(n)
  % Apply planar gradient transformation and calculate TF representations of
  % MEG power for each of two gradiometers per sensor/trial, using a single
  % Hanning taper for freq range 3-35Hz (window length: 400 ms, step size:
  % 50 ms, freq resolution: 2.5 Hz, bin size: 1 Hz). After TFR calculation,
  % power estimates from the two planar gradiometers for each sensor are
  % combined by taking their sum.

  
% ==================================================================
% SPECIFY PATHS AND GET SUBJECT-SPECIFIC FILES
% ==================================================================
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
ft_defaults 

megpath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/PreprocessedMotor/Data/';  % path of preprocessed MEG data
savepath = '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/PreprocessedMotor/';

% Pull desired MEG dataset
files = dir([megpath,'*_PreprocessedMotor.mat']);
subjfile = files(n).name;  % pull current meg filename
ID = subjfile(1:5);

% ==================================================================
% LOAD CURRENT MEG DATA
% ==================================================================
fprintf('\nLoading meg file: %s...\n',subjfile)
load([megpath,subjfile])  % load meg data

% ==================================================================
% PULL ONLY FULL-LENGTH TRIALS - consider skipping this to maximize trial counts
% ==================================================================
fprintf('Keeping only full-length trials...\n')
assert(length(data.trial)==size(trl,1),'ERROR: Mismatch in MEG/behaviour number of trials')
ts=find(~isnan(trl(:,4)));  % useable trials based on behavioural data

resps = trl(ts,4);
RTs = trl(ts,5);

cfg             = [];
cfg.trials      = ts;
data = ft_selectdata(cfg, data);

% ==================================================================
% PLANAR GRADIENT TRANSFORMATION
% ==================================================================
fprintf('\nRunning planar gradient transformation...\n')

% define neighbours based on CTF template
cfg                 = [];
cfg.method          = 'template';
cfg.layout          = 'CTF275';
neighbours          = ft_prepare_neighbours(cfg);

% compute planar gradiometers for MEG sensors
cfg                 = [];
cfg.feedback        = 'no';
cfg.method          = 'template';
cfg.planarmethod    = 'sincos';
cfg.channel         = 'MEG';
cfg.neighbours      = neighbours;
data                = ft_megplanar(cfg, data);

% ==================================================================
% TIME-FREQUENCY DECOMPOSITION
% ==================================================================
fprintf('\nRunning time-frequency decomposition...\n')

cfg                 = [];
cfg.output          = 'pow';
cfg.channel         = 'MEG';
cfg.method          = 'mtmconvol';   % specifying multi-taper method
cfg.taper           = 'dpss';     % with Hanning taper
cfg.keeptrials      = 'yes';
cfg.keeptapers      = 'no';
cfg.precision       = 'single'; % saves disk space
%cfg.feedback        = 'none'; % improve readability of logfiles

% make nice timebins at each 50 ms, will include 0 point of pre-mask onset; last time point will correspond to go cue time for trial with shortest dot-to-go interval
mintime = data.time{1}(1);
maxtime = data.time{1}(end);
toi = floor(mintime) : 0.05 : ceil(maxtime);
toi(toi < mintime) = []; toi(toi > maxtime) = [];

cfg.toi             = toi; % all time within each locking
cfg.pad             = 4; % pad to a fixed number of seconds before TFR

cfg.foi             = 36:4:120;   % frequencies of interest
cfg.t_ftimwin       = ones(1, length(cfg.foi)) .* 0.25;
% cfg.t_ftimwin       = 3*(1./cfg.foi);  % adaptive time window of 3 cycles per estimated frequency
cfg.tapsmofrq       = ones(1, length(cfg.foi)) .* 6;  % specifies half the spectral concentration, which should be between 3 and 11 times the Rayleigh frequency (1/t_ftimwin)

freq                 = ft_freqanalysis(cfg, data);

% ==================================================================
% COMBINE PLANAR GRADIOMETERS
% ==================================================================
freq = ft_combineplanar([], freq);

% ==================================================================
% SAVE
% ==================================================================
freq.resps = resps;
freq.RT = RTs;

save([savepath,ID,'_TF_HiFreq.mat'],'freq')
    
end