%% setting up the environment for using SPM 12
which spm
addpath('/Users/kubrafatulla/Documents/fMRI/MoAEpilot/')
spm
%% 
% Matlab script for realignment

% Load SPM12 and initialize paths
spmDir = '/Users/kubrafatulla/Documents/spm12';
addpath(spmDir);

% Load fMRI data and specify realignment options
dataDir = '/Users/kubrafatulla/Documents/fMRI/MoAEpilot/fM00223';
funcFiles = spm_select('FPList', dataDir, '^fM00223_.*\.img$');
realignment = struct('quality', 0.9, 'sep', 4, 'fwhm', 5);

% Run realignment using SPM12
spm_realign(funcFiles, realignment);

% Save realignment parameters to file
save(fullfile(dataDir, 'realign_params.mat'), 'realignment');