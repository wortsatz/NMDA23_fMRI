
%% ============= NMDA 23 - EEG MASTER BATCH SCRIPT =======================

% this script is the master batch script for EEG


%% -------------------- Parameter Changes -------------------------------

% Select which path is the basis Path
whichComp=1;

if whichComp==1
    spmPath='/Users/ttli/Documents/spm12/';
    data_path='/Users/ttli/Documents/GitHub/NMDA23/98_fMRI_data';
    project_path='/Users/ttli/Documents/GitHub/NMDA23/';
elseif whichComp==2
    spmPath='/Users/kubrafatulla/Documents/spm12';
    data_path='/Users/kubrafatulla/Documents/fMRI/MoAEpilot/';
    project_path='/Users/kubrafatulla/Documents/GitHub/NMDA23/';
elseif whichComp==3 
    spmPath='/Users/lp1/Documents/spm12';
    data_path='/Users/lp1/Nextcloud/FU/NMDA_data';
    project_path='/Users/lp1/Documents/GitHub/NMDA23/';
else
    spmPath='/Users/Nastya/программы/spm12';
    data_path='/Users/Nastya/Documents/Freie Universitat/NMDA';
    project_path='/Users/Nastya/Documents/GitHub/NMDA23/';
end
cd(data_path)
addpath((spmPath))
addpath(genpath(project_path))



%% ------------------------- Preprocessing --------------------------------


% List of open inputs
nrun = 1; % enter the number of runs here
preprocessing_jobfile = {'C:\Users\Nastya\Documents\Freie Universitat\EEG_preprocessing_batch_job.m'};
jobs = preprocessing_jobfile;
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'EEG');
spm_jobman('run', jobs, inputs{:});
