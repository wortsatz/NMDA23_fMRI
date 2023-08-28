
%% ============= NMDA 23 - EEG MASTER BATCH SCRIPT with BIDS data =======================

% this script is the master batch script for EEG


%% -------------------- Parameter Changes -------------------------------

% Select which path is the basis Path
whichComp=1;

if whichComp==1
    spmPath='/Users/ttli/Documents/spm12/';
    bidsPath = '/Users/ttli/Documents/GitHub/NMDA23/EEG/bids';
    project_path='/Users/ttli/Documents/GitHub/NMDA23/';
elseif whichComp==2
    spmPath='/Users/kubrafatulla/Documents/spm12';
    bidsPath = '/Users/kubrafatulla/Documents/GitHub/NMDA23/EEG/bids';
    project_path='/Users/kubrafatulla/Documents/GitHub/NMDA23/';
elseif whichComp==3 
    spmPath='/Users/lp1/Documents/spm12/';
    bidsPath = '/Users/lp1/Documents/GitHub/NMDA23/EEG/bids';
    project_path='/Users/lp1/Documents/GitHub/NMDA23/';
else
    spmPath='/Users/Nastya/программы/spm12/';
    bidsPath = '/Users/Nastya/Documents/GitHub/NMDA23/EEG/bids';
    project_path='/Users/Nastya/Documents/GitHub/NMDA23/';
end
cd(data_path)
addpath((spmPath))
addpath(genpath(project_path))

% Define BIDS variables
BIDS = bids.layout(bidsDir);
subjectList=BIDS.subjects;
sessionList = BIDS.sessions;

% Iterate over subjects and sessions
for s = 1:length(subjectList)
    for ses = 1:length(sessionList)
        subject = subjectList{s};
        session = sessionList{ses};
  
        % Construct paths to your BIDS-formatted data
        data_path = fullfile(bidsDir, '/sub-',num2str(s), '/session-',num2str(ses),'eeg/');

        % Load and process data using SPM
        jobs = EEG_preprocessing_batch_job(data_path, spmPath, length(subjectList), length(sessionList));
        spm('defaults', 'EEG');
        spm_jobman('run', jobs);
        
     
    end
end







