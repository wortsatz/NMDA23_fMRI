%% ===================== EEG TO BIDS SCRIPT ==============================

%% Path and Parameters

data_path='/Users/ttli/Documents/GitHub/NMDA23/EEG/EEG_dat/';
%data_path='/Users/ttli/Documents/GitHub/NMDA23/EEG/';


sub = {'1'};
age = [11];
sex = {'f'};

%% Conversion 

for subindx=1:numel(sub)

  cfg = [];
  cfg.method    = 'copy';
  cfg.datatype  = 'eeg';
  
  fileName=[data_path,'spmeeg_subject',num2str(subindx),'.mat'];
  %fileName=[data_path,'subject',num2str(subindx),'.bdf'];

  % specify the input file name, here we are using the same file for every subject
  cfg.dataset   = fileName;

  % specify the output directory
  cfg.bidsroot  = 'bids3';
  cfg.sub       = sub{subindx};
  

  % specify the information for the participants.tsv file
  % this is optional, you can also pass other pieces of info
  cfg.participants.age = age(subindx);
  cfg.participants.sex = sex{subindx};

  % specify the information for the scans.tsv file
  % this is optional, you can also pass other pieces of info
  cfg.scans.acq_time = datestr(now, 'yyyy-mm-ddThh:MM:SS'); % according to RFC3339

  % specify some general information that will be added to the eeg.json file
  cfg.InstitutionName             = 'Freie University Berlin';
  cfg.InstitutionalDepartmentName = 'Cognitive Neuroscience';
 
  % provide the mnemonic and long description of the task
  cfg.TaskName        = 'spm12MMN';
  cfg.TaskDescription = '128-channel EEG data set acquired from a study of mismatch negativity in the auditory system. The experiment comprised an auditory oddball paradigm in which subjects heard standard (500Hz) and deviant (550Hz) tones, occuring 80% (480 trials) and 20% (120 trials) of the time, respectively, in a pseudo-random sequence subject to the constraint that two deviant tones did not occur together.';

  % these are EEG specific
  cfg.eeg.PowerLineFrequency = 60;   % since recorded in the USA
  cfg.eeg.EEGReference       = 'C1'; % actually I do not know
  
  data2bids(cfg);

end



%%

% Define paths
%bdf_file = '/Users/ttli/Documents/GitHub/NMDA23/EEG/EEG_dat/subject1.bdf';
bdf_file='/Users/ttli/Documents/GitHub/NMDA23/EEG/EEG_dat/spmeeg_subject1.mat';
%data=load(bdf_file);

bids_root = '/Users/ttli/Documents/GitHub/NMDA23/EEG/bids';
subject_id = '1';
session_label = '1';
task_label = 'mmn';


% Load the .bdf data using FieldTrip
cfg = [];
cfg.dataset = bdf_file;
data = ft_preprocessing(cfg);

% Create the BIDS directory structure
bids_dir = fullfile(bids_root, ['sub-' subject_id], ['ses-' session_label], 'eeg');
mkdir(bids_dir);

% Save the EEG data in BIDS format
bids_file = fullfile(bids_dir, ['sub-' subject_id '_ses-' session_label '_task-' task_label '_eeg.bdf']);
cfg = [];
cfg.datafile = bdf_file;
cfg.headerfile = bdf_file;
cfg.dataset = bids_file;
ft_write_data(cfg.datafile, data);

% Create and save the JSON sidecar file
json_data = struct('TaskDescription', 'Description of your task');
json_file = fullfile(bids_dir, ['sub-' subject_id '_ses-' session_label '_task-' task_label '_eeg.json']);
savejson('', json_data, json_file);















