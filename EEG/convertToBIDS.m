sub = {'1'};
data_path='/Users/ttli/Dropbox/FreieU/EEG_dat/';
FT_path='/Users/ttli/Dropbox/fieldtrip-20230118/';

addpath(genpath(FT_path))

for subindx=1:numel(sub)

  cfg = [];
  cfg.method    = 'copy';
  cfg.datatype  = 'eeg';
  
  fileName=[data_path,'spmeeg_subject',num2str(subindx),'.mat'];

  % specify the input file name, here we are using the same file for every subject
  cfg.dataset   = fileName;

  % specify the output directory
  cfg.bidsroot  = 'bids';
  cfg.sub       = sub{subindx};


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