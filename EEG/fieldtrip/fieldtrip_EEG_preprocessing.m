%% =================== FIELDTRIP EEG PREPROCESSING =======================

% This script lets us do preprocessing using the fieldtrip functions


%% ------------------------ Set paths -------------------------------------

% this script uses .bdf for now - a extension for BIDS is in the works

data_path = '/Users/lp1/Documents/GitHub/NMDA23/EEG/SPM/subject1.bdf';
sensor_path = '/Users/lp1/Documents/GitHub/NMDA23/EEG/SPM/sensors.pol';
montage_path = '/Users/lp1/Documents/GitHub/NMDA23/EEG/SPM/avref_eog_montage.mat';
trialdef_path = '/Users/lp1/Documents/GitHub/NMDA23/EEG/SPM/trialdef.mat';
spmdir = '/Users/lp1/Documents/spm12';

% add data and SPM patht to the Matlab dir
addpath(genpath(data_path), genpath(spmdir));

%% ------------------------ Set Basic Parameters --------------------------

% Number of Subj
sub = 1;

% Here you clarify which Preprocessing steps you want to take
% 1 = Montage
% 2 = High-Pass Filter
% 3 = Downsample
% 4 = Low-Pass Filter
% 5 = Epoch
% 6 = Artefacts
prepro_switch = [1,2,3,4]; 

% Initialize FieldTrip
ft_defaults; 

%% --------------------------- Preprocessing ------------------------------

for s = 1:numel(sub)
    % Load data
    cfg = [];
    cfg.dataset = data_path;
    data = ft_preprocessing(cfg);

    for n = prepro_switch
        switch n
            case 1 % Montage
                % Montage the EEG data to the sensor locations of the EEG
                % Device 

                % Load & Apply montage to the channels
                load(montage_path);
                cfg = [];
                cfg.montage.labelnew = montage.labelnew;
                cfg.montage.labelorg = montage.labelorg;
                cfg.montage.tra = montage.tra;
                cfg.channel = 'all'; 
                data = ft_preprocessing(cfg, data);

                % Load sensor layout and add channal labels 
                cfg.layout = sensor_path;
                layout = ft_prepare_layout(cfg);
                channel_labels = layout.label;

            case 2 % High Pass Filter
                % High Pass Filter everything above 1 Hz and save it
                cfg = [];
                cfg.hpfilter = 'yes';
                cfg.hpfreq = 1;
                cfg.hpfilttype = 'firws'; % Filter design
                data = ft_preprocessing(cfg, data);

            case 3 %  Downsampling
                % Downsample the data with 200 samples and store it
                cfg = [];
                cfg.resamplefs = 200;
                data = ft_resampledata(cfg, data); 


            case 4 % Filtering (low-pass)
                % Low Pass Filter so that everything above 60 Hz is
                % disregarded

                cfg = [];
                cfg.lpfilter = 'yes';
                cfg.lpfreq = 60;
                cfg.lpfilttype = 'firws'; % Filter design
                data = ft_preprocessing(cfg, data);


            case 5 % Epoching
                % Epoch the data around the trials and save the epochs

                trialdef = load(trialdef_path);
                cfg = [];
                cfg.trl = trialdef.trl;
                data = ft_redefinetrial(cfg, data);


            case 6 % Artifacts
                % Do Artifact removal of bad channels and bad trials 


                % Visualise data and reject bad channels 
                cfg.method = 'summary';
                cfg.keepchannel = 'no';
                cfg.keeptrial = 'nan';
                cfg.summary.channel = 'all';
                cfg.summary.threshold = 80;
                data = ft_rejectvisual(cfg, data);

                % Delete the rejected trials also in the trialdef.conditionlabels file
                rejected_trials = [46 60 61 65 66 67 68 100 165 224 237]; % indices of the rejected trials
                for i = 1:length(rejected_trials)
                    trialdef.conditionlabels{rejected_trials(i)} = []; % replace rejected trials with an empty array
                end
                
                % Find and delete the empty rows 
                trialdef.conditionlabels(cellfun(@isempty,trialdef.conditionlabels))=[]; 
             
        end

    end

    % save the preprocessed data for each sub and after all steps were done
    cleaned_data = data;

    save('cleaned data.mat', 'cleaned_data');

end

