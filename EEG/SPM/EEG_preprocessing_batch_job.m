function [] = EEG_preprocessing_batch_job(datapath, spmPath)

%--------------------------JOB SCRIPT FOR EEG PREPROCESSING--------------------------------

% be aware that all the raw data and the created files should be stored in
% one folder on your computer


% Number of Subj
sub = 1;

% Number of runs
run = 1;

% Here you clarify which Preprocessing Step you want to take
% 1 = Convert
% 2 = Montage
% 3 = Prepare
% 4 = High-Pass Filter
% 5 = Downsample
% 6 = Low-Pass Filter
% 7 = Epoch
% 8 = Artefacts
% 8 = Averaging

prepro_switch = [1,2,3,4]; 

% predetermine the filename
bdf_filename = {'subject1.bdf'};
spm_filename = {'spmeeg_subject1.mat'};


% give the initial prefix
prefix = '';


%% -------------------------- PreProcessing ------------------------------

for s = 1:numel(sub)
    for r = 1:numel(run)
        for n = prepro_switch
            switch n
                case 1 % Convert
                    matlabbatch{1}.spm.meeg.convert.dataset = [ datapath prefix bdf_filename];
                    matlabbatch{1}.spm.meeg.convert.mode.continuous.readall = 1;
                    matlabbatch{1}.spm.meeg.convert.channels{1}.chanfile = [datapath 'channelselection.mat'];
                    matlabbatch{1}.spm.meeg.convert.outfile = '';
                    matlabbatch{1}.spm.meeg.convert.eventpadding = 0;
                    matlabbatch{1}.spm.meeg.convert.blocksize = 3276800;
                    matlabbatch{1}.spm.meeg.convert.checkboundary = 1;
                    matlabbatch{1}.spm.meeg.convert.saveorigheader = 0;
                    matlabbatch{1}.spm.meeg.convert.inputformat = 'autodetect';

                case 2 % Montage

                    matlabbatch{2}.spm.meeg.preproc.montage.D = [datapath prefix spm_filename];
                    matlabbatch{2}.spm.meeg.preproc.montage.mode.write.montspec.montage.montagefile = [datapath 'avref_eog_montage.mat'];
                    matlabbatch{2}.spm.meeg.preproc.montage.mode.write.montspec.montage.keepothers = 0;
                    matlabbatch{2}.spm.meeg.preproc.montage.mode.write.blocksize = 655360;
                    matlabbatch{2}.spm.meeg.preproc.montage.mode.write.prefix = 'M';
                    prefix = ['M' prefix]; % updates the prefix to show that the Montage was done

                case 3 % prepare

                    matlabbatch{3}.spm.meeg.preproc.prepare.D = [datapath prefix spm_filename];
                    matlabbatch{3}.spm.meeg.preproc.prepare.task{1}.loadeegsens.eegsens = [datapath 'sensors.pol'];
                    matlabbatch{3}.spm.meeg.preproc.prepare.task{1}.loadeegsens.megmatch.nomatch = 1;

                case 4 % Highpass Filter

                    matlabbatch{4}.spm.meeg.preproc.filter.D = [datapath prefix spm_filename];
                    matlabbatch{4}.spm.meeg.preproc.filter.type = 'butterworth';
                    matlabbatch{4}.spm.meeg.preproc.filter.band = 'high';
                    matlabbatch{4}.spm.meeg.preproc.filter.freq = 0.1;
                    matlabbatch{4}.spm.meeg.preproc.filter.dir = 'twopass';
                    matlabbatch{4}.spm.meeg.preproc.filter.order = 5;
                    matlabbatch{4}.spm.meeg.preproc.filter.prefix = 'f';
                    prefix = ['f' prefix]; % updates the prefix

                case 5 % Downsample

                    matlabbatch{5}.spm.meeg.preproc.downsample.D = [datapath prefix spm_filename];
                    matlabbatch{5}.spm.meeg.preproc.downsample.fsample_new = 200;
                    matlabbatch{5}.spm.meeg.preproc.downsample.method = 'resample';
                    matlabbatch{5}.spm.meeg.preproc.downsample.prefix = 'd';
                    prefix = ['d' prefix]; % updates the prefix

                case 6 % Lowpass filter

                    matlabbatch{6}.spm.meeg.preproc.filter.D = [datapath prefix spm_filename];
                    matlabbatch{6}.spm.meeg.preproc.filter.type = 'butterworth';
                    matlabbatch{6}.spm.meeg.preproc.filter.band = 'low';
                    matlabbatch{6}.spm.meeg.preproc.filter.freq = 30;
                    matlabbatch{6}.spm.meeg.preproc.filter.dir = 'twopass';
                    matlabbatch{6}.spm.meeg.preproc.filter.order = 5;
                    matlabbatch{6}.spm.meeg.preproc.filter.prefix = 'f';
                    prefix = ['f' prefix]; % updates the prefix

                case 7 % epoch

                    matlabbatch{7}.spm.meeg.preproc.epoch.D = [datapath prefix spm_filename];
                    matlabbatch{7}.spm.meeg.preproc.epoch.trialchoice.trlfile = [datapath 'trialdef.mat'];
                    matlabbatch{7}.spm.meeg.preproc.epoch.bc = 1;
                    matlabbatch{7}.spm.meeg.preproc.epoch.eventpadding = 0;
                    matlabbatch{7}.spm.meeg.preproc.epoch.prefix = 'e';
                    prefix = ['e' prefix]; % updates the prefix

                case 8 % artefacts

                    matlabbatch{8}.spm.meeg.preproc.prepare.D = [datapath prefix spm_filename];
                    matlabbatch{8}.spm.meeg.preproc.prepare.task{1}.setbadchan.channels{1}.chan = 'A14';
                    matlabbatch{8}.spm.meeg.preproc.prepare.task{1}.setbadchan.status = 1;
                    matlabbatch{9}.spm.meeg.preproc.artefact.D = [datapath prefix spm_filename];
                    matlabbatch{9}.spm.meeg.preproc.artefact.mode = 'mark';
                    matlabbatch{9}.spm.meeg.preproc.artefact.badchanthresh = 0.2;
                    matlabbatch{9}.spm.meeg.preproc.artefact.append = true;
                    matlabbatch{9}.spm.meeg.preproc.artefact.methods.channels{1}.all = 'all';
                    matlabbatch{9}.spm.meeg.preproc.artefact.methods.fun.threshchan.threshold = 80;
                    matlabbatch{9}.spm.meeg.preproc.artefact.methods.fun.threshchan.excwin = 1000;
                    matlabbatch{9}.spm.meeg.preproc.artefact.prefix = 'a';
                    prefix = ['a' prefix]; % updates the prefix

                case 9 % averaging

                    matlabbatch{10}.spm.meeg.averaging.average.D = [datapath prefix spm_filename];
                    matlabbatch{10}.spm.meeg.averaging.average.userobust.robust.ks = 3;
                    matlabbatch{10}.spm.meeg.averaging.average.userobust.robust.bycondition = true;
                    matlabbatch{10}.spm.meeg.averaging.average.userobust.robust.savew = true;
                    matlabbatch{10}.spm.meeg.averaging.average.userobust.robust.removebad = false;
                    matlabbatch{10}.spm.meeg.averaging.average.plv = false;
                    matlabbatch{10}.spm.meeg.averaging.average.prefix = 'm';
                    prefix = ['m' prefix]; % updates the prefix

            end
        end
    end
end
