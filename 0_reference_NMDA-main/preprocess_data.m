all_codes = [11 12 13 21 22 23 31 32 33 62];
all_labels = {'rhythm_std' 'rhythm_dev_word' 'rhythm_dev_syll' 'jitter_word_std' 'jitter_word_dev_word' 'jitter_word_dev_syll' 'jitter_syll_std' 'jitter_syll_dev_word' 'jitter_syll_dev_syll' 'rep'};

spmdir = 'C:\Users\45040\Documents\MATLAB\ana_test\spm12';
eeglabdir = 'C:\Users\45040\Documents\MATLAB\ana_test\eeglab2023.0';
noisedir = 'C:\Users\45040\Documents\MATLAB\ana_test\NoiseTools';

addpath(genpath(spmdir))

fiducial = {{[0.0 83.1 -43.0] [-74.7 -6.4 -40.9] [72.5 -8.6 -45.2]} ... % nas lpa rpa, mm
    {[0.0 94.1 -6.1] [-71.2 -7.7 -58.8] [72.2 -2.6 -58.8]} ...
    {[-7.3 83.0 -34.8] [-72.0 -25.2 -44.8] [72.9 -22.8 -44.8]} ...
    {[-0.0 83.7 -16.5] [-67.1 -9.8 -65.0] [73.2 -9.3 -60.9]} ...
    {[-2.2 81.9 -25.5] [-71.4 -14.9 -51.0] [66.0 -12.8 -65.9]} ...
    {[-1.1 80.2 -11.3] [-73.1 -20.6 -52.5] [66.8 -22.6 -55.6]} ...
    {[2.1 92.9 -32.3] [-79.3 -12.5 -47.9] [64.7 -17.7 -57.3]} ... 
    {[-2.6 86.7 -49.6] [-74.4 -24.2 -46.5] [71.3 -20.6 -56.8]} ...
    {[-0.0 82.6 -17.7] [-71.2 -39.7 -45.8] [77.4 -26.1 -38.5]} ...
    {[-0.0 74.2 -48.7] [-88.3 -27.6 -41.3] [74.4 -42.4 -54.0]} ...
    {[-3.2 78.7 -40.4] [-71.2 -29.8 -64.8] [79.8 -21.3 -55.3]} ...
    {[6.6 86.8 -24.6] [-70.9 -4.7 -53.9] [77.5 -10.3 -40.6]} ...
    {[0.0 90.1 -16.6] [-71.8 -30.1 -62.0] [71.8 -25.9 -62.0]} ...
    {[2.0 88.5 -9.3] [-75.9 -18.5 -51.4] [81.1 -16.5 -51.4]} ...
    {[2.1 77.3 -37.0] [-69.1 -24.8 -37.0] [69.1 -27.9 -37.0]} ...
    {[-2.1 85.0 -21.7] [-71.2 -20.7 -60.0] [73.2 -18.6 -51.7]} ...
    {[6.4 91.6 -8.5] [-80.1 -19.2 -39.3] [74.7 -25.5 -53.2]} ...
    {[0.0 80.9 -17.6] [-71.2 -14.0 -43.6] [65.6 -10.3 -56.5]} ...
    {[-0.0 81.9 -37.8] [-78.3 -23.6 -59.3] [78.3 -32.8 -53.2]} ...
    {[0.0 86.7 -12.7] [-71.8 -11.6 -59.2] [72.9 -3.1 -58.2]} ...
    {[-0.0 86.2 -37.1] [-82.0 -11.3 -47.4] [75.8 -19.5 -54.6]} ...
    {[-6.6 83.6 -14.1] [-66.7 -12.2 -65.6] [73.3 -8.4 -45.9]}};

datadir = 'C:\Users\45040\Documents\MATLAB\ana_test\continuous_data';
cd(datadir)

subs = dir('*_*');

for s = 1:length(subs)
% for s = [17] % s=17 i=1
    cd(subs(s).name)

    load StimVarsSpeechFFM.mat;

    temp=dir('*.ds');
% 
    for i = 1:length(temp)

        % convert
        clear S
        S.dataset = temp(i).name;
        S.mode = 'continuous';
        if i<10
            S.outfile = strcat('block_0',num2str(i));
        else
            S.outfile =   strcat('block_',num2str(i));
        end
        spm_eeg_convert(S);
        
        if s==8 && i==7 % missing header info
            D = spm_eeg_load(strcat(S.outfile,'.mat'));
            D2 = spm_eeg_load('block_06.mat');
            D = fiducials(D,D2.fiducials);
            D = sensors(D,'MEG',D2.sensors('MEG'));
            save(D);
        end

        % crop
        clear S
        if i<10
            S.D = strcat('block_0',num2str(i));
        else
            S.D = strcat('block_',num2str(i));
        end
        D = spm_eeg_load(S.D);
        D = chantype(D,indchannel(D,'EEG057'),'EOG');
        D = chantype(D,indchannel(D,'EEG058'),'EOG');
        D = chantype(D,indchannel(D,'EEG059'),'ECG');
        save(D);
        evs = D.events;
        S.timewin = [0 evs(max(find(strcmp({evs.type},'UPPT001_down')))).time+15]*1000;
        S.channels = {'MEG' 'Other' 'EOG' 'ECG'};
        D = spm_eeg_crop(S);

        if ~(s==8 && i==6) % missing header info - don't delete
            delete(strcat(S.D(1:end-3),'*'))
        end

        % high pass
        clear S
        S.D = D.fname;
        S.band = 'high';
        S.order = 4;
        S.freq = .1;
        D = spm_eeg_filter(S);
        delete(strcat(S.D(1:end-3),'*'))   

        % notch
        clear S
        S.D = D.fname;
        S.band = 'stop';
        S.freq = [48 52];
        D = spm_eeg_filter(S);
        delete(strcat(S.D(1:end-3),'*'))


        % downsample
        clear S
        S.D = D.fname;
        S.fsample_new = 300;
        D = spm_eeg_downsample(S);
        delete(strcat(S.D(1:end-3),'*'))

        % low pass
        clear S
        S.D = D.fname;
        S.band = 'low';
        S.freq = 90;
        D = spm_eeg_filter(S);
        delete(strcat(S.D(1:end-3),'*'))
        
        % remove movement artefacts % https://doi.org/10.1016/j.neuroimage.2012.11.047

        clear cc
        coil1=[D(indchannel(D,'HLC0011'),:); D(indchannel(D,'HLC0012'),:); D(indchannel(D,'HLC0013'),:)];
        coil2=[D(indchannel(D,'HLC0021'),:); D(indchannel(D,'HLC0022'),:); D(indchannel(D,'HLC0023'),:)];
        coil3=[D(indchannel(D,'HLC0031'),:); D(indchannel(D,'HLC0032'),:); D(indchannel(D,'HLC0033'),:)];
        N = size(coil1,2);
        xba = coil2(1,:) - coil1(1,:);
        yba = coil2(2,:) - coil1(2,:);
        zba = coil2(3,:) - coil1(3,:);
        xca = coil3(1,:) - coil1(1,:);
        yca = coil3(2,:) - coil1(2,:);
        zca = coil3(3,:) - coil1(3,:);
        balength = xba .* xba + yba .* yba + zba .* zba;
        calength = xca .* xca + yca .* yca + zca .* zca;
        xcrossbc = yba .* zca - yca .* zba;
        ycrossbc = zba .* xca - zca .* xba;
        zcrossbc = xba .* yca - xca .* yba;
        denominator = 0.5 ./ (xcrossbc .* xcrossbc + ycrossbc .* ycrossbc + zcrossbc .* zcrossbc);
        xcirca = ((balength .* yca - calength .* yba) .* zcrossbc - (balength .* zca - calength .* zba) .* ycrossbc) .* denominator;
        ycirca = ((balength .* zca - calength .* zba) .* xcrossbc - (balength .* xca - calength .* xba) .* zcrossbc) .* denominator;
        zcirca = ((balength .* xca - calength .* xba) .* ycrossbc - (balength .* yca - calength .* yba) .* xcrossbc) .* denominator;
        cc(1,:) = xcirca + coil1(1,:);
        cc(2,:) = ycirca + coil1(2,:);
        cc(3,:) = zcirca + coil1(3,:);
        v = [cc(1,:)', cc(2,:)', cc(3,:)'];
        vx = [zeros(1,N)', cc(2,:)', cc(3,:)']; % on the x-axis
        vy = [cc(1,:)', zeros(1,N)', cc(3,:)']; % on the y-axis
        vz = [cc(1,:)', cc(2,:)', zeros(1,N)']; % on the z-axis
        for j = 1:N
            thetax(j) = acos(dot(v(j,:),vx(j,:))/(norm(v(j,:))*norm(vx(j,:))));
            thetay(j) = acos(dot(v(j,:),vy(j,:))/(norm(v(j,:))*norm(vy(j,:))));
            thetaz(j) = acos(dot(v(j,:),vz(j,:))/(norm(v(j,:))*norm(vz(j,:))));
            cc(4,j) = (thetax(j) * (180/pi));
            cc(5,j) = (thetay(j) * (180/pi));
            cc(6,j) = (thetaz(j) * (180/pi));
        end   
        for j=1:size(cc,1)
            cc(j,find(cc(j,:)>nanmedian(cc(j,:))+3*nanstd(cc(j,:)) | cc(j,:)<nanmedian(cc(j,:))-3*nanstd(cc(j,:))))=nanmedian(cc(j,:));
        end
        cc_rel = [cc - repmat(cc(:,1),1,size(cc,2))]';
        for j=1:size(cc_rel,2)
            cc_rel(isnan(cc_rel(:,j)),j)=nanmedian(cc_rel(:,j));
        end
        corrected_data = D(:,:,:);
        choi=indchantype(D,{'MEG'});
        for j=1:length(choi)
            Y = D(choi(j),:)';
            [b,bint,r] = regress(Y,[ones(size(cc_rel,1),1) cc_rel]);
            corrected_data(choi(j),:) = r + b(1);
        end
        D(:,:)=corrected_data;
        save(D);

        % replace bad channels with 0 (if not, noise will be folded back
        % into data while removing artefacts; if channels removed
        % completely, errors during removing artefacts)

        std_over_time=std(D(indchantype(D,'MEG'),:),1,2);
        badchans_MEG = find(std_over_time>prctile(std_over_time,99));
        if s == 10
            badchans_MEG = unique([badchans_MEG; 121; 127; 128; 129]);
        end
        tempchans = indchantype(D,'MEG');
        badchans_MEG = tempchans(badchans_MEG);
        D(badchans_MEG,:,:)=0;
        save(D);

        structfile = dir(strcat('../../../structural_scans/',subs(s).name,'/s*.nii'));
        structfile = fullfile(pwd,'../../../structural_scans/',subs(s).name,structfile(1).name);

        clear matlabbatch
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.D={D.fname};
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.val=1;
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.comment='';
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshes.mri = {strcat(structfile,',1')};
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshres = 2;
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(1).fidname = 'nas';
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(1).specification.type = fiducial{s}{1};
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(2).specification.type = fiducial{s}{2};
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(3).specification.type = fiducial{s}{3};
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.useheadshape=1;
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.forward.eeg='EEG BEM';
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.forward.meg='Single Shell';
        spm_jobman('run',matlabbatch,cell(0,1));

        % eyeblink

        clear S
        S.D=D.fname;
        if s>=4
            if s == 19
                if i < 8
                    S.methods.channels='MLF13';
                else
                    S.methods.channels='MRF13';
                end
            else
                S.methods.channels='EEG058';
            end
        else
            S.methods.channels='EEG057';
        end
        S.methods.settings.chanind=indchannel(D,S.methods.channels);
        S.methods.fun='eyeblink';
        S.mode='mark';
        S.append=1;
        S.methods.settings.excwin=400;
        S.badchanthresh=.2;
        S.prefix='a';
        eyeblink_completed = 0;
        initial_threshold = 3;
        while eyeblink_completed == 0
            try
                S.methods.settings.threshold = initial_threshold;
                spm_eeg_artefact(S);
                eyeblink_completed = 1;
            catch            
                initial_threshold = initial_threshold + 1;
            end
        end
            

        clear S
        S.D=strcat('a',D.fname);
        S.bc=0;
        S.timewin=[-200 400];
        S.trialdef(1).conditionlabel='eyeblink';
        S.trialdef(1).eventtype='artefact_eyeblink';
        if s>=4
            if s == 19
                if i < 8
                    S.trialdef(1).eventvalue='MLF13';
                else
                    S.trialdef(1).eventvalue='MRF13';
                end
            else
                S.trialdef(1).eventvalue='EEG058';
            end
        else
            S.trialdef(1).eventvalue='EEG057';
        end
        S.reviewtrials=0;
        S.save=1;
        spm_eeg_epochs(S);

        clear S
        S.D=strcat('ea',D.fname);
        S.method='SVD';
        S.timewin=[-200 400];
        S.ncomp=2;
        spm_eeg_spatial_confounds(S);

        clear S
        S.D=strcat('a',D.fname);
        S.method='SPMEEG';
        S.conffile=strcat('ea',D.fname);
        spm_eeg_spatial_confounds(S);

        clear S
        S.D=strcat('a',D.fname);
        S.correction='Berg';
        S.save=1;
        D=spm_eeg_correct_sensor_data(S);
        delete(strcat(S.D(1:end-3),'*'))

        clear matlabbatch
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.D={D.fname};
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.val=1;
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.comment='';
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshes.template=1;
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshres=2;
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshes.mri = {strcat(structfile,',1')};
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshres = 2;
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(1).fidname = 'nas';
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(1).specification.type = fiducial{s}{1};
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(2).specification.type = fiducial{s}{2};
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(3).specification.type = fiducial{s}{3};
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.useheadshape=1;
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.forward.eeg='EEG BEM';
        matlabbatch{1}.spm.meeg.source.headmodelhelmet.forward.meg='Single Shell';
        spm_jobman('run',matlabbatch,cell(0,1));

        % heartbeat

        if ~(s==17 && i==6) % can't get it to work for this one

    
            addpath(genpath(eeglabdir))
    
            clear S
            S.D=D.fname;
            S.methods.channels='EEG059';
            S.methods.settings.chanind=indchannel(D,S.methods.channels);
            S.methods.settings.threshold=3;
            S.methods.fun='heartbeat';
            S.mode='mark';
            S.append=1;
            S.methods.settings.excwin=200;
            S.badchanthresh=.2;
            S.prefix='a';
            spm_eeg_artefact(S);
            
            clear S
            S.D=strcat('a',D.fname);
            S.bc=0;
            S.timewin=[-100 200];
            S.trialdef(1).conditionlabel='heartbeat';
            S.trialdef(1).eventtype='artefact_heartbeat';
            S.trialdef(1).eventvalue='EEG059';
            S.reviewtrials=0;
            S.save=1;
            spm_eeg_epochs(S);
    
            clear S
            S.D=strcat('ea',D.fname);
            S.method='SVD';
            S.timewin=[-100 200];
            S.ncomp=2;
            spm_eeg_spatial_confounds(S);
    
            clear S
            S.D=strcat('a',D.fname);
            S.method='SPMEEG';
            S.conffile=strcat('ea',D.fname);
            spm_eeg_spatial_confounds(S);
    
            clear S
            S.D=strcat('a',D.fname);
            S.correction='Berg';
            S.save=1;
            D=spm_eeg_correct_sensor_data(S);
            delete(strcat(S.D(1:end-3),'*'))
    
            rmpath(genpath(eeglabdir))

        else

            D2 = clone(D,strcat('Ta',D.fname));
            D2(:,:) = D(:,:);
            save(D2)

        end


        % update stimulus info

        evs = D.events;
        triggers = find(strcmp({evs.type},'UPPT001_down'));
        evstypes = cell2mat({evs(find(strcmp({evs.type},'UPPT001_down'))).value});
        which_evs = find(evstypes==60);
        triggers = triggers(which_evs);

        when_cnd = permute(when_conds((i-1)*20+1:i*20,:,:),[3 2 1]);
        when_cnd = when_cnd(:);
        what_cnd = permute(what_conds((i-1)*20+1:i*20,:,:),[3 2 1]);
        what_cnd = what_cnd(:);
        dev_all = permute(odd_times((i-1)*20+1:i*20,:,:),[3 2 1]);
        dev_all = dev_all(:);

        std_conds = when_cnd(find(dev_all==-1))*10+what_cnd(find(dev_all==-1))*1;
        dev_conds = when_cnd(find(dev_all==1))*10+what_cnd(find(dev_all==1))*1;
        all_conds = [std_conds dev_conds]';
        all_conds = all_conds(:);
        
        if length(all_conds) ~= length(triggers) % missing triggers
            all_conds = all_conds(end-length(triggers)+1:end);
        end
        
        for j=1:length(triggers)
            evs(triggers(j)).value = all_conds(j);
        end
        D = events(D,1,evs);
        save(D);

        % 30: when jittered syllables, predictable words
        % 20: when jittered words, predictable syllables 
        % 10: when rhythmic
        % 3: what deviant syllable
        % 2: what deviant word
        % 1: what standard

        % epoch

        trialcodes = unique(all_conds);

        clear S
        S.D = D.fname;
        S.bc = 1;
        S.timewin = [-50 250];
        for j=1:length(trialcodes)
            S.trialdef(j).conditionlabel = all_labels{find(all_codes==trialcodes(j))};
            S.trialdef(j).eventtype = 'UPPT001_down';
            S.trialdef(j).eventvalue = trialcodes(j);
        end
        D = spm_eeg_epochs(S);    

    end

    temp=dir('eTaTa*.mat');
    clear S
    for i=1:length(temp)
        S.D{i}=spm_eeg_load(temp(i).name);
    end
    S.recode = 'same';
    S.prefix = 'c';
    spm_eeg_merge(S)

    % denoise

    addpath(genpath(noisedir))
    temp=dir('ceTaTa*.mat');
    D = spm_eeg_load(temp(1).name);
    
%     % detrend MEG
%     xx = permute(D(indchantype(D,{'MEG'}),:,:),[2 1 3]);
%     ORDER=1;
%     [x,w]=nt_detrend(xx,ORDER);
%     x=permute(x,[2 1 3]);
%     D(indchantype(D,{'MEG'}),:,:)=x;
    % DSS
    xx = permute(D(indchantype(D,{'MEG'}),:,:),[2 1 3]);
    badtrls = [find(mean(std(xx,[],1),2)>median(mean(std(xx,[],1),2))+3*std(mean(std(xx,[],1),2))) find(mean(std(xx,[],1),2)<median(mean(std(xx,[],1),2))-3*std(mean(std(xx,[],1),2)))];
%     badtrls = [];
    % baseline correct
    xx = xx - repmat(mean(xx(min(find(D.time>=-.025)):max(find(D.time<.025)),:,:),1),[size(D,2) 1 1]);
    clear data
    for j=1:length(D.condlist)
        data{j}=xx(:,:,setdiff(indtrial(D,D.condlist{j}),badtrls));
    end
    nchans = size(xx,2);
    c0=zeros(nchans); c1=zeros(nchans);
    for iCondition=1:length(D.condlist)
        c0=c0+nt_cov(data{iCondition});
        c1=c1+nt_cov(mean(data{iCondition},3));
    end
    [todss,pwr0,pwr1]=nt_dss0(c0,c1);
    fromdss=pinv(todss);
    NKEEP = 8;
    for iCondition=1:length(D.condlist)
        z{iCondition}=nt_mmat(data{iCondition},todss(:,1:NKEEP)*fromdss(1:NKEEP,:));
    end
    for j=1:length(D.condlist)
        xx(:,:,setdiff(indtrial(D,D.condlist{j}),badtrls))=z{j};
    end
    xx = permute(xx,[2 1 3]);
    D(indchantype(D,{'MEG' 'EEG'}),:,:)=xx;
    D=badtrials(D,badtrls,ones(1,length(badtrls)));
    save(D);

    rmpath(genpath(noisedir))
% 
%     temp=dir('ceTaTa*.mat');
%     clear S
%     S.D=temp(1).name;
%     S.robust.savew=0;
%     S.robust.bycondition=1;
%     S.robust.ks=3;
%     S.prefix='m';
%     D=spm_eeg_average(S);
% 

     
%     clear S
%     S.D = D.fname;
%     S.band = 'low';
%     S.freq = [48];
%     D = spm_eeg_filter(S);
% 
%     clear S
%     S.D = D.fname;
%     S.c = zeros(6,9);
%     S.c(1,indtrial(D,'rhythm_dev_syll'))=1;
%     S.c(1,indtrial(D,'rhythm_std'))=-1;
%     S.c(2,indtrial(D,'rhythm_dev_word'))=1;
%     S.c(2,indtrial(D,'rhythm_std'))=-1;
%     S.c(3,indtrial(D,'jitter_syll_dev_syll'))=1;
%     S.c(3,indtrial(D,'jitter_syll_std'))=-1;
%     S.c(4,indtrial(D,'jitter_syll_dev_word'))=1;
%     S.c(4,indtrial(D,'jitter_syll_std'))=-1;
%     S.c(5,indtrial(D,'jitter_word_dev_syll'))=1;
%     S.c(5,indtrial(D,'jitter_word_std'))=-1;
%     S.c(6,indtrial(D,'jitter_word_dev_word'))=1;
%     S.c(6,indtrial(D,'jitter_word_std'))=-1;
%     S.label = {'rhythm_mismatch_syll' 'rhythm_mismatch_word' 'jitter_syll_mismatch_syll' 'jitter_syll_mismatch_word'  'jitter_word_mismatch_syll' 'jitter_word_mismatch_word'};
%     D = spm_eeg_contrast(S);
%     
%     datafile = dir('fmce*.mat');
%     datafile = datafile(1).name;
% %     D = spm_eeg_load(datafile);
% %     fidold = fiducials(D);
% %     fidold.fid.pnt = [105 0 0; 0 75 -75; 0 0 0]';
% %     D = fiducials(D,fidold);
% %     save(D);
%     structfile = dir(strcat('../../structural_scans/',subs(s).name,'/rs*.nii'));
%     structfile = fullfile(pwd,'../../structural_scans/',subs(s).name,structfile(1).name);
%     
%     clear matlabbatch
%     matlabbatch{1}.spm.meeg.source.headmodelhelmet.D = {fullfile(pwd,datafile)};
%     matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshes.template=1;
%     matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshres=2;
% %     matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshes.mri = {strcat(structfile,',1')};
% %     matlabbatch{1}.spm.meeg.source.headmodelhelmet.meshing.meshres = 2;
%     matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(1).fidname = 'nas';
% %     matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(1).specification.type = fiducial{s}{1};
%     matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(1).specification.select = 'nas';
%     matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
% %     matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(2).specification.type = fiducial{s}{2};
%     matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(2).specification.select = 'lpa';
%     matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
% %     matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(3).specification.type = fiducial{s}{3};
%     matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.fiducial(3).specification.select = 'rpa';
%     matlabbatch{1}.spm.meeg.source.headmodelhelmet.coregistration.coregspecify.useheadshape=1;
%     matlabbatch{1}.spm.meeg.source.headmodelhelmet.forward.eeg = 'EEG BEM';
%     matlabbatch{1}.spm.meeg.source.headmodelhelmet.forward.meg = 'Single Shell';
%     matlabbatch{2}.spm.meeg.source.invert.D(1) = cfg_dep('MEG helmet head model specification: M/EEG dataset(s) with a forward model', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','D'));    
%     matlabbatch{2}.spm.meeg.source.invert.val = 1;
%     matlabbatch{2}.spm.meeg.source.invert.whatconditions.all = 1;
% %     matlabbatch{2}.spm.meeg.source.invert.whatconditions.condlabel = {
% %                                                                   'dev'
% %                                                                   'std'
% %                                                                   }';    
%     matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.invtype = 'GS';
%     matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.woi = [0 250];
%     matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.foi = [];
%     matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.hanning = 1;
%     matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.priors.priorsmask = {''};
%     matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.priors.space = 1;
%     matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.restrict.locs = zeros(0, 3);
%     matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.restrict.radius = 32;
%     matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.restrict.mask = {''};
%     matlabbatch{2}.spm.meeg.source.invert.modality = {'MEG'};
%     spm_jobman('run',matlabbatch,cell(0,1));
% %     matlabbatch{1}.spm.meeg.source.headmodelhelmet.D = {fullfile(pwd,strcat('w',datafile))};
% %     spm_jobman('run',matlabbatch,cell(0,1));
% 
% %    
    cd ..
end
% 
% L = size(xAudio,2);
% for i=1:240
%     i
%     Y = fft(envelope(xAudio(i,:),30,'rms'));
%     P2 = abs(Y/L);
%     P1 = P2(1:L/2+1);
%     P1(2:end-1) = 2*P1(2:end-1);
%     spctrm(i,:) = P1;
% end
% meanspctrm(1,:)=mean(spctrm(find(when_conds==1),:),1);
% meanspctrm(2,:)=mean(spctrm(find(when_conds==2),:),1);
% meanspctrm(3,:)=mean(spctrm(find(when_conds==3),:),1);
% 
% Fs = 24414;
% f = Fs*(0:(L/2))/L;
% figure;
% plot(f,meanspctrm) 
% xlim([0.5 5])
% hold on
% line([4 4],[0 .25],'linestyle','--','color','k')
% line([2 2],[0 .25],'linestyle','--','color','k')
% 
% to_plot = xAudio(81,:);
% to_plot = envelope(to_plot,300,'rms');
% t_ax = (1:length(to_plot))/Fs;
% figure;plot(t_ax,to_plot)
% hold on
% for i=1:56
%     line([i*.25 i*.25],[0 1.2],'color','r','linestyle','-')
% end