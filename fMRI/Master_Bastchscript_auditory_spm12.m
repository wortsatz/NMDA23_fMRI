% Master Batch Script

% =================== SAVE AND ACCESS THE DATA ======================

%% (TL) if you have github local directory set up 
% then you don't need to run this first cell
% **** SKIP *****
%--------------------------------------------------------------------------
data_path = fileparts(mfilename('fullpath'));
if isempty(data_path), data_path = pwd; end %(TL) just navigate to github NMDA23 folder 
fprintf('%-40s:', 'Downloading Auditory dataset...');
urlwrite('http://www.fil.ion.ucl.ac.uk/spm/download/data/MoAEpilot/MoAEpilot.zip','MoAEpilot.zip');
unzip(fullfile(data_path,'MoAEpilot.zip'));
fprintf(' %30s\n', '...done');

 
%% (TL) add path to spm and cd into directory contraining Auditory data
clc
clear
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


%%

% Initialise SPM
%--------------------------------------------------------------------------
spm('Defaults','fMRI');
spm_jobman('initcfg');
% spm_get_defaults('cmdline',true);

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PREAMBLE: DUMMY SCANS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% **** SKIP *****
% don't need to run this part either 
% this excludes the first session a priori and puts the excluded images
% into the dummy folder

f = spm_select('FPList', fullfile(data_path,'fM00223'), '^f.*\.img$') ;

clear matlabbatch

matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = cellstr(data_path);
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'dummy';

matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_move.files = cellstr(f(1:12,:));
matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_move.action.moveto = cellstr(fullfile(data_path,'dummy'));

spm_jobman('run',matlabbatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SPATIAL PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%######################### auditory tutorial ########################
% f = spm_select('FPList', fullfile(data_path,'fM00223'), '^f.*\.img$');
% a = spm_select('FPList', fullfile(data_path,'sM00223'), '^s.*\.img$');
% 
% clear matlabbatch
% 
% % Realign
% %--------------------------------------------------------------------------
% matlabbatch{1}.spm.spatial.realign.estwrite.data = {cellstr(f)};
% matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1];
% 
% % Coregister
% %--------------------------------------------------------------------------
% matlabbatch{2}.spm.spatial.coreg.estimate.ref    = cellstr(spm_file(f(1,:),'prefix','mean'));
% matlabbatch{2}.spm.spatial.coreg.estimate.source = cellstr(a);
% 
% % Segment
% %--------------------------------------------------------------------------
% matlabbatch{3}.spm.spatial.preproc.channel.vols  = cellstr(a);
% matlabbatch{3}.spm.spatial.preproc.channel.write = [0 1];
% matlabbatch{3}.spm.spatial.preproc.warp.write    = [0 1];
% 
% % Normalise: Write
% %--------------------------------------------------------------------------
% matlabbatch{4}.spm.spatial.normalise.write.subj.def      = cellstr(spm_file(a,'prefix','y_','ext','nii'));
% matlabbatch{4}.spm.spatial.normalise.write.subj.resample = cellstr(f);
% matlabbatch{4}.spm.spatial.normalise.write.woptions.vox  = [3 3 3];
% 
% matlabbatch{5}.spm.spatial.normalise.write.subj.def      = cellstr(spm_file(a,'prefix','y_','ext','nii'));
% matlabbatch{5}.spm.spatial.normalise.write.subj.resample = cellstr(spm_file(a,'prefix','m','ext','nii'));
% matlabbatch{5}.spm.spatial.normalise.write.woptions.vox  = [1 1 3];
% 
% % Smooth
% %--------------------------------------------------------------------------
% matlabbatch{6}.spm.spatial.smooth.data = cellstr(spm_file(f,'prefix','w'));
% matlabbatch{6}.spm.spatial.smooth.fwhm = [6 6 6];
% 
% spm_jobman('run',matlabbatch);
%##################################################  


% create variables for the functional data (f) and the structural data (a)

f = spm_select('FPList', fullfile(data_path,'fM00223'), '^f.*\.img$'); % gives the full path of all data files 
a = spm_select('FPList', fullfile(data_path,'sM00223'), '^s.*\.img$');

clear matlabbatch
%--------------------------------------------------------------------------
% Realign
% initialize the parameters for the realignment
r = realignment(f);
matlabbatch{1} = r;
%--------------------------------------------------------------------------
% Coregister
%--------------------------------------------------------------------------
% initialize the parameters for the coregistration
c = coregistration(a);
matlabbatch{2} = c;
%--------------------------------------------------------------------------
% Segment
% initialize the parameters for the segmentation
s = segmentation(spmPath,a);
matlabbatch{3} = s; 
%--------------------------------------------------------------------------
% Normalise: Write
% initialize the parameters for the normalization
[F,S] = normalization(a, f);
matlabbatch{4} = F;
matlabbatch{5} = S;
%--------------------------------------------------------------------------
% Smooth
% initialize the parameters for the smoothing
sm = smoothing(f);
matlabbatch{6} = sm;
%--------------------------------------------------------------------------
%%
 %Job script 
spm_jobman('run',matlabbatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GLM SPECIFICATION, ESTIMATION, INFERENCE, RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = spm_select('FPList', fullfile(data_path,'fM00223'), '^swf.*\.img$');

clear matlabbatch

% Output Directory
%--------------------------------------------------------------------------
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = cellstr(data_path);
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'GLM';

% Model Specification
%--------------------------------------------------------------------------
matlabbatch{2}.spm.stats.fmri_spec.dir = cellstr(fullfile(data_path,'GLM'));
matlabbatch{2}.spm.stats.fmri_spec.timing.units = 'scans';
matlabbatch{2}.spm.stats.fmri_spec.timing.RT = 7;
matlabbatch{2}.spm.stats.fmri_spec.sess.scans = cellstr(f);
matlabbatch{2}.spm.stats.fmri_spec.sess.cond.name = 'active';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond.onset = 6:12:84;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond.duration = 6;

% Model Estimation
%--------------------------------------------------------------------------
matlabbatch{3}.spm.stats.fmri_est.spmmat = cellstr(fullfile(data_path,'GLM','SPM.mat'));

% Contrasts
%--------------------------------------------------------------------------
matlabbatch{4}.spm.stats.con.spmmat = cellstr(fullfile(data_path,'GLM','SPM.mat'));
matlabbatch{4}.spm.stats.con.consess{1}.tcon.name = 'Listening > Rest';
matlabbatch{4}.spm.stats.con.consess{1}.tcon.weights = [1 0];
matlabbatch{4}.spm.stats.con.consess{2}.tcon.name = 'Rest > Listening';
matlabbatch{4}.spm.stats.con.consess{2}.tcon.weights = [-1 0];

% Inference Results
%--------------------------------------------------------------------------
matlabbatch{5}.spm.stats.results.spmmat = cellstr(fullfile(data_path,'GLM','SPM.mat'));
matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{5}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{5}.spm.stats.results.conspec.extent = 0;
matlabbatch{5}.spm.stats.results.print = false;

% Rendering
%--------------------------------------------------------------------------
matlabbatch{6}.spm.util.render.display.rendfile = {fullfile(spm('Dir'),'canonical','cortex_20484.surf.gii')};
matlabbatch{6}.spm.util.render.display.conspec.spmmat = cellstr(fullfile(data_path,'GLM','SPM.mat'));
matlabbatch{6}.spm.util.render.display.conspec.contrasts = 1;
matlabbatch{6}.spm.util.render.display.conspec.threshdesc = 'FWE';
matlabbatch{6}.spm.util.render.display.conspec.thresh = 0.05;
matlabbatch{6}.spm.util.render.display.conspec.extent = 0;

spm_jobman('run',matlabbatch);
