% This batch script analyses the Auditory fMRI dataset available from the 
% SPM website:
%   http://www.fil.ion.ucl.ac.uk/spm/data/auditory/
% as described in the SPM manual:
%   http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#Chap:data:auditory
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: auditory_spm12_batch.m 8 2014-09-29 18:11:56Z guillaume $

% Directory containing the Auditory data


% =================== SAVE AND ACCESS THE DATA ======================

%% (TL) if you have github local directory set up then you don't need to run this first cell
%--------------------------------------------------------------------------
data_path = fileparts(mfilename('fullpath'));
if isempty(data_path), data_path = pwd; end %(TL) just navigate to github NMDA23 folder 
fprintf('%-40s:', 'Downloading Auditory dataset...');
urlwrite('http://www.fil.ion.ucl.ac.uk/spm/download/data/MoAEpilot/MoAEpilot.zip','MoAEpilot.zip');
unzip(fullfile(data_path,'MoAEpilot.zip'));
fprintf(' %30s\n', '...done');

 
%% (TL) add path to spm and cd into directory contraining Auditory data
whichComp=3;

if whichComp==1
    spmPath='/Users/ttli/Dropbox/spm12';
    data_path='/Users/ttli/Dropbox/Mac (2)/Documents/GitHub/NMDA23/fMRI_data';
elseif whichComp==2
    spmPath='/Users/kubrafatulla/Documents/spm12';
    data_path='/Users/kubrafatulla/Documents/fMRI/MoAEpilot/';
elseif whichComp==3 
    spmPath='/Users/lp1/Documents/spm12';
    data_path='/Users/lp1/Nextcloud/FU/NMDA_data';
else
    spmPath='/Users/Nastya/программы/spm12';
    data_path='/Users/Nastya/Documents/Freie Universitat/NMDA';
end
cd(data_path)
addpath(spmPath)

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

% create variables for the functional data (f) and the structural data (a)

f = spm_select('FPList', fullfile(data_path,'fM00223'), '^f.*\.img$'); % gives the full path of all data files 
a = spm_select('FPList', fullfile(data_path,'sM00223'), '^s.*\.img$');

clear matlabbatch
%--------------------------------------------------------------------------
% Realign
% initialize the parameters for the realignment
r = realignment(f);
matlabbatch{1} = r{1};
%--------------------------------------------------------------------------
% Coregister
%--------------------------------------------------------------------------
% initialize the parameters for the coregistration
c = coregistration(a);
matlabbatch{2} = c{2};
%--------------------------------------------------------------------------
% Segment
% initialize the parameters for the segmentation
s = segmentation(a);
matlabbatch{3} = s{3}; 
%--------------------------------------------------------------------------
% Normalise: Write
% initialize the parameters for the normalization
n = normalization(a, f);
matlabbatch{4} = n{4};
matlabbatch{5} = n{5};
%--------------------------------------------------------------------------
% Smooth
% initialize the parameters for the smoothing
sm = smoothing(f);
matlabbatch{6} = sm{6};
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
