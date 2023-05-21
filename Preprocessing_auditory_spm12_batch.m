% This batch script analyses the Auditory fMRI dataset available from the 
% SPM website:
%   http://www.fil.ion.ucl.ac.uk/spm/data/auditory/
% as described in the SPM manual:
%   http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#Chap:data:auditory
%__________________________________________________________________________

% =================== SAVE AND ACCESS THE DATA ======================

%% (TL) if you have github local directory set up then you don't need to run this first cell
%--------------------------------------------------------------------------
data_path = fileparts(mfilename('fullpath'));
if isempty(data_path), data_path = pwd; end %(TL) just navigate to github NMDA23 folder 
fprintf('%-40s:', 'Downloading Auditory dataset...');
urlwrite('http://www.fil.ion.ucl.ac.uk/spm/download/data/MoAEpilot/MoAEpilot.zip','MoAEpilot.zip');
unzip(fullfile(data_path,'MoAEpilot.zip'));
fprintf(' %30s\n', '...done');

<<<<<<< Updated upstream
 
%% (TL) add path to spm and cd into directory contraining Auditory data
whichComp=2;
=======
%%
%(TL) add path to spm and cd into directory contraining Auditory data
whichComp=2;

if whichComp==1
    spmPath='/Users/ttli/Dropbox/spm12';
    data_path='/Users/ttli/Dropbox/Mac (2)/Documents/GitHub/NMDA23/fMRI_data';
elseif whichComp==2
    spmPath='/Users/kubrafatulla/Documents/spm12';
    data_path='/Users/kubrafatulla/Documents/MoAEpilot';
elseif whichComp==3 
    spmPath='/Users/lp1/Documents/spm12';
    data_path='/Users/lp1/Nextcloud/FU/NMDA_data';
else
    spmPath='/Users/USERNAME/WHERE/spm12';
    data_path='ADD';
end
cd(data_path)
addpath(spmPath)

%%

% Initialise SPM
%--------------------------------------------------------------------------
spm('Defaults','fMRI');
spm_jobman('initcfg');
spm_get_defaults('cmdline',true);

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

% Realign
% data

%--------------------------------------------------------------------------

% initialize the parameters for the realignment

matlabbatch{1}.spm.spatial.realign.estwrite.data = {cellstr(f)};
mmatlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9; % Highest quality (1) gives most precise results, whereas lower qualities gives faster realignment.
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4; % The separation (in mm) between the points sampled in the reference image. Smaller sampling distances gives more accurate results, but will be slower.
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5; % The FWHM of the Gaussian smoothing kernel (mm) applied to the images before estimating the realignment parameters.
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1; % Register to first: Images are registered to the first image in the series. Register to mean: A two pass procedure is used in order to register the images to the mean of the images after the first realignment
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4; matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';


% Coregister
%--------------------------------------------------------------------------

% initialize the parameters for the coregistration

matlabbatch{2}.spm.spatial.coreg.estimate.ref    = cellstr(spm_file(f(1,:),'prefix','mean'));
matlabbatch{2}.spm.spatial.coreg.estimate.source = cellstr(a);
matlabbatch{2}.spm.spatial.coreg.estimate.other = {''};
matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{2}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
spm_jobman('run',matlabbatch);
% Segment
%--------------------------------------------------------------------------

% initialize the parameters for the segmentation

matlabbatch{3}.spm.spatial.preproc.channel.vols  = cellstr(a);
matlabbatch{3}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{3}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{3}.spm.spatial.preproc.channel.write = [0 1];
matlabbatch{3}.spm.spatial.preproc.tissue(1).tpm =  {'/Users/kubrafatulla/Documents/spm12/tpm/TPM.nii,1'};
matlabbatch{3}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{3}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{3}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(2).tpm =  {'/Users/kubrafatulla/Documents/spm12/tpm/TPM.nii,2'};
matlabbatch{3}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{3}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{3}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(3).tpm =  {'/Users/kubrafatulla/Documents/spm12/tpm/TPM.nii,3'};
matlabbatch{3}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{3}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{3}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(4).tpm =  {'/Users/kubrafatulla/Documents/spm12/tpm/TPM.nii,4'};
matlabbatch{3}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{3}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{3}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(5).tpm =  {'/Users/kubrafatulla/Documents/spm12/tpm/TPM.nii,5'};
matlabbatch{3}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{3}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{3}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(6).tpm =  {'/Users/kubrafatulla/Documents/spm12/tpm/TPM.nii,6'};
matlabbatch{3}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{3}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{3}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{3}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{3}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{3}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{3}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{3}.spm.spatial.preproc.warp.write = [0 1];
matlabbatch{3}.spm.spatial.preproc.warp.vox = NaN;
matlabbatch{3}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                                           NaN NaN NaN];


% Normalise: Write
%--------------------------------------------------------------------------

% initialize the parameters for the normalization

% This is the normalization of the functional image 
matlabbatch{4}.spm.spatial.normalise.write.subj.def      = cellstr(spm_file(a,'prefix','y_','ext','nii'));
matlabbatch{4}.spm.spatial.normalise.write.subj.resample = cellstr(f);
matlabbatch{4}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
matlabbatch{4}.spm.spatial.normalise.write.woptions.vox = [3 3 3];
matlabbatch{4}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{4}.spm.spatial.normalise.write.woptions.prefix = 'w';
spm_jobman('run',matlabbatch);

% This is for the Normalization of the structural image 
matlabbatch{5}.spm.spatial.normalise.write.subj.def      = cellstr(spm_file(a,'prefix','y_','ext','nii'));
matlabbatch{5}.spm.spatial.normalise.write.subj.resample = cellstr(spm_file(a,'prefix','m','ext','nii'));
matlabbatch{5}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
matlabbatch{5}.spm.spatial.normalise.write.woptions.vox = [1 1 3];
matlabbatch{5}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{5}.spm.spatial.normalise.write.woptions.prefix = 'w';
spm_jobman('run',matlabbatch);
% Smooth
%--------------------------------------------------------------------------

% initialize the parameters for the smoothing

matlabbatch{6}.spm.spatial.smooth.data = cellstr(spm_file(f,'prefix','w'));
matlabbatch{6}.spm.spatial.smooth.fwhm = [6 6 6];
matlabbatch{6}.spm.spatial.smooth.dtype = 0;
matlabbatch{6}.spm.spatial.smooth.im = 0;
matlabbatch{6}.spm.spatial.smooth.prefix = 's';

spm_jobman('run',matlabbatch);
