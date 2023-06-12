function [matlabbatchF,matlabbatchS] = normalization(a, f)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% This is the normalization of the functional image 
matlabbatchF.spm.spatial.normalise.write.subj.def      = cellstr(spm_file(a,'prefix','y_','ext','nii'));
matlabbatchF.spm.spatial.normalise.write.subj.resample = cellstr(f);
matlabbatchF.spm.spatial.normalise.write.woptions.vox  = [3 3 3];

matlabbatchF.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
matlabbatchF.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatchF.spm.spatial.normalise.write.woptions.prefix = 'w';

% This is for the Normalization of the structural image 
matlabbatchS.spm.spatial.normalise.write.subj.def      = cellstr(spm_file(a,'prefix','y_','ext','nii'));
matlabbatchS.spm.spatial.normalise.write.subj.resample = cellstr(spm_file(a,'prefix','m','ext','nii'));
matlabbatchS.spm.spatial.normalise.write.woptions.vox  = [1 1 3];


end

