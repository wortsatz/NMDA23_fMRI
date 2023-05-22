
function [matlabbatch] = realignment(f)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
matlabbatch{1}.spm.spatial.realign.estwrite.data = {cellstr(f)};
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9; 
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4; 
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5; 
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1; 
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4; 
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

end 