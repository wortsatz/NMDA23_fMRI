
function [matlabbatch] = realignment(f)

matlabbatch.spm.spatial.realign.estwrite.data = {cellstr(f)};
matlabbatch.spm.spatial.realign.estwrite.eoptions.quality = 0.9; 
matlabbatch.spm.spatial.realign.estwrite.eoptions.sep = 4; 
matlabbatch.spm.spatial.realign.estwrite.eoptions.fwhm = 5; 
matlabbatch.spm.spatial.realign.estwrite.eoptions.rtm = 1; 
matlabbatch.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch.spm.spatial.realign.estwrite.roptions.which = [0 1];
matlabbatch.spm.spatial.realign.estwrite.roptions.interp = 4; 
matlabbatch.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch.spm.spatial.realign.estwrite.roptions.prefix = 'r';

end 