function [matlabbatch] = smoothing(f)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
matlabbatch.spm.spatial.smooth.data = cellstr(spm_file(f,'prefix','w'));
matlabbatch.spm.spatial.smooth.fwhm = [6 6 6];
matlabbatch.spm.spatial.smooth.dtype = 0;
matlabbatch.spm.spatial.smooth.im = 0;
matlabbatch.spm.spatial.smooth.prefix = 's';

end

