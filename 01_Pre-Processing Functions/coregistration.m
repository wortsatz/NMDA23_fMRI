
function [matlabbatch] = coregistration(a)

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
matlabbatch.spm.spatial.coreg.estimate.ref    = cellstr(a(1,:));
matlabbatch.spm.spatial.coreg.estimate.source = cellstr(a);
matlabbatch.spm.spatial.coreg.estimate.other = {''};
matlabbatch.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

end
