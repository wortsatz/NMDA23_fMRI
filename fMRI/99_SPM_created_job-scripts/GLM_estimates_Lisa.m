% List of open inputs
nrun = X; % enter the number of runs here
jobfile = {'/Users/lp1/Documents/GitHub/NMDA23/99_SPM_created_job-scripts/GLM_estimates_Lisa_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
