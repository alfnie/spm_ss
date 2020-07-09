% example script.
% This script will perform GcSS fROI analyses, using the contrast 'S-N'
% both as localizer contrast and contrast of interest (the toolbox will
% automatically break down these contrasts by sessions and perform
% cross-validation with the resulting session-specific estimates)
%

experiments=struct(...
    'name','SWJN_SWJN_and_SWJNV2',...% SWJN subjects, SWJNV2 subjects
    'pwd1','/groups/swjn/data/',...
    'pwd2','firstlevel_SWJNno',...
    'data',{{'SWJN_01','SWJN_05', 'SWJN_07', 'SWJN_09', 'SWJN_10', 'SWJN_11', 'SWJN_12','SWJN_13', 'SWJN_101', 'SWJN_102', 'SWJN_104', 'SWJN_105',...
    'SWJNV2_01', 'SWJNV2_02', 'SWJNV2_03', 'SWJNV2_04', 'SWJNV2_05', 'SWJNV2_06','SWJNV2_07', 'SWJNV2_08', 'SWJNV2_09', 'SWJNV2_10', 'SWJNV2_12', 'SWJNV2_13','SWJNV2_15'}});

spmfiles={};
for nsub=1:length(experiments.data),
    spmfiles{nsub}=fullfile(experiments.pwd1,experiments.data{nsub},experiments.pwd2,'SPM.mat');
end

ss=struct(...
    'swd',fullfile(pwd,'test'),...                          % output directory
    'files_spm',{spmfiles},...                              % first-level SPM.mat files
    'EffectOfInterest_contrasts',{{'S-N'}},...              % contrast of interest
    'Localizer_contrasts',{{'S-N'}},...                     % localizer contrast (note: if these contrasts are not orthogonal the toolbox will automatically partition theses contrasts by sessions and perform cross-validation) 
    'Localizer_thr_type','FDR',...
    'Localizer_thr_p',.05,...
    'type','GcSS',...                                       % can be 'GcSS' (for automatically defined ROIs), 'mROI' (for manually defined ROIs), or 'voxel' (for voxel-based analyses)
    'smooth',6,...                                          % (FWHM mm)
    'model',1,...                                           % can be 1 (one-sample t-test), 2 (two-sample t-test), or 3 (multiple regression)
    'ask','none');                                          % can be 'none' (any missing information is assumed to take default values), 'missing' (any missing information will be asked to the user), 'all' (it will ask for confirmation on each parameter)
ss=spm_ss_design(ss);                                          % see help spm_ss_design for additional information
ss=spm_ss_estimate(ss);





