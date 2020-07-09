% example script.
% This script will perform GcSS fROI analyses, using the contrast 'S-N'
% as localizer contrast and estimating the contrasts of interest for the
% individual effects S, W, J, and N (the toolbox will automatically break 
% down the non-orthogonal contrasts by sessions and perform cross-validation 
% with the resulting session-specific estimates)
%

experiments=struct(...
    'name','test',...
    'pwd1','/users/alfnie/test17',... %'/mindhive/evlab/u/Shared/SUBJECTS',...
    'pwd2','model_langlocSN',...
    'data',{{'288_FED_20150613a_3T1_PL2017','325_FED_20150827b_3T2_PL2017','331_FED_20151111c_3T1_PL2017','340_FED_20151217b_3T2_PL2017','594_FED_20170423d_3T2_PL2017'}});

spmfiles={};
for nsub=1:length(experiments.data),
    spmfiles{nsub}=fullfile(experiments.pwd1,experiments.data{nsub},experiments.pwd2,'SPM.mat');
end

ss=struct(...
    'swd','/users/alfnie/test17/testlanglocSN',...                                   % output directory
    'files_spm',{spmfiles},...                              % first-level SPM.mat files
    'EffectOfInterest_contrasts',{{'S','N','S-N'}},...    % contrasts of interest
    'Localizer_contrasts',{{'S-N'}},...                     % localizer contrast (note: if these contrasts are not orthogonal the toolbox will automatically partition theses contrasts by sessions and perform cross-validation) 
    'Localizer_thr_type','FDR',...
    'Localizer_thr_p',.05,...
    'type','GcSS',...                                       % can be 'GcSS' (for automatically defined ROIs), 'mROI' (for manually defined ROIs), or 'voxel' (for voxel-based analyses)
    'smooth',6,...                                          % (FWHM mm)
    'ExplicitMasking','',...
    'model',1,...                                           % can be 1 (one-sample t-test), 2 (two-sample t-test), or 3 (multiple regression)
    'ask','none');                                       % can be 'none' (any missing information is assumed to take default values), 'missing' (any missing information will be asked to the user), 'all' (it will ask for confirmation on each parameter)
ss=spm_ss_design(ss);                                          % see help spm_ss_design for additional information
ss=spm_ss_estimate(ss);





