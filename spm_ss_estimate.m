function ss=spm_ss_estimate(ss)
% SPM_SS_ESTIMATE subject-specific ROI- or voxel- based model estimation
% 
% ss=spm_ss_estimate(ss)
% see SPM_SS_DESIGN

% spm_ss initialization
spm_ss init silent;

if nargin<1, 
    str='Select spm_ss*.mat analysis file';
    disp(str);
    Pdefault='';objname=findobj('tag','spm_ss');if numel(objname)==1,objdata=get(objname,'userdata');if isfield(objdata,'files_spm_ss'),Pdefault=objdata.files_spm_ss;end;end;
    P=spm_select(1,'^SPM_ss.*\.mat$',str,{Pdefault});
    if numel(objname)==1&&~isempty(P),objdata.files_spm_ss=P;set(objname,'userdata',objdata);end;
    load(P);
    ss.swd=fileparts(P);
end

ss=spm_ss_design(ss); % checks that information is complete
if ss.typen==1,  ss=spm_ss_estimate_voxel(ss);
else            ss=spm_ss_estimate_ROI(ss); end


