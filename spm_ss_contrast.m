function ss=spm_ss_contrast(ss,varargin)
% SPM_SS_CONTRAST evaluates contrast in subject-specific ROI- or voxel- based analyses
% 
% ss=spm_ss_contrast(ss,Ic)
% see SPM_SS_DESIGN

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
if ~isfield(ss,'estimate'), % checks that model has been estimated
    str={'This model has not been estimated. Estimate now?',...
        ['pwd = ',ss.swd]};
    if spm_input(str,1,'bd','stop|continue',[1,0],0),return;end
    ss=spm_ss_estimate(ss);
end
if ss.typen==1,  ss=spm_ss_contrast_voxel(ss,varargin{:});
else            ss=spm_ss_contrast_ROI(ss,varargin{:}); end


