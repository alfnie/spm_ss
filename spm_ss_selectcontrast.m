function [ss,Ic]=spm_ss_selectcontrast(ss)
% SPM_SS_SELECTCONTRAST gui to select (or define new) contrasts
%
% [ss,Ic]=spm_ss_selectcontrast(ss);
%

posstr=1;
if nargin<1, 
    str='Select spm_ss*.mat analysis file';
    disp(str);
    Pdefault='';objname=findobj('tag','spm_ss');if numel(objname)==1,objdata=get(objname,'userdata');if isfield(objdata,'files_spm_ss'),Pdefault=objdata.files_spm_ss;end;end;
    P=spm_select(1,'^SPM_ss.*\.mat$',str,{Pdefault});
    if numel(objname)==1&&~isempty(P),objdata.files_spm_ss=P;set(objname,'userdata',objdata);end;
    load(P);
    ss.swd=fileparts(P);
end

if numel(ss.C)>=1,Cnames=cat(2,{ss.C(:).name});else Cnames={};end;Cnames=cat(2,Cnames,{'<define new contrast>'});
str='Select 2nd-level contrast';
Ic=listdlg('promptstring',str,'selectionmode','single','liststring',Cnames);

if Ic==numel(Cnames),
    % define new contrast
    if size(ss.X,2)>1,
        %figure('units','norm','position',[0.5,0.05,0.2,0.4],'color','w','menubar','none','numbertitle','off','name','design matrix'); imagesc(ss.X);colormap gray;set(gca,'units','norm','position',[.25,.4,.6,.5],'xtick',1:size(ss.X,2),'ytick',1:size(ss.X,1),'xticklabel',[]);text(1:size(ss.X,2),size(ss.X,1)+ones(1,size(ss.X,2)),ss.Xname,'rotation',90,'horizontalalignment','right','fontsize',10);ylabel('\fontsize{16}Subjects');title('\fontsize{16}Regressors');
        if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
        str='Regressors: ';for n1=1:numel(ss.Xname),str=[str,'#',num2str(n1),'=',ss.Xname{n1},'    '];end;disp(str);
        str='Between-subjects Contrast vector/matrix (number of contrasts x number of regressors)?';disp(str);
        newC.between=spm_input(str,posstr,'r','',[inf,size(ss.X,2)]);posstr='+1';
    else
        newC.between=1;
    end
    if numel(ss.EffectOfInterest_contrasts)>1 || size(ss.X,2)==1,
        if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
        str='Effects of interest: ';for n1=1:numel(ss.EffectOfInterest_contrasts),str=[str,'#',num2str(n1),'=',ss.EffectOfInterest_contrasts{n1},'    '];end;disp(str);
        str='Within-subjects Contrast vector/matrix (number of contrasts x number of effects of interest)?';disp(str);
        newC.within=spm_input(str,posstr,'r','',[inf,numel(ss.EffectOfInterest_contrasts)]);posstr='+1';
    else
        newC.within=1;
    end
    if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
    newC.name='';while isempty(newC.name), newC.name=spm_input('Contrast name?',posstr,'s'); end; posstr='+1';
    ss.C(Ic)=newC;
end

