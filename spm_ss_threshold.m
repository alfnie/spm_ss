function [thr_p,thr_type,filenames]=spm_ss_threshold(option,spm_data,Ic,Ec,options,optionROIs,overwrite)
% SPM_SS_THRESHOLD
% creates localizer masks
% estimates localizer contrast optimal first-level FDR-corrected threshold across multiple subjects
%

persistent OptimalThreshold Threshold_type Threshold_p Conjunction_type AnyAutomatic SPMfiles ContrastIndexes ExpIndexes MaskFilenames

if nargin<7, overwrite=0; end
if nargin<6, optionROIs=''; end
if nargin<5, options=''; end
if nargin<1, %gui option
    str='Select first-level SPM.mat(s) (one per subject)';
    disp(str);
    Pdefault={''};objname=findobj('tag','spm_ss');if numel(objname)==1,objdata=get(objname,'userdata');if isfield(objdata,'files_spm'),Pdefault=objdata.files_spm;end;end;
    P=cellstr(spm_select(inf,'^SPM\.mat$',str,Pdefault));
    if numel(objname)==1&&~isempty(P),objdata.files_spm=P;set(objname,'userdata',objdata);end;
    
    spm_ss_threshold('begin');
    for np=1:numel(P),
        disp(['Loading subject #',num2str(np),'...']);
        load(P{np},'SPM');
        SPM.swd=fileparts(P{np});
        Cnames={SPM.xCon(:).name};
        if np==1,
            str='Select contrast(s)';
            disp(str);
            Ic=listdlg('promptstring',str,'selectionmode','multiple','liststring',Cnames);
            Icnames={Cnames{Ic}};
        end
        ic=[];ok=1;for n1=1:length(Icnames),temp=strmatch(Icnames{n1},Cnames,'exact');if numel(temp)~=1,ok=0;break;else ic(n1)=temp;end;end
        if ~ok, error(['the target contrasts are not found inside ',P{np}]); end
        spm_ss_threshold('subject',SPM,ic(:));
    end
    thr_p=spm_ss_threshold('end');
    txt=['Optimal FDR-corrected threshold: FDR-p < ',num2str(thr_p)];
    disp(txt);
    msgbox(txt);
    return;
end

switch(lower(option))
    case 'begin',
        SPMfiles={};
        ContrastIndexes={};
        ExpIndexes={};
        MaskFilenames={};
        OptimalThreshold={};
        if nargin>1,
            Threshold_type=spm_data;
            Threshold_p=Ic;
            if nargin>3, Conjunction_type=Ec;
            else Conjunction_type='and';
            end
        else
            Threshold_type={'automatic'};
            Threshold_p=nan;
            Conjunction_type='and';
        end
        AnyAutomatic=~isempty(strmatch('automatic',Threshold_type,'exact'));
        
    case 'subject',
        if size(Ic,2)>numel(Threshold_type),Threshold_type={Threshold_type{min(numel(Threshold_type),1:size(Ic,2))}}; end
        if size(Ic,2)>numel(Threshold_p),Threshold_p=Threshold_p(min(numel(Threshold_type),1:size(Ic,2))); end
        if AnyAutomatic
            % estimate one threshold value separately for each column of Ic
            for nic2=1:size(Ic,2),
                SPM=spm_data.SPM{Ec(nic2)};
                if strcmpi(Threshold_type{nic2},'automatic')
                    for nic1=1:size(Ic,1),
                        a=spm_vol(fullfile(SPM.swd,SPM.xCon(abs(Ic(nic1,nic2))).Vspm.fname));
                        b=spm_read_vols(a);
                        idx=find(~isnan(b)&b~=0);
                        dof=[SPM.xCon(abs(Ic(nic1,nic2))).eidf,SPM.xX.erdf];
                        STAT=SPM.xCon(abs(Ic(nic1,nic2))).STAT;
                        Y=nan+zeros(size(b));
                        switch(STAT),
                            case 'Z',Y(idx)=1-spm_Ncdf(b(idx));
                            case 'T',Y(idx)=1-spm_Tcdf(b(idx),dof(2));
                            case 'X',Y(idx)=1-spm_Xcdf(b(idx),dof(2));
                            case 'F',Y(idx)=1-spm_Fcdf(b(idx),dof);
                            otherwise, error('null');
                        end
                        Y(:)=spm_ss_fdr(Y(:));
                        N=numel(Y);
                        [nill,idx]=sort(Y(:));
                        Q=zeros(N,1);Q(idx)=((1:N)'/N).*((1-Y(idx)).^1);
                        [nill,idxt]=max(Q);
                        if numel(OptimalThreshold)<nic2,OptimalThreshold{nic2}=[];end
                        OptimalThreshold{nic2}(end+1)=Y(idxt);
                    end
                end
            end
        else
            MaskFilenames{end+1}=spm_ss_createlocalizermask(spm_data.SPM,Ic,Ec,overwrite,Threshold_type,Threshold_p,Conjunction_type,options,optionROIs);
        end
        for nexp=1:numel(spm_data.SPM),
            idx=find(Ec==nexp);
            spm_data.SPM{nexp}.xCon=spm_data.SPM{nexp}.xCon(abs(Ic(:,idx)));
            Ic(:,idx)=sign(Ic(:,idx)).*reshape(1:numel(spm_data.SPM{nexp}.xCon),[size(Ic,1),numel(idx)]);
        end
        SPMfiles{end+1}=spm_data.SPM;%fullfile(SPM.swd,'SPM.mat');
        ContrastIndexes{end+1}=Ic;
        ExpIndexes{end+1}=Ec;
    case 'end',
        for nic2=1:numel(Threshold_type),
            if strcmpi(Threshold_type{nic2},'automatic')
                Threshold_p(nic2)=mean(OptimalThreshold{nic2}(~isnan(OptimalThreshold{nic2})));
                Threshold_type{nic2}='FDR';
                disp(['Optimal FDR-corrected threshold: FDR-p < ',num2str(Threshold_p(nic2))]);
            end
        end
        if nargout>1&&AnyAutomatic
            for np=1:numel(SPMfiles),
                %load(SPMfiles{np},'SPM');
                MaskFilenames{np}=spm_ss_createlocalizermask(SPMfiles{np},ContrastIndexes{np},ExpIndexes{np},overwrite,Threshold_type,Threshold_p,Conjunction_type);
            end
        end
        thr_p=Threshold_p;
        thr_type=Threshold_type;
        filenames=MaskFilenames;        
end

