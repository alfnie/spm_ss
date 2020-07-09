function ss=spm_ss_overlap(ss)
% SPM_SS_overlap computes overlap measures pair-wise among multiple contrasts 
% ss=spm_ss_overlap;     Defines a new analysis entirely from the gui
% ss=spm_ss_overlap(ss); Defines a new analysis with the gui asking only those aspects of the design that have not been specified
%
% Fields of ss structure:
%
% ss.swd                Directory for output analysis files
% ss.n                  Number of subjects
% (when selecting by contrast names; default)
%    ss.files_spm                cell array of SPM.mat files (one per subject; ss.files_spm{nsubject} is a string pointing to a 1st-level SPM.mat file); In addition, in the case of multiple effect of interest contrasts resulting from different experiments ss.files_spm{neffect,nsubject} can be used to point to each contrast-specific subject-specific SPM.mat file.
%    ss.Localizer_contrasts      cell array of contrasts-of-interest names (ss.Localizer_contrasts{neffect} is a char array containing one or several 1st-level contrast name(s); for manual cross-validation use ss.Localizer_contrast{ncrossvalid,neffect} and the toobox will perform cross-validation across these partitions);
%    ss.Localizer_thr_p          false positive threshold
%    ss.Localizer_thr_type       multiple comparisons correction type ('FDR','FWE','none','automatic') (default 'FDR')
% (when selecting by contrast files; i.e. ss.files_selectmanually=1)
%    ss.Localizer                cell array of contrasts-of-interest contrast filenames (ss.Localizer{nsubject}{neffect} is a char array containing a con*.img contrast file)  
% ss.ExplicitMasking    explicit mask file name (only voxels where the mask takes values above 0 will be considered in any analysis; default [])
% ss.ManualROIs         manually-defined ROI file name (ROI image should contain integer numbers, from 1 to m, where m is the number of ROIs; default [])
%                       when using subject-specific ROI files, ss.ManualROIs is a cell array with ss.ManualROIs{i} defining the ROI file name for subject i 
% ss.ask                gui interaction level: 'none' (any missing information is assumed to take default values), 'missing' (any missing information will be asked to the user), 'all' (it will ask for confirmation on each parameter)
% 

cwd=pwd;
initestimation=0;
posstr=1;
if nargin<1, 
    ss=[]; 
    if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
    str='Specify a new model?';
    disp(str);
    newmodel=spm_input(str,posstr,'m','Specify a new model|Modify an existing model',[1,2], 1);posstr='+1';
    initestimation=1;
    if newmodel==2,
        str={'Warning!: Modifying a model removes any estimated model parameters on the new model',...
            'Select a different output folder if you do not want to loose any currently estimated model parameters',...
            'Are you sure you want to continue?'};
        if spm_input(str,posstr,'bd','stop|continue',[1,0],1),return;end;posstr='+1';
        str='Select spm_ss_overlap.mat analysis file';
        disp(str);
        Pdefault='';objname=findobj('tag','spm_ss');if numel(objname)==1,objdata=get(objname,'userdata');if isfield(objdata,'files_spm_ss'),Pdefault=objdata.files_spm_ss;end;end;
        P=spm_select(1,'^SPM_ss_overlap.*\.mat$',str,{Pdefault});
        if numel(objname)==1&&~isempty(P),objdata.files_spm_ss=P;set(objname,'userdata',objdata);end;
        load(P);
        ss.swd=fileparts(P);
        ss.ask='all';
    end
end
if ~isfield(ss,'files_selectmanually')||isempty(ss.files_selectmanually), ss.files_selectmanually=0;end
if ~isfield(ss,'ask')||isempty(ss.ask), 
    ss.ask='missing'; ss.askn=1;
else
    types={'none','missing','all'};typesn=[0,1,2];
    if isnumeric(ss.ask), sstype=ss.ask;
    else sstype=strmatch(lower(ss.ask),lower(types),'exact'); end
    ss.ask=types{sstype};
    ss.askn=typesn(sstype);
end

ss.type='overlap';
if ss.askn>1||~isfield(ss,'swd')||isempty(ss.swd), 
    if ~isfield(ss,'swd')||isempty(ss.swd), ss.swd=''; end
    if ss.askn,ss.swd=spm_select(1,'Dir','Select directory for output analysis files',{ss.swd}); end
else
    if ~isdir(ss.swd),
        [swdp,swdn,swde]=fileparts(ss.swd);
        [ok,nill]=mkdir(swdp,[swdn,swde]);
    end
end


if ~ss.files_selectmanually,
    if isfield(ss,'files_spm')&&~isempty(ss.files_spm)&&(~isfield(ss,'Localizer_spm')||isempty(ss.Localizer_spm)),% batch files_spm use
        ss.Localizer_spm=ss.files_spm;
    end
    if ss.askn>1||~isfield(ss,'Localizer_spm')||isempty(ss.Localizer_spm),
        if ~isfield(ss,'Localizer_spm')||isempty(ss.Localizer_spm), ss.Localizer_spm={}; end
        str='Select first-level SPM.mat(s) (one per subject, containing LOCALIZER contrasts) - or hit cancel if you prefer to select mask files manually';
        disp(str);
        Pdefault={''};objname=findobj('tag','spm_ss');if numel(objname)==1,objdata=get(objname,'userdata');if isfield(objdata,'files_spm'),Pdefault=objdata.files_spm;end;end;
        if ~isempty(ss.Localizer_spm),Pdefault=ss.Localizer_spm;end
        if ss.askn,P=cellstr(spm_select(inf,'^SPM\.mat$',str,Pdefault));else P=Pdefault; end
        if numel(objname)==1&&~isempty(P)&&~isempty(P{1}),objdata.files_spm=P;set(objname,'userdata',objdata);end;
        if isempty(P)||isempty(P{1}), ss.files_selectmanually=1; 
        else
            if numel(ss.Localizer_spm)~=numel(P), ss.Localizer={}; ss.n=numel(P); else for n1=1:numel(ss.Localizer_spm), if ~strcmp(ss.Localizer_spm{n1},P{n1}), ss.Localizer={}; ss.n=numel(P); end; end; end
            ss.Localizer_spm=P; 
            %ss.n=numel(P); 
            ss.files_selectmanually=0; 
        end
    end
end
if ~ss.files_selectmanually&&isfield(ss,'n')&&~isempty(ss.n),
    if ~rem(numel(ss.Localizer_spm),ss.n),ss.Localizer_spm=reshape(ss.Localizer_spm,numel(ss.Localizer_spm)/ss.n,ss.n);end
end
if ~ss.files_selectmanually,%&&(~isfield(ss,'n')||isempty(ss.n)),
    ss.n=size(ss.Localizer_spm,2); 
end

if ss.files_selectmanually,
    if ss.askn>1||~isfield(ss,'n')||isempty(ss.n), 
        if ~isfield(ss,'n')||isempty(ss.n), ss.n=[]; end
        if ss.askn, 
            if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
            str='Number of subjects?';
            disp(str);
            ss.n=spm_input(str,posstr,'n',ss.n,1);posstr='+1'; 
        end
    end
end

if ss.askn>1||~isfield(ss,'Localizer')||isempty(ss.Localizer),
    if ~isfield(ss,'Localizer')||isempty(ss.Localizer),ss.Localizer={};end
    if ss.files_selectmanually,
        for np=1:ss.n,
            if numel(ss.Localizer)<np,ss.Localizer{np}={''};end
            str=['Select LOCALIZER mask volumes for subject #',num2str(np),'(one mask file, or multiple mask files if using crossvalidation)'];
            disp(str);
            ss.Localizer{np}=cellstr(spm_select(inf,'image',str,ss.Localizer{np}));
        end
    else
        if ~isfield(ss,'Localizer_contrasts')||isempty(ss.Localizer_contrasts),ss.Localizer_contrasts={};end
        nnc=size(ss.Localizer_spm,1);
        if nnc>1,
            %ss.Localizer_contrasts={};
            spm_data=[];Ec=[];
            for nc=1:nnc,
                current_spm=ss.Localizer_spm{min(size(ss.Localizer_spm,1),nc),1};
                [spm_data,SPM,Ec(nc)]=spm_ss_importspm(spm_data,current_spm);
%                 load(ss.Localizer_spm{nc},'SPM');
%                 SPM.swd=fileparts(ss.Localizer_spm{nc});
                Cnames={SPM.xCon(:).name};
                for ncr=1:max(1,size(ss.Localizer_contrasts,1)),
                    ic=[];for n1=nc:min(nc,size(ss.Localizer_contrasts,2)),temp=strmatch(ss.Localizer_contrasts{ncr,n1},Cnames,'exact');if numel(temp)~=1&&~isempty(strmatch('-not',ss.Localizer_contrasts{ncr,n1})),temp=-strmatch(ss.Localizer_contrasts{ncr,n1}(5:end),Cnames,'exact');end;if numel(temp)~=1,break;else ic=temp;end;end
                    if ss.askn>1||isempty(ic),
                        str=['Select LOCALIZER contrast #',num2str(nc)];
                        disp(str);
                        Ic=listdlg('promptstring',str,'selectionmode','single','liststring',Cnames,'initialvalue',ic); %Ic=spm_conman(SPM,'T|F',inf,str,'',0);
                        str='Conjunction: Inclusive (1) or exclusive (0) for this contrast?';
                        disp(str);
                        disp(char({Cnames{Ic}}));
                        signIc=spm_input(str,posstr,'r',num2str(sign(ic(:)')>0),[1,numel(Ic)],[0,1]);posstr='+1';
                        Ic(signIc<.5)=-Ic(signIc<.5);
                    else
                        Ic=ic;
                    end
                    if numel(ic)~=numel(Ic) || any(ic(:)~=Ic(:)), ss.Localizer={}; end
                    ss.Localizer_contrasts{ncr,nc}=Cnames{abs(Ic)}; % note: the gui does not allow to manually cross-validate when using this option (i.e. multiple contrasts are interpreted as multiple effects of interest contrasts, not as multiple 'sessions'), use batch scripts instead
                    if sign(Ic)<0, ss.Localizer_contrasts{ncr,nc}=['-not',ss.Localizer_contrasts{ncr,nc}];end %adds -not prefix to exclusive contrasts
                end
            end
        else
            load(ss.Localizer_spm{1},'SPM');
            SPM.swd=fileparts(ss.Localizer_spm{1});
            Cnames={SPM.xCon(:).name};
            for ncr=1:max(1,size(ss.Localizer_contrasts,1)),
                ic=[];for n1=1:size(ss.Localizer_contrasts,2),temp=strmatch(ss.Localizer_contrasts{ncr,n1},Cnames,'exact');if numel(temp)~=1&&~isempty(strmatch('-not',ss.Localizer_contrasts{ncr,n1})),temp=-strmatch(ss.Localizer_contrasts{ncr,n1}(5:end),Cnames,'exact');end;if numel(temp)~=1,ic=[];break;else ic(n1)=temp;end;end
                if ss.askn>1||isempty(ic),
                    str='Select LOCALIZER contrast(s)';
                    disp(str);
                    Ic=listdlg('promptstring',str,'selectionmode','multiple','liststring',Cnames,'initialvalue',abs(ic)); %Ic=spm_conman(SPM,'T|F',inf,str,'',0);
                    if numel(Ic)>1,
                        str='Conjunction: Inclusive (1) or exclusive (0) mask for each contrast?';
                        disp(str);
                        disp(char({Cnames{Ic}}));
                        signIc=spm_input(str,posstr,'r',num2str(sign(ic(:)')>0),[1,numel(Ic)],[0,1]);posstr='+1';
                        Ic(signIc<.5)=-Ic(signIc<.5);
                    end
                else Ic=ic; end
                for n1=1:numel(Ic), ss.Localizer_contrasts{ncr,n1}=Cnames{abs(Ic(n1))}; end % note: the gui does not allow to manually cross-validate when using this option (i.e. multiple contrasts are interpreted as a conjunction of the selected localizer contrasts, not as multiple 'sessions'), use batch scripts instead
                for n1=1:numel(Ic),if sign(Ic(n1))<0, ss.Localizer_contrasts{ncr,n1}=['-not',ss.Localizer_contrasts{ncr,n1}];end;end %adds -not prefix to exclusive contrasts
                if numel(ic)~=numel(Ic) || any(ic(:)~=Ic(:)), ss.Localizer={}; end
            end
        end
    end
end

if ss.files_selectmanually&&(ss.askn>1||~isfield(ss,'Localizer_contrasts')||isempty(ss.Localizer_contrasts)),
    for nc=1:numel(ss.Localizer{1}),
        if numel(ss.Localizer_contrasts)<nc,ss.Localizer_contrasts{nc}=['Effect #',num2str(nc)];end
        if ss.askn,
            if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
            str=['Effect #',num2str(nc),'name?'];
            disp(str);
            ss.Localizer_contrasts{nc}=spm_input(str,posstr,'s',ss.Localizer_contrasts{nc});posstr='+1';
        end
    end
end

if isfield(ss,'Localizer_thr_type')&&ischar(ss.Localizer_thr_type),ss.Localizer_thr_type=cellstr(ss.Localizer_thr_type);end
if ~ss.files_selectmanually&&(ss.askn>1||~isfield(ss,'Localizer_thr_type')||numel(ss.Localizer_thr_type)<size(ss.Localizer_contrasts,2)||~isfield(ss,'Localizer_thr_p')||numel(ss.Localizer_thr_p)<size(ss.Localizer_contrasts,2)),
    if ~isfield(ss,'Localizer_thr_type')||isempty(ss.Localizer_thr_type), ss.Localizer_thr_type=repmat({'FDR'},[1,size(ss.Localizer_contrasts,2)]); end
    if numel(ss.Localizer_thr_type)~=size(ss.Localizer_contrasts,2), ss.Localizer_thr_type={ss.Localizer_thr_type{min(numel(ss.Localizer_thr_type),1:size(ss.Localizer_contrasts,2))}}; end
    if ~isfield(ss,'Localizer_thr_p')||isempty(ss.Localizer_thr_p),
        ss.Localizer_thr_p=[];for n1=1:numel(ss.Localizer_thr_type), if strcmpi(ss.Localizer_thr_type{n1},'none'),ss.Localizer_thr_p(n1)=.001;elseif strcmpi(ss.Localizer_thr_type{n1},'automatic'), ss.Localizer_thr_p(n1)=nan; else ss.Localizer_thr_p(n1)=.05;end; end
    end
    if numel(ss.Localizer_thr_p)~=size(ss.Localizer_contrasts,2), ss.Localizer_thr_p=ss.Localizer_thr_p(min(numel(ss.Localizer_thr_p),1:size(ss.Localizer_contrasts,2))); end
    if ss.askn,
        if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
        sstype1={};sstype2=[];
        for n1=1:size(ss.Localizer_contrasts,2)
            types={'FDR','FWE','none','automatic'};
            sstype=strmatch(lower(ss.Localizer_thr_type{n1}),lower(types),'exact');
            str=['Contrast #',num2str(n1),' localizer p value adjustment to control? (',ss.Localizer_contrasts{1,n1},')'];
            disp(str);
            sstype1{n1}=types{spm_input(str,posstr,'m','FDR|FWE|none|automatic',[],sstype)};posstr='+1';
            if strcmpi(sstype1{n1},'automatic'), sstype2(n1)=nan;
            else
                str=['Contrast #',num2str(n1),' localizer p value threshold? (',ss.Localizer_contrasts{1,n1},')'];
                disp(str);
                sstype2(n1)=spm_input(str,posstr,'r',ss.Localizer_thr_p(n1)); posstr='+1';
            end
        end
        if numel(ss.Localizer_thr_type)~=numel(sstype1), ss.Localizer={}; else for n1=1:numel(ss.Localizer_thr_type),if ~strcmpi(ss.Localizer_thr_type{n1},sstype1{n1}), ss.Localizer={}; end; end; end
        if any(ss.Localizer_thr_p~=sstype2), ss.Localizer={}; end
        ss.Localizer_thr_type=sstype1;
        ss.Localizer_thr_p=sstype2;
    end
end

if ss.askn>1||~isfield(ss,'ExplicitMasking'),
    if ~isfield(ss,'ExplicitMasking')||isempty(ss.ExplicitMasking), ss.ExplicitMasking=''; end
    str='Explicit masking? (this is an optional aditional mask used across all subjects)';
    disp(str);
    sstype=spm_input(str,posstr,'b','none|select',[1,2],1+~isempty(ss.ExplicitMasking));posstr='+1';
    if sstype==1,ss.ExplicitMasking='';
    else
        str='Select mask file';
        disp(str);
        ss.ExplicitMasking=spm_select([0,1],'image',str,{ss.ExplicitMasking});
    end
end

if (ss.askn>1||~isfield(ss,'ManualROIs')),
    if ~isfield(ss,'ManualROIs')||isempty(ss.ManualROIs), ss.ManualROIs=''; end
    str='Type of ROI file?';
    disp(str);
    sstype=spm_input(str,posstr,'m','Single ROI file|Subject-specific ROI file',[1,2], sstype); posstr='+1';
    if sstype==1
        str='Select ROI file (ROIs are labeled as integer numbers within this volume)';
        disp(str);
        ss.ManualROIs=spm_select(1,'image',str,cellstr(ss.ManualROIs));
    else
        str='Select ROI files (one file per subject; ROIs are labeled as integer numbers within each volume)';
        disp(str);
        ss.ManualROIs=cellstr(spm_select(ss.n,'image',str,cellstr(ss.ManualROIs)));
    end
end
if ~isfield(ss,'ManualROIs'), ss.ManualROIs=''; end
if ~iscell(ss.ManualROIs), ss.ManualROIs=cellstr(ss.ManualROIs); end

% imports contrast files
if ~ss.files_selectmanually&&(~isfield(ss,'Localizer')||isempty(ss.Localizer)),
    initestimation=1;
    ss.Localizer={};
    disp('Importing contrast information from SPM.mat files. Please wait...');
    spm_ss_threshold('begin',ss.Localizer_thr_type,ss.Localizer_thr_p);
    for np=1:size(ss.Localizer_spm,2),%subjects
        spm_data=[];
        Ic2=[];Ec2=[];ok=1;
        for ncontrast=1:size(ss.Localizer_contrasts,2),%contrasts
            current_spm=ss.Localizer_spm{min(size(ss.Localizer_spm,1),ncontrast),np};
            [spm_data,SPM,Ec2(ncontrast)]=spm_ss_importspm(spm_data,current_spm);
            Cnames={SPM.xCon(:).name};
            for ncrossvalid=1:size(ss.Localizer_contrasts,1),%crossvalidation partitions
                temp=strmatch(ss.Localizer_contrasts{ncrossvalid,ncontrast},Cnames,'exact');if numel(temp)~=1&&~isempty(strmatch('-not',ss.Localizer_contrasts{ncrossvalid,ncontrast})),temp=-strmatch(ss.Localizer_contrasts{ncrossvalid,ncontrast}(5:end),Cnames,'exact');end;if numel(temp)~=1,ok=0;else Ic2(ncrossvalid,ncontrast)=temp;end
                if ~ok, error(['the target contrasts are not found inside ',current_spm]); end
            end
        end
        ncross=size(Ic2,1);
        nconjunction=size(Ic2,2);
        spm_ss_threshold('subject',spm_data,Ic2,Ec2,'noconjunction',ss.ManualROIs{min(numel(ss.ManualROIs),np)});
        %spm_ss_threshold('subject',spm_data,Ic2,Ec2,'noconjunction');
        fprintf('.');
    end
    fprintf('\n');
    [ss.Localizer_thr_p,ss.Localizer_thr_type,ss.Localizer]=spm_ss_threshold('end');
end

if initestimation, ss.overlap={}; end
ss.ask=[];ss.askn=[];

objname=findobj('tag','spm_ss');if numel(objname)==1,objdata.files_spm_ss=fullfile(ss.swd,['SPM_ss_',ss.type,'.mat']);set(objname,'userdata',objdata);end;
save(fullfile(ss.swd,['SPM_ss_',ss.type,'.mat']),'ss');
disp(['Analysis file saved: ',fullfile(ss.swd,['SPM_ss_',ss.type,'.mat'])]);


%% overlap computation

% % explicit mask
if ~isempty(ss.ExplicitMasking),XM=spm_vol(ss.ExplicitMasking);else XM=[];end

% creates transformed measures
k=0;
ss.PV=[];
neffects=[];
ncross=[];
idxPV={};
for n=1:ss.n,
    [ncross(n),temp]=size(ss.Localizer{n});
    if ~isempty(neffects)&&neffects~=temp, error(['mismatched number of contrasts in subject',num2str(n)]); return;
    else neffects=temp; end
    for ne=1:neffects,
        for nc=1:ncross(n),
            [pth1,nm1,ext1,num1]=spm_fileparts(ss.Localizer{n}{nc,ne});
            Nvolume=[nm1,ext1,num1];
            k=k+1;
            ss.PN{k}=fullfile(pth1,Nvolume);
            ss.PV(n,k)=ne;
        end
        temp=find(ss.PV(n,:)>0);
        idxPV{n,ne}=find(ss.PV(n,temp)==ne);
    end
end

cd(ss.swd);
ss.VN=spm_vol(char(ss.PN));

ss.PM=ss.ManualROIs;
ss.VM=spm_vol(char(ss.PM));

% analysis
extname=['_',ss.type];

fprintf('Performing overlap estimation...');
for n=1:ss.n,
    fprintf(1,'.');
    
    idxk=ss.PV(n,:)>0;
    idxk1=find(idxk,1);
    clear XYZ; [XYZ{1},XYZ{2},XYZ{3}]=ndgrid(1:ss.VN(idxk1).dim(1),1:ss.VN(idxk1).dim(2),1:ss.VN(idxk1).dim(3));
    XYZ=reshape(cat(4,XYZ{:}),[],3)';XYZ=ss.VN(idxk1).mat*cat(1,XYZ,ones(1,size(XYZ,2))); % voxel-coordinates for first functional volume for subject n
    vm=ss.VM(min(numel(ss.VM),n));
    frois=reshape(round(spm_get_data(vm,pinv(vm.mat)*XYZ)),ss.VN(idxk1).dim(1:3));
    if n==1, 
        nrois=max(frois(:));
        Nplane=nan+zeros([ss.n,nrois]);
        Qplane=nan+zeros([ss.n,neffects,neffects,nrois]); 
    else
        if max(frois(:))~=nrois, error('All subject-specific ROI files should contain the same number of ROIs'); end
    end
    if ~isempty(XM),
        tXM=reshape(double(spm_get_data(XM,pinv(XM.mat)*XYZ)>0),ss.VN(idxk1).dim(1:3));
        frois=tXM.*frois;
    end
    
    for nroi=1:nrois,
        idx=find(frois==nroi);
        [idx1,idx2,idx3]=ind2sub(size(frois),idx);
        xyz=[idx1,idx2,idx3,ones(numel(idx1),1)]';
        Nplane(n,nroi)=numel(idx);
    
        N=spm_get_data(ss.VN(idxk),xyz); % note: all volumes for each subject are expected to be in the same space
        N=N>0;
        for ne1=1:neffects,
            idx1=idxPV{n,ne1};
            for ne2=1:neffects,
                idx2=idxPV{n,ne2};
                Qplane(n,ne1,ne2,nroi)=mean(sum(N(idx1,:)&N(idx2,:),2),1);
            end
        end
    end
end
fprintf(1,'\n');
cd(cwd);

% save files
ss.estimate=struct('voxels',Nplane,'overlap',Qplane); 
save(fullfile(ss.swd,['SPM_ss',extname,'.mat']),'ss');
disp(['Analysis file saved: ',fullfile(ss.swd,['SPM_ss',extname,'.mat'])]);

fname=['spm_ss',extname,'_data.csv'];
fh=fopen(fullfile(ss.swd,fname),'wt');
fprintf(fh,'Overlap analyses (number-of-voxels intersection for each pair of localizer contrasts)\n');
fprintf(fh,'\n');
for ne=1:size(ss.Localizer_contrasts,2), fprintf(fh,'Contrast#%d',ne); if ne<size(ss.Localizer_contrasts,2), fprintf(fh,','); else fprintf(fh,'\n'); end; end
for nc=1:size(ss.Localizer_contrasts,1), 
    for ne=1:size(ss.Localizer_contrasts,2), fprintf(fh,'%s',ss.Localizer_contrasts{nc,ne}); if ne<size(ss.Localizer_contrasts,2), fprintf(fh,','); else fprintf(fh,'\n'); end; end
end
fprintf(fh,'\n');
fprintf(fh,'ROI#,');
for ns=1:ss.n,fprintf(fh,'Subject#%d size,',ns); for ne1=1:neffects,for ne2=ne1:neffects,fprintf(fh,'Subject#%d[%d;%d]',ns,ne1,ne2); if ne1<neffects||ne2<neffects||ns<ss.n, fprintf(fh,','); else fprintf(fh,'\n'); end; end; end; end
for nroi=1:nrois,
    fprintf(fh,'%d,',nroi);
    for ns=1:ss.n, fprintf(fh,'%d,',round(ss.estimate.voxels(ns,nroi))); for ne1=1:neffects,for ne2=ne1:neffects,fprintf(fh,'%f',ss.estimate.overlap(ns,ne1,ne2,nroi)); if ne1<neffects||ne2<neffects||ns<ss.n, fprintf(fh,','); else fprintf(fh,'\n'); end; end; end; end
end
fclose(fh);

end


