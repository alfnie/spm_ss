function bcc=spm_bcc_design(bcc)
% SPM_BCC_DESIGN gui defining a multi-subject between-condition-correlation analysis 
% bcc=spm_bcc_design;     Defines a new analysis entirely from the gui
% bcc=spm_bcc_design(bcc); Defines a new analysis with the gui asking only those aspects of the design that have not been specified
%
% Fields of bcc structure:
%
% bcc.swd                Directory for output analysis files
% bcc.n                  Number of subjects
% (when selecting by contrast names; default)
%    bcc.files_spm                    cell array of SPM.mat files (one per subject; bcc.files_spm{nsubject} is a string pointing to a 1st-level SPM.mat file); alternatively you can use the fields: bcc.EffectOfInterest1_spm and bcc.EffectOfInterest2_spm if the effect of interest contrasts have been defined in different SPM.mat files (for the same subjects). 
%    bcc.EffectOfInterest1_contrasts  cell array of Effects-of-interest contrast names (bcc.EffectOfInterest_contrasts{neffect} is a char array containing one or several 1st-level contrast name(s); for manual cross-validation use bcc.EffectOfInterest_contrasts{ncrossvalid,neffect} and the toobox will perform cross-validation across these partitions) 
%    bcc.EffectOfInterest2_contrasts  cell array of Effects-of-interest contrast names (bcc.EffectOfInterest_contrasts{neffect} is a char array containing one or several 1st-level contrast name(s); for manual cross-validation use bcc.EffectOfInterest_contrasts{ncrossvalid,neffect} and the toobox will perform cross-validation across these partitions) 
%      note: (automatic cross-validation) the toolbox will automatically check the orthogonality between the EOI1 and EOI2. If these are found not to be orthogonal the toolbox will create and estimate new (now orthogonal) contrasts by partitioning the selected contrasts across sessions
% (when selecting by contrast files; e.g. bcc.files_selectmanually(1)==1 (for EffectOfInterest1) or bcc.files_selectmanually(2)= 1 (for EffectOfInterest2) 
%    bcc.EffectOfInterest1   cell array of Effects-of-interest contrast files (bcc.EffectOfInterest1{nsubject}{neffect,ncrossvalid} is a char array containing a con*.img contrast file)  
%    bcc.EffectOfInterest2   cell array of Effects-of-interest contrast files (bcc.EffectOfInterest2{nsubject}{neffect,ncrossvalid} is a char array containing a con*.img contrast file)  
% bcc.ExplicitMasking    explicit mask file name (only voxels where the mask takes values above 0 will be considered in any analysis; default [])
% bcc.ManualROIs         manually-defined ROI file name; if left empty correlations will be computed across the entire brain (ROI image should contain integer numbers, from 1 to m, where m is the number of ROIs);
%                        when using subject-specific ROI files, bcc.ManualROIs is a cell array with bcc.ManualROIs{i} defining the ROI file name for subject i 
% bcc.model              Between-subjects model type (1: one-sample t-test; 2: two-sample t-test; 3: multiple regression) (note: this field is disregarded if the design matrix bcc.X below is directly defined) 
% bcc.X                  Between-subjects design matrix [nsubjects,nregressors]
% bcc.Xname              Cell array of Between-subjects regressor names [nregressors]
% bcc.C                  Array of structures containing contrast information
% bcc.C(k).between       Between-subjects contrast(s) vector/matrix [m,nregressors]
% bcc.C(k).within        Within-subjects contrast(s) vector/matrix [n,neffects1*neffects2]
% bcc.C(k).name          Contrast name 
% bcc.ask                gui interaction level: 'none' (any missing information is assumed to take default values), 'missing' (any missing information will be asked to the user), 'all' (it will ask for confirmation on each parameter)
% 

% spm_ss initialization
spm_ss init silent;

initestimation=0;
posstr=1;
if nargin<1, 
    bcc=[]; 
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
        str='Select spm_bcc*.mat analysis file';
        disp(str);
        Pdefault='';objname=findobj('tag','spm_bcc');if numel(objname)==1,objdata=get(objname,'userdata');if isfield(objdata,'files_spm_bcc'),Pdefault=objdata.files_spm_bcc;end;end;
        P=spm_select(1,'^SPM_bcc.*\.mat$',str,{Pdefault});
        if numel(objname)==1&&~isempty(P),objdata.files_spm_bcc=P;set(objname,'userdata',objdata);end;
        load(P);
        bcc.swd=fileparts(P);
        bcc.ask='all';
    end
end
if ~isfield(bcc,'files_selectmanually')||isempty(bcc.files_selectmanually), bcc.files_selectmanually=[0 0];end
if numel(bcc.files_selectmanually)==1, bcc.files_selectmanually=bcc.files_selectmanually([1 1]); end

if ~isfield(bcc,'ask')||isempty(bcc.ask), 
    bcc.ask='missing'; bcc.askn=1;
else
    types={'none','missing','all'};typesn=[0,1,2];
    if isnumeric(bcc.ask), sstype=bcc.ask;
    else sstype=strmatch(lower(bcc.ask),lower(types),'exact'); end
    bcc.ask=types{sstype};
    bcc.askn=typesn(sstype);
end

bcc.type='';

if bcc.askn>1||~isfield(bcc,'swd')||isempty(bcc.swd), 
    if ~isfield(bcc,'swd')||isempty(bcc.swd), bcc.swd=''; end
    if bcc.askn,bcc.swd=spm_select(1,'Dir','Select directory for output analysis files',{bcc.swd}); end
else
    if ~isdir(bcc.swd),
        [swdp,swdn,swde]=fileparts(bcc.swd);
        [ok,nill]=mkdir(swdp,[swdn,swde]);
    end
end

if all(~bcc.files_selectmanually),
    if isfield(bcc,'files_spm')&&~isempty(bcc.files_spm)&&(~isfield(bcc,'EffectOfInterest1_spm')||isempty(bcc.EffectOfInterest1_spm)||~isfield(bcc,'EffectOfInterest2_spm')||isempty(bcc.EffectOfInterest2_spm)),% batch files_spm use
        bcc.EffectOfInterest1_spm=bcc.files_spm;
        bcc.EffectOfInterest2_spm=bcc.files_spm;
    end
end

if ~bcc.files_selectmanually(1),
    if bcc.askn>1||~isfield(bcc,'EffectOfInterest1_spm')||isempty(bcc.EffectOfInterest1_spm),
        if ~isfield(bcc,'EffectOfInterest1_spm')||isempty(bcc.EffectOfInterest1_spm), bcc.EffectOfInterest1_spm={}; end
        str='Select first-level SPM.mat files (one per subject, containing EFFECTS OF INTEREST1 contrasts) - or hit cancel if you prefer to select contrast files manually';
        disp(str);
        Pdefault={''};objname=findobj('tag','spm_bcc');if numel(objname)==1,objdata=get(objname,'userdata');if isfield(objdata,'files_spm'),Pdefault=objdata.files_spm;end;end;
        if ~isempty(bcc.EffectOfInterest1_spm),Pdefault=bcc.EffectOfInterest1_spm;end
        if bcc.askn,P=cellstr(spm_select(inf,'^SPM\.mat$',str,Pdefault))';else P=Pdefault; end
        if numel(objname)==1&&~isempty(P)&&~isempty(P{1}),objdata.files_spm=P;set(objname,'userdata',objdata);end;
        if isempty(P)||isempty(P{1}), bcc.files_selectmanually(1)=1; 
        else
            if numel(bcc.EffectOfInterest1_spm)~=numel(P), bcc.EffectOfInterest2={}; bcc.EffectOfInterest1={}; bcc.n=numel(P); else for n1=1:numel(bcc.EffectOfInterest1_spm), if ~strcmp(bcc.EffectOfInterest1_spm{n1},P{n1}), bcc.EffectOfInterest2={}; bcc.EffectOfInterest1={}; bcc.n=numel(P); end; end; end
            bcc.EffectOfInterest1_spm=P; 
            %bcc.n=numel(P); 
            bcc.files_selectmanually(1)=0; 
            if ~bcc.files_selectmanually(2)&&(~isfield(bcc,'EffectOfInterest2_spm')||isempty(bcc.EffectOfInterest2_spm)),
                bcc.EffectOfInterest2_spm=bcc.EffectOfInterest1_spm; % (exception; when using gui first time assume same EffectOfInterest1/2 SPM.mat files)
            end
        end
    end
end
if ~bcc.files_selectmanually(2),
    if bcc.askn>1||~isfield(bcc,'EffectOfInterest2_spm')||isempty(bcc.EffectOfInterest2_spm),
        if ~isfield(bcc,'EffectOfInterest2_spm')||isempty(bcc.EffectOfInterest2_spm), bcc.EffectOfInterest2_spm={}; end
        str='Select first-level SPM.mat(s) (one per subject, containing EFFECTS OF INTEREST2 contrasts) - or hit cancel if you prefer to select mask files manually';
        disp(str);
        Pdefault={''};objname=findobj('tag','spm_bcc');if numel(objname)==1,objdata=get(objname,'userdata');if isfield(objdata,'files_spm'),Pdefault=objdata.files_spm;end;end;
        if ~isempty(bcc.EffectOfInterest2_spm),Pdefault=bcc.EffectOfInterest2_spm;end
        if bcc.askn,P=cellstr(spm_select(inf,'^SPM\.mat$',str,Pdefault));else P=Pdefault; end
        if numel(objname)==1&&~isempty(P)&&~isempty(P{1}),objdata.files_spm=P;set(objname,'userdata',objdata);end;
        if isempty(P)||isempty(P{1}), bcc.files_selectmanually(2)=1; 
        else
            if numel(bcc.EffectOfInterest2_spm)~=numel(P), bcc.EffectOfInterest2={}; bcc.EffectOfInterest1={}; bcc.n=numel(P); else for n1=1:numel(bcc.EffectOfInterest2_spm), if ~strcmp(bcc.EffectOfInterest2_spm{n1},P{n1}), bcc.EffectOfInterest2={}; bcc.EffectOfInterest1={}; bcc.n=numel(P); end; end; end
            bcc.EffectOfInterest2_spm=P; 
            %bcc.n=numel(P); 
            bcc.files_selectmanually(2)=0; 
        end
    end
end
if ~any(bcc.files_selectmanually)&&isfield(bcc,'n')&&~isempty(bcc.n),
    if ~rem(numel(bcc.EffectOfInterest2_spm),bcc.n),bcc.EffectOfInterest2_spm=reshape(bcc.EffectOfInterest2_spm,numel(bcc.EffectOfInterest2_spm)/bcc.n,bcc.n);end
    if ~rem(numel(bcc.EffectOfInterest1_spm),bcc.n),bcc.EffectOfInterest1_spm=reshape(bcc.EffectOfInterest1_spm,numel(bcc.EffectOfInterest1_spm)/bcc.n,bcc.n);end
end
if ~any(bcc.files_selectmanually),%&&(~isfield(bcc,'n')||isempty(bcc.n)),
    if size(bcc.EffectOfInterest1_spm,2)~=size(bcc.EffectOfInterest2_spm,2), error('Number of SPM.mat files must match in EffectOfInterest1_spm and EffectOfInterest2_spm (one pair of SPM.mat files per subject)'); end
    bcc.n=size(bcc.EffectOfInterest2_spm,2); 
end

if any(bcc.files_selectmanually),
    if bcc.askn>1||~isfield(bcc,'n')||isempty(bcc.n), 
        if ~isfield(bcc,'n')||isempty(bcc.n), bcc.n=[]; end
        if bcc.askn, 
            if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
            str='Number of subjects?';
            disp(str);
            bcc.n=spm_input(str,posstr,'n',bcc.n,1);posstr='+1'; 
        end
    end
end

default_cRcontrast={};
default_cRname={};
if bcc.askn>1||~isfield(bcc,'EffectOfInterest1')||isempty(bcc.EffectOfInterest1),
    if ~isfield(bcc,'EffectOfInterest1')||isempty(bcc.EffectOfInterest1),bcc.EffectOfInterest1={};end
    if bcc.files_selectmanually(1),
        for np=1:bcc.n,
            if numel(bcc.EffectOfInterest1)<np,bcc.EffectOfInterest1{np}={''};end
            str=['Select EFFECT OF INTEREST1 contrast volumes for subject #',num2str(np),'(one or several contrast files)'];
            disp(str);
            bcc.EffectOfInterest1{np}=cellstr(spm_select(inf,'image',str,bcc.EffectOfInterest1{np}));
        end
    else
        if ~isfield(bcc,'EffectOfInterest1_contrasts')||isempty(bcc.EffectOfInterest1_contrasts),bcc.EffectOfInterest1_contrasts={};end
        nc=1;
        spm_data=[];Ec=[];
        EffectOfInterest1_contrasts=bcc.EffectOfInterest1_contrasts;
        bcc.EffectOfInterest1_contrasts={};
        while nc<=size(bcc.EffectOfInterest1_spm,1),
            current_spm=bcc.EffectOfInterest1_spm{min(size(bcc.EffectOfInterest1_spm,1),nc),1};
            [spm_data,SPM,Ec(nc)]=spm_ss_importspm(spm_data,current_spm);
            Cnames={SPM.xCon(:).name};
            if size(bcc.EffectOfInterest1_spm,1)==1
                ic=[];ok=1;for n1=1:length(EffectOfInterest1_contrasts),temp=strmatch(EffectOfInterest1_contrasts{n1},Cnames,'exact');if numel(temp)==1,ic(n1)=temp; elseif EffectOfInterest1_contrasts{n1}(1)=='*', temp=strmatch(EffectOfInterest1_contrasts{n1}(2:end),Cnames,'exact'); if numel(temp)==1, ic(n1)=-temp; else ok=0;break; end; else ok=0;break; end; end
            else
                ic=[];ok=1;for n1=nc:min(nc,length(EffectOfInterest1_contrasts)),temp=strmatch(EffectOfInterest1_contrasts{n1},Cnames,'exact');if numel(temp)==1,ic=temp; elseif EffectOfInterest1_contrasts{n1}(1)=='*', temp=strmatch(EffectOfInterest1_contrasts{n1}(2:end),Cnames,'exact'); if numel(temp)==1, ic=-temp; else ok=0;break; end; else ok=0;break; end; end
            end
            if bcc.askn>1||~ok||isempty(ic),
                if size(bcc.EffectOfInterest1_spm,1)==1
                    str='Select EFFECT OF INTEREST1 contrasts';
                    disp(str);
                    Ic=listdlg('promptstring',str,'selectionmode','multiple','liststring',Cnames,'initialvalue',abs(ic)); %Ic=spm_conman(SPM,'T|F',inf,str,'',0);
                else
                    str=['Select EFFECT OF INTEREST1 contrast #',num2str(nc)];
                    disp(str);
                    Ic=listdlg('promptstring',str,'selectionmode','single','liststring',Cnames,'initialvalue',abs(ic)); %Ic=spm_conman(SPM,'T|F',inf,str,'',0);
                end
            else
                Ic=ic;
            end
            if numel(ic)==numel(Ic) && all(abs(ic)==abs(Ic)), Ic=abs(Ic).*sign(ic); end
            if numel(ic)~=numel(Ic) || any(ic(:)~=Ic(:)), bcc.EffectOfInterest1={}; end
            if any(Ic<0)
                error('This option not implemented in bcc analyses');
            end
            bcc.EffectOfInterest1_contrasts(nc+(1:numel(Ic))-1)=Cnames(Ic); % note: the gui does not allow to manually cross-validate when using this option (i.e. multiple contrasts are interpreted as multiple effects of interest contrasts, not as multiple 'sessions'), use batch scripts instead
            nc=nc+numel(Ic);
        end
    end
end
if bcc.files_selectmanually(1)&&(bcc.askn>1||~isfield(bcc,'EffectOfInterest1_contrasts')||isempty(bcc.EffectOfInterest1_contrasts)),
    for nc=1:numel(bcc.EffectOfInterest1{1}),
        if numel(bcc.EffectOfInterest1_contrasts)<nc,bcc.EffectOfInterest1_contrasts{nc}=['Effect #',num2str(nc)];end
        if bcc.askn,
            if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
            str=['Effect #',num2str(nc),'name?'];
            disp(str);
            bcc.EffectOfInterest1_contrasts{nc}=spm_input(str,posstr,'s',bcc.EffectOfInterest1_contrasts{nc});posstr='+1';
        end
    end
end

if bcc.askn>1||~isfield(bcc,'EffectOfInterest2')||isempty(bcc.EffectOfInterest2),
    if ~isfield(bcc,'EffectOfInterest2')||isempty(bcc.EffectOfInterest2),bcc.EffectOfInterest2={};end
    if bcc.files_selectmanually(2),
        for np=1:bcc.n,
            if numel(bcc.EffectOfInterest2)<np,bcc.EffectOfInterest2{np}={''};end
            str=['Select EFFECT OF INTEREST2 contrast volumes for subject #',num2str(np),'(one or several contrast files)'];
            disp(str);
            bcc.EffectOfInterest2{np}=cellstr(spm_select(inf,'image',str,bcc.EffectOfInterest2{np}));
        end
    else
        if ~isfield(bcc,'EffectOfInterest2_contrasts')||isempty(bcc.EffectOfInterest2_contrasts),bcc.EffectOfInterest2_contrasts={};end
        nc=1;
        spm_data=[];Ec=[];
        EffectOfInterest2_contrasts=bcc.EffectOfInterest2_contrasts;
        bcc.EffectOfInterest2_contrasts={};
        while nc<=size(bcc.EffectOfInterest2_spm,1),
            current_spm=bcc.EffectOfInterest2_spm{min(size(bcc.EffectOfInterest2_spm,1),nc),1};
            [spm_data,SPM,Ec(nc)]=spm_ss_importspm(spm_data,current_spm);
            Cnames={SPM.xCon(:).name};
            if size(bcc.EffectOfInterest2_spm,1)==1
                ic=[];ok=1;for n1=1:length(EffectOfInterest2_contrasts),temp=strmatch(EffectOfInterest2_contrasts{n1},Cnames,'exact');if numel(temp)==1,ic(n1)=temp; elseif EffectOfInterest2_contrasts{n1}(1)=='*', temp=strmatch(EffectOfInterest2_contrasts{n1}(2:end),Cnames,'exact'); if numel(temp)==1, ic(n1)=-temp; else ok=0;break; end; else ok=0;break; end; end
            else
                ic=[];ok=1;for n1=nc:min(nc,length(EffectOfInterest2_contrasts)),temp=strmatch(EffectOfInterest2_contrasts{n1},Cnames,'exact');if numel(temp)==1,ic=temp; elseif EffectOfInterest2_contrasts{n1}(1)=='*', temp=strmatch(EffectOfInterest2_contrasts{n1}(2:end),Cnames,'exact'); if numel(temp)==1, ic=-temp; else ok=0;break; end; else ok=0;break; end; end
            end
            if bcc.askn>1||~ok||isempty(ic),
                if size(bcc.EffectOfInterest2_spm,1)==1
                    str='Select EFFECT OF INTEREST2 contrasts';
                    disp(str);
                    Ic=listdlg('promptstring',str,'selectionmode','multiple','liststring',Cnames,'initialvalue',abs(ic)); %Ic=spm_conman(SPM,'T|F',inf,str,'',0);
                else
                    str=['Select EFFECT OF INTEREST2 contrast #',num2str(nc)];
                    disp(str);
                    Ic=listdlg('promptstring',str,'selectionmode','single','liststring',Cnames,'initialvalue',abs(ic)); %Ic=spm_conman(SPM,'T|F',inf,str,'',0);
                end
            else
                Ic=ic;
            end
            if numel(ic)==numel(Ic) && all(abs(ic)==abs(Ic)), Ic=abs(Ic).*sign(ic); end
            if numel(ic)~=numel(Ic) || any(ic(:)~=Ic(:)), bcc.EffectOfInterest2={}; end
            if any(Ic<0)
                error('This option not implemented in bcc analyses');
            end
            bcc.EffectOfInterest2_contrasts(nc+(1:numel(Ic))-1)=Cnames(Ic); % note: the gui does not allow to manually cross-validate when using this option (i.e. multiple contrasts are interpreted as multiple effects of interest contrasts, not as multiple 'sessions'), use batch scripts instead
            nc=nc+numel(Ic);
        end
    end
end
if bcc.files_selectmanually(2)&&(bcc.askn>1||~isfield(bcc,'EffectOfInterest2_contrasts')||isempty(bcc.EffectOfInterest2_contrasts)),
    for nc=1:numel(bcc.EffectOfInterest2{1}),
        if numel(bcc.EffectOfInterest2_contrasts)<nc,bcc.EffectOfInterest2_contrasts{nc}=['Effect #',num2str(nc)];end
        if bcc.askn,
            if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
            str=['Effect #',num2str(nc),'name?'];
            disp(str);
            bcc.EffectOfInterest2_contrasts{nc}=spm_input(str,posstr,'s',bcc.EffectOfInterest2_contrasts{nc});posstr='+1';
        end
    end
end

if bcc.askn>1||~isfield(bcc,'ExplicitMasking'),
    if ~isfield(bcc,'ExplicitMasking')||isempty(bcc.ExplicitMasking), bcc.ExplicitMasking=''; end
    str='Explicit masking? (this is an optional aditional mask used across all subjects)';
    disp(str);
    sstype=spm_input(str,posstr,'b','none|select',[1,2],1+~isempty(bcc.ExplicitMasking));posstr='+1';
    if sstype==1,bcc.ExplicitMasking='';
    else
        str='Select mask file';
        disp(str);
        bcc.ExplicitMasking=spm_select([0,1],'image',str,{bcc.ExplicitMasking});
    end
end

if bcc.askn>1||~isfield(bcc,'ManualROIs')||isempty(bcc.ManualROIs),
    if ~isfield(bcc,'ManualROIs')||isempty(bcc.ManualROIs), bcc.ManualROIs=''; end
    str='Type of ROI file?';
    disp(str);
    sstype=spm_input(str,posstr,'m','Single ROI file|Subject-specific ROI file',[1,2], sstype); posstr='+1';
    if sstype==1
        str='Select ROI file (ROIs are labeled as integer numbers within this volume)';
        disp(str);
        bcc.ManualROIs=spm_select(1,'image',str,cellstr(bcc.ManualROIs));
    elseif sstype==2
        str='Select ROI files (one file per subject; ROIs are labeled as integer numbers within each volume)';
        disp(str);
        bcc.ManualROIs=cellstr(spm_select(bcc.n,'image',str,cellstr(bcc.ManualROIs)));
    end
end
if ~isfield(bcc,'ManualROIs'), bcc.ManualROIs=''; end
if ~iscell(bcc.ManualROIs)&&~isempty(bcc.ManualROIs), bcc.ManualROIs=cellstr(bcc.ManualROIs); end

if bcc.askn>1||~isfield(bcc,'model')||isempty(bcc.model), 
    if ~isfield(bcc,'model')||isempty(bcc.model), bcc.model=1; end
    if bcc.askn,
        if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
        str='Select model for second-level analysis';
        disp(str);
        bcc.model=spm_input(str,posstr,'m','one-sample t-test|two-sample t-test|multiple regression',[],bcc.model);posstr='+1';
    end
end

if bcc.askn>1||~isfield(bcc,'X')||isempty(bcc.X),
    switch(bcc.model),
        case 1,
            bcc.X=ones(bcc.n,1);
            default_xname={'Group'};
            default_cLcontrast={1};default_cLname={'Group'};
        case 2,
           if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
           if ~isfield(bcc,'X')||isempty(bcc.X), bcc.X=[]; end
            done=0;
            idx1=find(bcc.X(:,1));
            while ~done,
                str='Subjects in group 1?';disp(str);
                idx1=spm_input(str,posstr,'n',mat2str(idx1(:)'),[1,inf],bcc.n);posstr='+1';
                str='Subjects in group 2?';disp(str);
                idx2=spm_input(str,posstr,'n',mat2str(setdiff(1:bcc.n,idx1(:)')),[1,inf],bcc.n);posstr='+1';
                if numel(idx1)+numel(idx2)==bcc.n&&isempty(intersect(idx1,idx2)),done=1;end
            end
            bcc.X=zeros(bcc.n,2);bcc.X(idx1,1)=1;bcc.X(idx2,2)=1;
            default_xname={'Group 1','Group 2'};
            default_cLcontrast={};
        case 3,
            if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
            if ~isfield(bcc,'X')||isempty(bcc.X), bcc.X=[]; end
            str='Regressor matrix?';
            disp(str);
            bcc.X=spm_input(str,posstr,'r',mat2str(bcc.X),[bcc.n,inf]);posstr='+1';
            default_xname=cellstr([repmat('Regressor #',[size(bcc.X,2),1]),num2str((1:size(bcc.X,2))')]);
            default_cLcontrast={};
        otherwise,
            default_xname={};
            default_cLcontrast={};
    end
    if ~isfield(bcc,'Xname')||numel(bcc.Xname)~=size(bcc.X,2), bcc.Xname=default_xname; end
else
    default_xname=repmat({'regressor-name'},[size(bcc.X,2),1]);
    default_cLcontrast={};
end

if bcc.askn>1||~isfield(bcc,'Xname')||numel(bcc.Xname)~=size(bcc.X,2), 
    if ~isfield(bcc,'Xname')||numel(bcc.Xname)~=size(bcc.X,2), bcc.Xname=default_xname; end
    for nc=1:size(bcc.X,2),
        if numel(bcc.Xname)<nc,bcc.Xname{nc}=['Regressor #',num2str(nc)];end
        if bcc.askn,
            if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
            str=['Regressor #',num2str(nc),' name?'];disp(str);
            bcc.Xname{nc}=spm_input(str,posstr,'s',bcc.Xname{nc});posstr='+1';
        end
    end
end

if ~isfield(bcc,'C'), bcc.C=repmat(struct('between',[],'within',[],'name',''),[1,0]); end 
if ~isempty(default_cLcontrast),
    nnc=1;
    neffects=size(bcc.EffectOfInterest1_contrasts,2)*size(bcc.EffectOfInterest2_contrasts,2);
    effect_names=cellfun(@(a,b)[a,'/',b],repmat(bcc.EffectOfInterest1_contrasts(:),1,numel(bcc.EffectOfInterest2_contrasts)),repmat(bcc.EffectOfInterest2_contrasts(:)',numel(bcc.EffectOfInterest1_contrasts),1),'uni',0);
    I=eye(neffects);
    for nc=1:numel(default_cLcontrast),
        if isempty(default_cRcontrast)
            for ne=1:neffects,
                bcc.C(nnc)=struct('between',default_cLcontrast{nc},'within',I(ne,:),'name',[default_cLname{nc},' (',effect_names{ne},')']);
                nnc=nnc+1;
            end
            bcc.C(nnc)=struct('between',default_cLcontrast{nc},'within',I,'name',[default_cLname{nc},' (',strcat(effect_names{:}),')']);
            nnc=nnc+1;
        else
            for ne1=1:numel(default_cRname),
                for ne2=1:numel(default_cRname{ne1}),
                    bcc.C(nnc)=struct('between',default_cLcontrast{nc},'within',cat(2,default_cRcontrast{ne1}(:,:,ne2),zeros(size(default_cRcontrast{ne1},1),neffects-size(default_cRcontrast{ne1},2))),'name',[default_cLname{nc},' (',default_cRname{ne1}{ne2},')']);
                    nnc=nnc+1;
                end
            end
        end
    end
end

% imports contrast files / orthogonalize if necessary
if ~any(bcc.files_selectmanually)&&(~isfield(bcc,'EffectOfInterest2')||isempty(bcc.EffectOfInterest2)||~isfield(bcc,'EffectOfInterest1')||isempty(bcc.EffectOfInterest1)),
    initestimation=1;
    okorth=1;
    bcc.EffectOfInterest1={};
    bcc.EffectOfInterest2={};
    disp('Importing contrast information from SPM.mat files. Please wait...');
    for np=1:size(bcc.EffectOfInterest1_spm,2),%subjects
        spm_data=[];
        Ic1=[];Ec1=[];ok=1;
        for ncontrast=1:size(bcc.EffectOfInterest1_contrasts,2),%contrasts
            current_spm=bcc.EffectOfInterest1_spm{min(size(bcc.EffectOfInterest1_spm,1),ncontrast),np};
            [spm_data,SPM,Ec1(ncontrast)]=spm_ss_importspm(spm_data,current_spm);
            Cnames={SPM.xCon(:).name};
            for ncrossvalid=1:size(bcc.EffectOfInterest1_contrasts,1),%crossvalidation partitions
                temp=strmatch(bcc.EffectOfInterest1_contrasts{ncrossvalid,ncontrast},Cnames,'exact');if numel(temp)~=1,ok=0;else Ic1(ncrossvalid,ncontrast)=temp;end
                if ~ok, error(['the target contrasts are not found inside ',current_spm]); end
            end
        end
        Ic2=[];Ec2=[];ok=1;
        for ncontrast=1:size(bcc.EffectOfInterest2_contrasts,2),%contrasts
            current_spm=bcc.EffectOfInterest2_spm{min(size(bcc.EffectOfInterest2_spm,1),ncontrast),np};
            [spm_data,SPM,Ec2(ncontrast)]=spm_ss_importspm(spm_data,current_spm);
            Cnames={SPM.xCon(:).name};
            for ncrossvalid=1:size(bcc.EffectOfInterest2_contrasts,1),%crossvalidation partitions
                temp=strmatch(bcc.EffectOfInterest2_contrasts{ncrossvalid,ncontrast},Cnames,'exact');if numel(temp)~=1,ok=0;else Ic2(ncrossvalid,ncontrast)=temp;end
                if ~ok, error(['the target contrasts are not found inside ',current_spm]); end
            end
        end
        
        ncross=size(Ic2,1);
        neffects2=size(Ic2,2);
        neffects1=size(Ic1,2);
        if size(Ic1,1)~=ncross, error('Mismatched number of partitions between EFFECT OF INTEREST1 and EFFECT OF INTEREST2 contrasts'); end

        Inewc1={};Inewc2={};automaticcrossvalidation=0;
        % checks orthogonality & create contrasts if necessary, separately for each experiment (SPM.mat file) involved
        for nexp=1:numel(spm_data.SPMfiles),
            iIc1=find(Ec1==nexp);iIc2=find(Ec2==nexp);
            Inewc1{nexp}=Ic1(:,iIc1);Inewc2{nexp}=Ic2(:,iIc2);
            if ~isempty(iIc1)&&~isempty(iIc2),        
                o=spm_SpUtil('ConO',spm_data.SPM{nexp}.xX.X,[spm_data.SPM{nexp}.xCon([Ic1(:,iIc1),abs(Ic2(:,iIc2))]).c]);
                o=permute(reshape(o(1:numel(iIc1),numel(iIc1)+1:end),[ncross,numel(iIc1),ncross,numel(iIc2)]),[2,4,1,3]);
                oo=o(:,:,1:ncross+1:ncross*ncross);
                if ~all(oo(:))&&okorth&&ncross>1,
                    str={'WARNING! You are choosing manual cross-validation, yet not all of the EFFECT OF INTEREST1 and EFFECT OF INTEREST2 contrast pairs selected are orthogonal',...
                        ['pwd = ',bcc.swd],...
                        'Are you sure you want to continue? (the results will likely be invalid)'};
                    fidx=find(~oo);[fneffects,fnconjunction,fncross]=ind2sub(size(oo),fidx);for n1=1:numel(fneffects),disp(['EFFECT OF INTEREST1 ',bcc.EffectOfInterest1_contrasts{(fneffects(n1)-1)*ncross+fncross(n1)},' not orthogonal to EFFECT OF INTEREST2 ',bcc.EffectOfInterest2_contrasts{(fnconjunction(n1)-1)*ncross+fncross(n1)}]); end
                    disp(char(str));
                    if bcc.askn,if spm_input(str,posstr,'bd','stop|continue',[1,0],0),return;end;posstr='+1';end
                    disp('...continuing anyway');
                    okorth=0;
                elseif ~all(oo(:)),
                    if okorth,
                        str={'EFFECT OF INTEREST1 and EFFECT OF INTEREST2 contrast pairs are not orthogonal',...
                            ['pwd = ',bcc.swd],...
                            'New contrasts (broken down by sessions) will be created now if they do not already exist'};
                        disp(char(str));
                        if bcc.askn,if spm_input(str,posstr,'bd','stop|continue',[1,0],0),return;end;posstr='+1';end
                    end
                    [spm_data.SPM{nexp},IcCV1,sess1]=spm_ss_crossvalidate_sessions(spm_data.SPM{nexp},Ic1(:,iIc1),0,'skip');
                    [spm_data.SPM{nexp},IcCV2,sess2]=spm_ss_crossvalidate_sessions(spm_data.SPM{nexp},Ic2(:,iIc2),0,'skip');
                    [nill,idxvalid1,idxvalid2]=intersect(sess1,sess2);
                    if isempty(nill), error(['Subject #',num2str(np),' (',spm_data.SPMfiles{nexp},') not possible to cross-validate.']);
                    else disp(['Subject ',num2str(np),' cross-validation across sessions ',num2str(nill(:)'),' in file ',spm_data.SPMfiles{nexp}]); end
                    Inewc1{nexp}=IcCV1(:,idxvalid1,1)';%Ic1=Ic1(:);
                    Inewc2{nexp}=(diag(sign(Ic2(:,iIc2)))*IcCV2(:,idxvalid2,2))';%Ic2=Ic2(:);
                    okorth=0;
                    automaticcrossvalidation=1;
                end
            end
        end
        if automaticcrossvalidation,
            % combines new cross-validated contrasts (possibly across multiple experiments)
            ic1=Ic1;ic2=Ic2;
            nexps=numel(spm_data.SPMfiles);
            nexpcross=zeros(1,nexps);for nexp=1:nexps,nexpcross(nexp)=size(Inewc1{nexp},1);end
            Ic1=zeros([size(ic1,2),nexpcross]);
            Ic2=zeros([size(ic2,2),nexpcross]);
            for nexp=1:nexps,
                iIc1=find(Ec1==nexp);iIc2=find(Ec2==nexp);
                if ~isempty(iIc1), Ic1(iIc1,:,:)=repmat(Inewc1{nexp}',[1,1,size(Ic1(:,:,:),3)]); end
                if ~isempty(iIc2), Ic2(iIc2,:,:)=repmat(Inewc2{nexp}',[1,1,size(Ic2(:,:,:),3)]); end
                Ic1=permute(Ic1,[1,3:nexps+1,2]);
                Ic2=permute(Ic2,[1,3:nexps+1,2]);
            end
            Ic1=Ic1(:,:)';
            Ic2=Ic2(:,:)';
        end
        ncross=size(Ic2,1);
        neffects2=size(Ic2,2);
        neffects1=size(Ic1,2);
        for nc=1:ncross,
            for ne=1:neffects1,
                bcc.EffectOfInterest1{np}{ne,nc}=fullfile(fileparts(bcc.EffectOfInterest1_spm{min(size(bcc.EffectOfInterest1_spm,1),ne),np}),['con_',num2str(Ic1((ne-1)*ncross+nc),'%04d'),'.nii']);
                if isempty(dir(bcc.EffectOfInterest1{np}{ne,nc})), fprintf('warning: file %s not found. Trying .img contrast files\n',bcc.EffectOfInterest1{np}{ne,nc}); bcc.EffectOfInterest1{np}{ne,nc}=regexprep(bcc.EffectOfInterest1{np}{ne,nc},'\.nii$','.img'); end
            end
            for ne=1:neffects2,
                bcc.EffectOfInterest2{np}{ne,nc}=fullfile(fileparts(bcc.EffectOfInterest2_spm{min(size(bcc.EffectOfInterest2_spm,1),ne),np}),['con_',num2str(Ic2((ne-1)*ncross+nc),'%04d'),'.nii']);
                if isempty(dir(bcc.EffectOfInterest2{np}{ne,nc})), fprintf('warning: file %s not found. Trying .img contrast files\n',bcc.EffectOfInterest2{np}{ne,nc}); bcc.EffectOfInterest2{np}{ne,nc}=regexprep(bcc.EffectOfInterest2{np}{ne,nc},'\.nii$','.img'); end
            end
        end
        fprintf('.');
    end
    fprintf('\n');
end

if initestimation, bcc.estimate={}; bcc.evaluate={}; end
bcc.ask=[];bcc.askn=[];

objname=findobj('tag','spm_bcc');if numel(objname)==1,objdata.files_spm_bcc=fullfile(bcc.swd,['SPM_bcc.mat']);set(objname,'userdata',objdata);end;
save(fullfile(bcc.swd,['SPM_bcc.mat']),'bcc');
disp(['Analysis file saved: ',fullfile(bcc.swd,['SPM_bcc.mat'])]);

end
    
    
