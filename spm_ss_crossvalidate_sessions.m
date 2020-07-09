function [SPM,allIc,validsessions,nsess,txt]=spm_ss_crossvalidate_sessions(SPM,Ic,overwrite,mask_type,mask_p)
% SPM_SS_CROSSVALIDATE_SESSIONS
% breaks down selected contrasts by sessions
% (one possible way to orthogonalize contrasts)
%
% For each initial contrast, SPM_SS_CROSSVALIDATE_SESSIONS creates and
% estimates several new contrasts:
%  a) SESSION{n}_contrastname (the original contrast estimated only using
%  data from session n)
%  b) ORTH_TO_SESSION{n}_contrastname (the original contrast estimated
%  using data from all sessions except n)
%
% Any contrast of the form SESSION{n}_* will in this way always be
% orthogonal to any contrast of the form ORTH_TO_SESSION{n}*
%

allIc=[];
validsessions=[]; 
nsess=[];
txt={};

if nargin<3, overwrite=1; end
if nargin<4, mask_type=[]; end
if nargin<5, mask_p=[]; end

if nargin<1||isempty(SPM)||ischar(SPM)||iscell(SPM),
    if nargin>0 && ~isempty(SPM), 
        if ischar(SPM),P=cellstr(SPM);else P=SPM; end
    else
        str='Select first-level SPM.mat(s) (one per subject)';
        disp(str);
        Pdefault={''};objname=findobj('tag','spm_ss');if numel(objname)==1,objdata=get(objname,'userdata');if isfield(objdata,'files_spm'),Pdefault=objdata.files_spm;end;end;
        P=cellstr(spm_select(inf,'^SPM\.mat$',str,Pdefault));
        if numel(objname)==1&&~isempty(P),objdata.files_spm=P;set(objname,'userdata',objdata);end;
    end
    if nargin>1&&(ischar(Ic)||iscell(Ic)),
        if ischar(Ic),ContrastNames=cellstr(Ic);
        else ContrastNames=Ic; end
    else
        ContrastNames={};
    end
    txt={};
    nsess=[];
    for np=1:length(P),
        load(P{np},'SPM');
        SPM.swd=fileparts(P{np});
        
        if ~isempty(ContrastNames),
            Cnames={SPM.xCon(:).name};
            ic=[];ok=1;for n1=1:length(ContrastNames),temp=strmatch(ContrastNames{n1},Cnames,'exact');if numel(temp)~=1,ok=0;break;else ic(n1)=temp;end;end
            if ~ok, error(['the target contrasts are not found inside ',P{np}]); end
            Ic=ic;
        end
        if isempty(Ic),
            Ic=spm_conman(SPM,'T|F',inf,'Select contrast(s) to cross-validate across sessions','',0);
        end
        if isempty(ContrastNames),ContrastNames={SPM.xCon(Ic).name}; end
        
        [nill,nill,nill,nsessions,txtn]=spm_ss_crossvalidate_sessions(SPM,Ic);
        nsess=cat(1,nsess,nsessions);
        txt=cat(1,txt,txtn{:});
    end
    if ~nargout, 
        %disp(strvcat(txt{:})); 
        disp('Number of sessions per subject cross-validated:'); 
        disp(nsess(:)');
    end
    return;
end

cwd=pwd;

if nargin<2||isempty(Ic),
    Ic=spm_conman(SPM,'T|F',inf,'Select contrast(s) to cross-validate across sessions','',0);
end

% build a list of session indices per regressor

if isfield(SPM,'Sess'),
    names=SPM.xX.name;
    sessions=zeros(numel(names),1);
    for i=1:numel(SPM.Sess),
        sessions(SPM.Sess(i).col)=i;
    end
else
    names=SPM.xX.name;
    sessions=zeros(numel(names),1);
    for i=1:numel(names),
        [a,v,w] = regexp(names{i},'Sn\(\d+\)','match','start','end'); % removes Sn(???) from name
        if ~isempty(a), sessions(i)=str2num(a{1}(4:end-1)); end
    end
end
nsess=max(sessions);
%disp(['Session information: ',num2str(sessions(:)')]);
if nsess<=1,disp(['Sorry, no multiple session information was found in this SPM file. pwd: ',SPM.swd]); return;end

% adds new Cross-validation contrasts for the selected contrasts of interest
ContrastNames={SPM.xCon(:).name};
okdata=ones(1,nsess);
for nIc=1:numel(Ic),
    for nses=1:nsess,%first checks whether we will have valid SESSION/ORTH_TO_SESSION contrasts for each session for all contrasts
        C=SPM.xCon(Ic(nIc)).c;
        C(sessions~=nses,:)=0;
        okdata(nses)=okdata(nses)&all(any(C~=0,1),2);
        C=SPM.xCon(Ic(nIc)).c;
        C(sessions==nses,:)=0;
        okdata(nses)=okdata(nses)&all(any(C~=0,1),2);
    end
end
validsessions=find(okdata);
allIc=nan+zeros(numel(Ic),numel(validsessions),2);
newIc=[];
for nIc=1:numel(Ic),
    for nvalid=1:numel(validsessions),%find and estimate if necessary the SESSION/ORTH_TO_SESSION contrasts for each session
        nses=validsessions(nvalid);
        cname=['SESSION',num2str(nses,'%02d'),'_',SPM.xCon(Ic(nIc)).name];
        idxmatch=strmatch(cname,ContrastNames,'exact');
        if isempty(idxmatch),
            C=SPM.xCon(Ic(nIc)).c;
            C(sessions~=nses,:)=0;
            for nc=1:size(C,2),nC=max(sum(C(C(:,nc)>0,nc)),sum(-C(C(:,nc)<0,nc)));C(:,nc)=C(:,nc)/nC;end
            SPM.xCon(end+1)=spm_FcUtil('Set',cname,SPM.xCon(Ic(nIc)).STAT,'c',C,SPM.xX.xKXs);
            newIc(end+1)=length(SPM.xCon);
            allIc(nIc,nvalid,1)=length(SPM.xCon);
            txt{end+1}=['Contrast #',num2str(newIc(end),'%04d'),': ',cname];
        else
            if overwrite,%overwrite existing contrast
                C=SPM.xCon(Ic(nIc)).c;
                C(sessions~=nses,:)=0;
                for nc=1:size(C,2),nC=max(sum(C(C(:,nc)>0,nc)),sum(-C(C(:,nc)<0,nc)));C(:,nc)=C(:,nc)/nC;end
                SPM.xCon(idxmatch(1))=spm_FcUtil('Set',cname,SPM.xCon(Ic(nIc)).STAT,'c',C,SPM.xX.xKXs);
                newIc(end+1)=idxmatch(1);
            end
            allIc(nIc,nvalid,1)=idxmatch(1);
            txt{end+1}=['Contrast #',num2str(idxmatch(1),'%04d'),': ',cname];
        end
    end
    for nvalid=1:numel(validsessions),%find and estimate if necessary the SESSION/ORTH_TO_SESSION contrasts for each session
        nses=validsessions(nvalid);
        cname=['ORTH_TO_SESSION',num2str(nses,'%02d'),'_',SPM.xCon(Ic(nIc)).name];
        idxmatch=strmatch(cname,ContrastNames,'exact');
        if isempty(idxmatch),
            C=SPM.xCon(Ic(nIc)).c;
            C(sessions==nses,:)=0;
            for nc=1:size(C,2),nC=max(sum(C(C(:,nc)>0,nc)),sum(-C(C(:,nc)<0,nc)));C(:,nc)=C(:,nc)/nC;end
            SPM.xCon(end+1)=spm_FcUtil('Set',cname,SPM.xCon(Ic(nIc)).STAT,'c',C,SPM.xX.xKXs);
            newIc(end+1)=length(SPM.xCon);
            allIc(nIc,nvalid,2)=length(SPM.xCon);
            txt{end+1}=['Contrast #',num2str(newIc(end),'%04d'),': ',cname];
        else
            if overwrite,%overwrite existing contrast
                C=SPM.xCon(Ic(nIc)).c;
                C(sessions==nses,:)=0;
                for nc=1:size(C,2),nC=max(sum(C(C(:,nc)>0,nc)),sum(-C(C(:,nc)<0,nc)));C(:,nc)=C(:,nc)/nC;end
                SPM.xCon(idxmatch(1))=spm_FcUtil('Set',cname,SPM.xCon(Ic(nIc)).STAT,'c',C,SPM.xX.xKXs);
                newIc(end+1)=idxmatch(1);
            end
            allIc(nIc,nvalid,2)=idxmatch(1);
            txt{end+1}=['Contrast #',num2str(idxmatch(1),'%04d'),': ',cname];
        end
    end
end

if ~isempty(newIc),
    % estimates new contrasts
    cwd=pwd;
    reduced=0;
    if isfield(SPM,'spm_ss_SPMfile'),
        reduced=1;
        spm_ss_SPMfile=SPM.spm_ss_SPMfile;
        spm_xCon=SPM.xCon;
        load(SPM.spm_ss_SPMfile,'SPM');
        SPM.swd=fileparts(spm_ss_SPMfile);
        SPM.xCon(newIc)=spm_xCon(newIc);
    else 
        spm_ss_SPMfile=fullfile(SPM.swd,'SPM.mat');
    end
    cd(SPM.swd);
    SPM=spm_contrasts(SPM,newIc);
    save(spm_ss_SPMfile,'SPM','-append');
    cd(cwd);
    if reduced
        [nill,SPM]=spm_ss_importspm(SPM,spm_ss_SPMfile);
    end
end
if ~isempty(allIc),
    % creates threhsolded volumes
    spm_ss_createlocalizermask({SPM},allIc(:)',[],overwrite,mask_type,mask_p);
end

if ~nargout, disp(strvcat(txt{:})); end

