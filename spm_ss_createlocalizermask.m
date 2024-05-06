function maskfilename=spm_ss_createlocalizermask(SPM,Ic,Ec,overwrite,thr_type,thr,conjunction_type,options,optionROIs)
% SPM_SS_CREATELOCALIZERMASK
% thresholds 1st-level contrast and create binary masks that can be used as
% localizer images.
%

if nargin<1,SPM={};end
if nargin<2,Ic=[];end
if nargin<3,Ec=[];end
if nargin<4,overwrite=1;end
if nargin<5,thr_type={};end
if nargin<6,thr=[];end
if nargin<7||isempty(conjunction_type),conjunction_type='and'; end
if nargin<8,options='';end
if nargin<9,optionROIs='';end
if ~isempty(thr_type)&&ischar(thr_type)&&isempty(strmatch(thr_type,{'FDR','FWE','none','percentile-whole-brain','percentile-ROI-level','Nvoxels-whole-brain','Nvoxels-ROI-level'},'exact')),return;end
maskfilename={};
if ~iscell(thr_type),thr_type=cellstr(thr_type); end
if ~isempty(optionROIs), optionROIs=spm_vol(optionROIs); end
filemask_options='0002';
try
    [ok,msg]=system('umask');
    if ~ok,
        msg1=base2dec(regexprep(msg,'[^012345678]',''),8);
        msg2=base2dec(filemask_options,8);
        if bitand(msg1,intmax-msg2)
            fprintf('warning: unexpected file mode creation mask. Mask set to %s. Use system-command ''umask %s'' to revert this change\n',dec2base(bitand(msg1,msg2),8),dec2base(msg1,8));
            [nill,msg]=system(sprintf('umask %s',dec2base(bitand(msg1,msg2),8)));
        end
    end
end

posstr=1;
if isempty(SPM),
    str='Select first-level SPM.mat(s) (one per subject)';
    disp(str);
    Pdefault={''};objname=findobj('tag','spm_ss');if numel(objname)==1,objdata=get(objname,'userdata');if isfield(objdata,'files_spm'),Pdefault=objdata.files_spm;end;end;
    P=cellstr(spm_select(inf,'^SPM\.mat$',str,Pdefault));
    if numel(objname)==1&&~isempty(P),objdata.files_spm=P;set(objname,'userdata',objdata);end;
    
    for np=1:numel(P),
        load(P{np},'SPM');
        SPM.swd=fileparts(P{np});
        maskfilename{np}=spm_ss_createlocalizermask({SPM},Ic,[],overwrite,thr_type,thr,conjunction_type);
    end
    return;
end

if ~isempty(Ic)&&(ischar(Ic)||iscell(Ic)),
    if ischar(Ic),ContrastNames=cellstr(Ic);
    else ContrastNames=Ic; end
else
    ContrastNames={};
end

if isempty(Ec), if ~isempty(ContrastNames), Ec=ones(size(ContrastNames)); else Ec=ones(size(Ic)); end; end

if ~isempty(ContrastNames),
    ic=[];ok=1;
    for n1=1:numel(ContrastNames),
        Cnames={SPM{Ec(n1)}.xCon(:).name};
        temp=strmatch(ContrastNames{n1},Cnames,'exact');
        if numel(temp)~=1&&~isempty(strmatch('-not',ContrastNames{n1})),temp=-strmatch(ContrastNames{n1}(5:end),Cnames,'exact');end;
        if numel(temp)~=1,ok=0;break;else ic(n1)=temp;end;
    end
    if ~ok, error(['the target contrasts are not found inside ',SPM{Ec(n1)}.swd]); end
    Ic=reshape(ic,size(ContrastNames));
end
if isempty(Ic),
    ic={};ec={};
    for np=1:numel(SPM),
        str=['Select LOCALIZER contrast(s) from ',SPM{np}.swd];
        disp(str);
        ic{np}=spm_conman(SPM{np},'T|F',inf,str,'',0);
        ic{np}=ic{np}(:)';
        ec{np}=np+zeros(1,numel(ic{np}));
    end
    Ic=cat(2,ic{:});
    Ec=cat(2,ec{:});
    if length(Ic)>1,
        if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
        str='Create a separate localizer volume per contrast|Create a single localizer volume from the conjunction of the selected contrasts';
        disp(str);
        conjunction=spm_input('how to handle multiple contrasts?',posstr,'m',str,[0,1],1);posstr='+1';
        if conjunction, 
            Ic=Ic(:)'; 
            Ec=Ec(:)';
            str='Inclusive (1) or exclusive (0) mask for each contrast?';
            disp(str);
            signIc=spm_input(str,posstr,'r',num2str(ones(1,numel(Ic))),[1,numel(Ic)],[0,1]);posstr='+1';
            Ic(signIc<.5)=-Ic(signIc<.5);
        else Ic=Ic(:); Ec=Ec(:); end
    end
end
%if isempty(ContrastNames),ContrastNames={SPM.xCon(abs(Ic)).name}; end
if isempty(thr_type),
    if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
    str='FDR|FWE|none|percentile-whole-brain|percentile-ROI-level|Nvoxels-whole-brain|Nvoxels-ROI-level';
    thr_type={};
    if numel(thr)>1
        for nic2=1:max(1,numel(thr)),
            thr_type{nic2}=spm_input(['contrast #',num2str(nic2),': p value adjustment to control'],posstr,'b',str,[],1);posstr='+1';
        end
    else
        thr_type={spm_input('p value adjustment to control',posstr,'b',str,[],1)};posstr='+1';
    end
end
if isempty(thr),
    if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
    thr=[];
    if numel(thr_type)>1
        for nic2=1:numel(thr_type)
            if strcmpi(thr_type{nic2},'none'),defp=.001; elseif strncmp(thr_type{nic2},'Nvoxels',7), defp=5; else defp=.05;end
            if strcmpi(thr_type{nic2},'automatic'), thr(nic2)=nan;
            else thr(nic2)=spm_input(['contrast #',num2str(nic2),': p value (',thr_type{nic2},')'],posstr,'r',num2str(defp),1,[0,1]);posstr='+1'; end
        end
    else
        if strcmpi(thr_type{1},'none'),defp=.001; elseif strncmp(thr_type{1},'Nvoxels',7), defp=5; else defp=.05;end
        if strcmpi(thr_type{1},'automatic'), thr=nan;
        else thr=spm_input(['p value (',thr_type{1},')'],posstr,'r',num2str(defp),1,[0,1]);posstr='+1'; end
    end
end
if isempty(optionROIs)&&(any(strcmp(thr_type,'percentile-ROI-level'))||any(strcmp(thr_type,'Nvoxels-ROI-level'))),
    str='Select ROI file (ROIs are labeled as integer numbers within this volume)';
    disp(str);
    optionROIs=spm_select(1,'image',str,{''});
    if ~isempty(optionROIs), optionROIs=spm_vol(optionROIs); end
end

if numel(thr)<size(Ic,2),thr=repmat(thr,[1,size(Ic,2)]);end
if numel(thr_type)<size(Ic,2), thr_type=repmat(thr_type,[1,size(Ic,2)]); end
if size(Ec,1)<size(Ic,1), Ec=Ec(min(size(Ec,1),1:size(Ic,1)),:); end

symbols={'','not'};signIc=1+(sign(Ic)<0);
if strcmp(options,'noconjunction')
    for nic1=1:size(Ic,1), % computes one separate thresholded volume per row and column of Ic (multiple columns are treated as separate contrasts)
        for nic2=1:size(Ic,2),
            filename='locT_';
            filename=[filename,symbols{signIc(nic1,nic2)},num2str(abs(Ic(nic1,nic2)),'%04d'),'_',thr_type{nic2},num2str(thr(nic2))];
            if ~isempty(optionROIs)&&(strcmpi(thr_type{nic2},'percentile-ROI-level')|strcmpi(thr_type{nic2},'Nvoxels-ROI-level')) % note: avoid same localizer file for different optionROIs files
                filename=[filename,'_',char(mlreportgen.utils.hash(fileread(optionROIs(1).fname)))]; % note: requires Matlab 18b or above
                parcelfile={optionROIs.fname};
                %[nill,fname,nill]=fileparts(optionROIs(1).fname);
                %filename=[filename,'_',fname]; 
            else parcelfile={};
            end
            filename=[filename,'.img'];
            maskfilename{nic1,nic2}=fullfile(SPM{Ec(nic1,nic2)}.swd,filename);
            if overwrite||isempty(dir(maskfilename{nic1,nic2})),
                Z=1;
                U={};
                if ~isempty(strmatch(thr_type{nic2},{'FDR','FWE','none','percentile-whole-brain','percentile-ROI-level','Nvoxels-whole-brain','Nvoxels-ROI-level'},'exact')),
                    %         try,
                    a=spm_vol(fullfile(SPM{Ec(nic1,nic2)}.swd,SPM{Ec(nic1,nic2)}.xCon(abs(Ic(nic1,nic2))).Vspm.fname));
                    b=spm_read_vols(a);
                    idx=find(~isnan(b)&b~=0);
                    dof=[SPM{Ec(nic1,nic2)}.xCon(abs(Ic(nic1,nic2))).eidf,SPM{Ec(nic1,nic2)}.xX.erdf];
                    STAT=SPM{Ec(nic1,nic2)}.xCon(abs(Ic(nic1,nic2))).STAT;
                    R=SPM{Ec(nic1,nic2)}.xVol.R;
                    S=SPM{Ec(nic1,nic2)}.xVol.S;
                    n=1;
                    switch(thr_type{nic2}),
                        case 'FWE',
                            u=spm_uc(thr(nic2),dof,STAT,R,n,S);
                            Y=double(b>u);
                        case 'FDR',
                            Y=nan+zeros(size(b));
                            switch(STAT),
                                case 'Z',Y(idx)=1-spm_Ncdf(b(idx));
                                case 'T',Y(idx)=1-spm_Tcdf(b(idx),dof(2));
                                case 'X',Y(idx)=1-spm_Xcdf(b(idx),dof(2));
                                case 'F',Y(idx)=1-spm_Fcdf(b(idx),dof);
                                otherwise, error('null');
                            end
                            Y(:)=spm_ss_fdr(Y(:));
                            Y=double(Y<thr(nic2));
                        case 'none'
                            u=spm_u(thr(nic2),dof,STAT);
                            Y=double(b>u);
                        case 'percentile-whole-brain'
                            Y=zeros(size(b));
                            tN=b(idx);
                            tN(tN==0|isnan(tN))=-inf;
                            tNthreshold=sort(tN);
                            tNthreshold=tNthreshold(max(1,min(numel(tN), ceil(numel(tN)*(1-thr(nic2))))));
                            Y(idx)=double(tN>tNthreshold);
                        case 'percentile-ROI-level'
                            if isempty(optionROIs), error('no ROI file provided for percentile-ROI-level option'); end
                            clear XYZ;
                            [XYZ{1},XYZ{2},XYZ{3}]=ndgrid(1:a(1).dim(1),1:a(1).dim(2),1:a(1).dim(3));
                            XYZ=reshape(cat(4,XYZ{:}),[],3)';XYZ=a(1).mat*cat(1,XYZ,ones(1,size(XYZ,2)));
                            frois=round(spm_get_data(optionROIs,pinv(optionROIs.mat)*XYZ));
                            nrois=max(frois(:));
                            Y=zeros(size(b));
                            for nroi=1:nrois,
                                idxroi=find(frois==nroi);
                                if ~isempty(idxroi)
                                    tN=b(idxroi);
                                    tN(tN==0|isnan(tN))=-inf;
                                    tNthreshold=sort(tN);
                                    tNthreshold=tNthreshold(max(1,min(numel(tN), ceil(numel(tN)*(1-thr(nic2))))));
                                    Y(idxroi)=double(tN>tNthreshold);
                                end
                            end
                        case 'Nvoxels-whole-brain'
                            Y=zeros(size(b));
                            tN=b(idx);
                            tN(tN==0|isnan(tN))=-inf;
                            [tNthreshold,tNidx]=sort(tN);
                            tNidx(tNidx)=1:numel(tNidx);
                            Y(idx)=double(tNidx>numel(tN)-thr(nic2)&~isinf(tN));
                        case 'Nvoxels-ROI-level'
                            if isempty(optionROIs), error('no ROI file provided for Nvoxels-ROI-level option'); end
                            clear XYZ;
                            [XYZ{1},XYZ{2},XYZ{3}]=ndgrid(1:a(1).dim(1),1:a(1).dim(2),1:a(1).dim(3));
                            XYZ=reshape(cat(4,XYZ{:}),[],3)';XYZ=a(1).mat*cat(1,XYZ,ones(1,size(XYZ,2)));
                            frois=round(spm_get_data(optionROIs,pinv(optionROIs.mat)*XYZ));
                            nrois=max(frois(:));
                            Y=zeros(size(b));
                            for nroi=1:nrois,
                                idxroi=find(frois==nroi);
                                if ~isempty(idxroi)
                                    tN=b(idxroi);
                                    tN(tN==0|isnan(tN))=-inf;
                                    [tNthreshold,tNidx]=sort(tN);
                                    tNidx(tNidx)=1:numel(tNidx);
                                    Y(idxroi)=double(tNidx>numel(tN)-thr(nic2)&~isinf(tN));
                                end
                            end
                    end
                    U{nic2}=sprintf('p=%f(%s),u=%f,STAT=%s,dof=[%f,%f],R=[%f,%f,%f,%f],S=%f,n=%d,I=%d;',thr(nic2),thr_type{nic2},mean([max(b(~Y)),min(b(Y>0))]),STAT,dof(1),dof(2),R(1),R(2),R(3),R(4),S,n,sign(Ic(nic1,nic2))>0); %p=spm_P_RF(1,0,u,dof,STAT,R,n);
                    
                    if sign(Ic(nic1,nic2))<0, Y=~Y; end % exclusion mask
                    if ~any(Y(:)>0), disp(['Warning! Output localizer volume #',num2str(nic1), ', contrast name ',symbols{signIc(nic1,nic2)},SPM{Ec(nic1,nic2)}.xCon(abs(Ic(nic1,nic2))).name,',  in contrast file ',fullfile(SPM{Ec(nic1,nic2)}.swd,filename),', contains no supra-threshold voxels']);
                    else disp([num2str(sum(Y(:)>0)),' voxels in output localizer volume #',num2str(nic1), ', contrast name ',symbols{signIc(nic1,nic2)},SPM{Ec(nic1,nic2)}.xCon(abs(Ic(nic1,nic2))).name,',  in contrast file ',fullfile(SPM{Ec(nic1,nic2)}.swd,filename)]); end
                    % computes conjunction across rows of Ic
                    Z=Z&Y;
                else
                    error(['Incorrect threshold type option',thr_type{nic2}]);
                end
                if numel(Z)>1
                    if size(Ic,2)>1,
                        if ~any(Z(:)>0), disp(['Warning! Output localizer volume #',num2str(nic1), ',  in contrast file ',maskfilename{nic1,nic2},', contains no supra-threshold voxels']);
                        else disp([num2str(sum(Z(:)>0)),' voxels in output localizer volume #',num2str(nic1), ',  in contrast file ',maskfilename{nic1,nic2}]); end
                    end
                    cwd0=pwd;
                    cd(fileparts(maskfilename{nic1,nic2}));
                    Vo=struct(...
                        'fname',    filename,...
                        'dim',      a.dim,...
                        'dt',       [spm_type('uint8') spm_platform('bigend')],...
                        'mat',      a.mat,...
                        'pinfo',    [1;0;0],...
                        'descrip',  sprintf('SPM_SS LOCALIZER{%s}',cat(2,U{:})));
                    Vo=spm_write_vol(Vo,Z);
                    spm_jsonwrite([regexprep(filename,'\.nii$|\.img$',''),'.json'],struct('threshold_value',thr(nic2),'threshold_type',{thr_type(nic2)},'threshold_roifile',{parcelfile},'output_nvoxels',nnz(Z),'full_filename',filename));
                    cd(cwd0);
                end
            end
        end
    end
elseif strcmp(conjunction_type,'and')|strcmp(conjunction_type,'or') % CONJUNCTION AND/OR
    for nic1=1:size(Ic,1), % computes one separate thresholded volume per row of Ic (multiple columns are treated as a conjunction)
        filename='locT_conjunction_';
        for nic2=1:size(Ic,2),
            if Ec(nic1,nic2)~=Ec(nic1,1), temp=SPM{Ec(nic1,nic2)}.swd; temp(temp==filesep)='_';filename=[filename,temp]; end
            filename=[filename,symbols{signIc(nic1,nic2)},num2str(abs(Ic(nic1,nic2)),'%04d'),'_',thr_type{nic2},num2str(thr(nic2))];
            if nic2<size(Ic,2), filename=[filename,'_']; if ~strcmp(conjunction_type,'and'),filename=[filename,conjunction_type,'_']; end; end
        end
        if ~isempty(optionROIs)&&(any(strcmp(thr_type,'percentile-ROI-level'))||any(strcmp(thr_type,'Nvoxels-ROI-level'))) % note: avoid same localizer file for different optionROIs files
            filename=[filename,'_',char(mlreportgen.utils.hash(fileread(optionROIs(1).fname)))]; % note: requires Matlab 18b or above
            parcelfile={optionROIs.fname};
            %[nill,fname,nill]=fileparts(optionROIs(1).fname);
            %filename=[filename,'_',fname];
        else parcelfile={};
        end
        filename=[filename,'.img'];
        maskfilename{nic1}=fullfile(SPM{Ec(nic1,1)}.swd,filename);
        if overwrite||isempty(dir(maskfilename{nic1})),
            switch(conjunction_type)
                case 'and',     Z=1;
                case 'or',      Z=0;
                otherwise,      error('unrecognized conjunction type %s',conjuntion_type);
            end
            U={};
            for nic2=1:size(Ic,2),
                if ~isempty(strmatch(thr_type{nic2},{'FDR','FWE','none','percentile-whole-brain','percentile-ROI-level','Nvoxels-whole-brain','Nvoxels-ROI-level'},'exact')),
                    %         try,
                    a=spm_vol(fullfile(SPM{Ec(nic1,nic2)}.swd,SPM{Ec(nic1,nic2)}.xCon(abs(Ic(nic1,nic2))).Vspm.fname));
                    b=spm_read_vols(a);
                    idx=find(~isnan(b)&b~=0);
                    dof=[SPM{Ec(nic1,nic2)}.xCon(abs(Ic(nic1,nic2))).eidf,SPM{Ec(nic1,nic2)}.xX.erdf];
                    STAT=SPM{Ec(nic1,nic2)}.xCon(abs(Ic(nic1,nic2))).STAT;
                    R=SPM{Ec(nic1,nic2)}.xVol.R;
                    S=SPM{Ec(nic1,nic2)}.xVol.S;
                    n=1;
                    switch(thr_type{nic2}),
                        case 'FWE',
                            u=spm_uc(thr(nic2),dof,STAT,R,n,S);
                            Y=double(b>u);
                        case 'FDR',
                            Y=nan+zeros(size(b));
                            switch(STAT),
                                case 'Z',Y(idx)=1-spm_Ncdf(b(idx));
                                case 'T',Y(idx)=1-spm_Tcdf(b(idx),dof(2));
                                case 'X',Y(idx)=1-spm_Xcdf(b(idx),dof(2));
                                case 'F',Y(idx)=1-spm_Fcdf(b(idx),dof);
                                otherwise, error('null');
                            end
                            Y(:)=spm_ss_fdr(Y(:));
                            Y=double(Y<thr(nic2));
                        case 'none'
                            u=spm_u(thr(nic2),dof,STAT);
                            Y=double(b>u);
                        case 'percentile-whole-brain'
                            Y=zeros(size(b));
                            tN=b(idx);
                            tN(tN==0|isnan(tN))=-inf;
                            tNthreshold=sort(tN);
                            tNthreshold=tNthreshold(max(1,min(numel(tN), ceil(numel(tN)*(1-thr(nic2))))));
                            Y(idx)=double(tN>tNthreshold);
                        case 'percentile-ROI-level'
                            if isempty(optionROIs), error('no ROI file provided for percentile-ROI-level option'); end
                            clear XYZ;
                            [XYZ{1},XYZ{2},XYZ{3}]=ndgrid(1:a(1).dim(1),1:a(1).dim(2),1:a(1).dim(3));
                            XYZ=reshape(cat(4,XYZ{:}),[],3)';XYZ=a(1).mat*cat(1,XYZ,ones(1,size(XYZ,2)));
                            frois=round(spm_get_data(optionROIs,pinv(optionROIs.mat)*XYZ));
                            nrois=max(frois(:));
                            Y=zeros(size(b));
                            for nroi=1:nrois,
                                idxroi=find(frois==nroi);
                                if ~isempty(idxroi)
                                    tN=b(idxroi);
                                    tN(tN==0|isnan(tN))=-inf;
                                    tNthreshold=sort(tN);
                                    tNthreshold=tNthreshold(max(1,min(numel(tN), ceil(numel(tN)*(1-thr(nic2))))));
                                    Y(idxroi)=double(tN>tNthreshold);
                                end
                            end
                        case 'Nvoxels-whole-brain'
                            Y=zeros(size(b));
                            tN=b(idx);
                            tN(tN==0|isnan(tN))=-inf;
                            [tNthreshold,tNidx]=sort(tN);
                            tNidx(tNidx)=1:numel(tNidx);
                            Y(idx)=double(tNidx>numel(tN)-thr(nic2)&~isinf(tN));
                        case 'Nvoxels-ROI-level'
                            if isempty(optionROIs), error('no ROI file provided for Nvoxels-ROI-level option'); end
                            clear XYZ;
                            [XYZ{1},XYZ{2},XYZ{3}]=ndgrid(1:a(1).dim(1),1:a(1).dim(2),1:a(1).dim(3));
                            XYZ=reshape(cat(4,XYZ{:}),[],3)';XYZ=a(1).mat*cat(1,XYZ,ones(1,size(XYZ,2)));
                            frois=round(spm_get_data(optionROIs,pinv(optionROIs.mat)*XYZ));
                            nrois=max(frois(:));
                            Y=zeros(size(b));
                            for nroi=1:nrois,
                                idxroi=find(frois==nroi);
                                if ~isempty(idxroi)
                                    tN=b(idxroi);
                                    tN(tN==0|isnan(tN))=-inf;
                                    [tNthreshold,tNidx]=sort(tN);
                                    tNidx(tNidx)=1:numel(tNidx);
                                    Y(idxroi)=double(tNidx>numel(tN)-thr(nic2)&~isinf(tN));
                                end
                            end
                    end
                    U{nic2}=sprintf('p=%f(%s),u=%f,STAT=%s,dof=[%f,%f],R=[%f,%f,%f,%f],S=%f,n=%d,I=%d;',thr(nic2),thr_type{nic2},mean([max(b(~Y)),min(b(Y>0))]),STAT,dof(1),dof(2),R(1),R(2),R(3),R(4),S,n,sign(Ic(nic1,nic2))>0); %p=spm_P_RF(1,0,u,dof,STAT,R,n);
                    
                    if sign(Ic(nic1,nic2))<0, Y=~Y; end % exclusion mask
                    if ~any(Y(:)>0), disp(['Warning! Output localizer volume #',num2str(nic1), ', contrast name ',symbols{signIc(nic1,nic2)},SPM{Ec(nic1,nic2)}.xCon(abs(Ic(nic1,nic2))).name,',  in contrast file ',fullfile(SPM{Ec(nic1,nic2)}.swd,filename),', contains no supra-threshold voxels']);
                    else disp([num2str(sum(Y(:)>0)),' voxels in output localizer volume #',num2str(nic1), ', contrast name ',symbols{signIc(nic1,nic2)},SPM{Ec(nic1,nic2)}.xCon(abs(Ic(nic1,nic2))).name,',  in contrast file ',fullfile(SPM{Ec(nic1,nic2)}.swd,filename)]); end
                    % computes conjunction across rows of Ic
                    switch(conjunction_type)
                        case 'and',     Z=Z&Y;
                        case 'or',      Z=Z|Y;
                        otherwise,      error('unrecognized conjunction type %s',conjuntion_type);
                    end
                else
                    error(['Incorrect threshold type option',thr_type{nic2}]);
                end
            end
            if numel(Z)>1
                if size(Ic,2)>1,
                    if ~any(Z(:)>0), disp(['Warning! Output localizer volume #',num2str(nic1), ',  in contrast file ',maskfilename{nic1},', contains no supra-threshold voxels']);
                    else disp([num2str(sum(Z(:)>0)),' voxels in output localizer volume #',num2str(nic1), ',  in contrast file ',maskfilename{nic1}]); end
                end
                cwd0=pwd;
                cd(fileparts(maskfilename{nic1}));
                fullfilename=filename; if numel(filename)>64, filename=['locT_conjunction_',char(mlreportgen.utils.hash(filename)),'.nii']; end
                Vo=struct(...
                    'fname',    filename,...
                    'dim',      a.dim,...
                    'dt',       [spm_type('uint8') spm_platform('bigend')],...
                    'mat',      a.mat,...
                    'pinfo',    [1;0;0],...
                    'descrip',  sprintf('SPM_SS LOCALIZER{%s}',cat(2,U{:})));
                Vo=spm_write_vol(Vo,Z);
                spm_jsonwrite([regexprep(filename,'\.nii$|\.img$',''),'.json'],struct('threshold_value',thr,'threshold_type',{thr_type},'conjunction_type',conjunction_type,'threshold_roifile',{parcelfile},'output_nvoxels',nnz(Z>0),'full_filename',fullfilename));
                cd(cwd0);
            end
        end
    end
else % CONJUNCTION MIN/MAX/PROD/SUM/OMNIBUS
    for nic1=1:size(Ic,1), % computes one separate thresholded volume per row of Ic (multiple columns are treated as a conjunction)
        filename='locT_conjunction_';
        for nic2=1:size(Ic,2),
            if Ec(nic1,nic2)~=Ec(nic1,1), temp=SPM{Ec(nic1,nic2)}.swd; temp(temp==filesep)='_';filename=[filename,temp]; end
            filename=[filename,symbols{signIc(nic1,nic2)},num2str(abs(Ic(nic1,nic2)),'%04d'),'_',thr_type{nic2},num2str(thr(nic2))];
            if nic2<size(Ic,2), filename=[filename,'_']; if ~strcmp(conjunction_type,'and'),filename=[filename,conjunction_type,'_']; end; end
        end
        if ~isempty(optionROIs)&&(any(strcmp(thr_type,'percentile-ROI-level'))||any(strcmp(thr_type,'Nvoxels-ROI-level'))) % note: avoid same localizer file for different optionROIs files
            filename=[filename,'_',char(mlreportgen.utils.hash(fileread(optionROIs(1).fname)))]; % note: requires Matlab 18b or above
            parcelfile={optionROIs.fname};
            %[nill,fname,nill]=fileparts(optionROIs(1).fname);
            %filename=[filename,'_',fname];
        else parcelfile={};
        end
        filename=[filename,'.img'];
        maskfilename{nic1}=fullfile(SPM{Ec(nic1,1)}.swd,filename);
        if overwrite||isempty(dir(maskfilename{nic1})),
            switch(conjunction_type)
                case 'min',     Z=nan;
                case 'max',     Z=nan;
                case 'prod',    Z=1;
                case 'sum',     Z=0;
                case 'omnibus', Z=0; dfZ=0;
                otherwise,      error('unrecognized conjunction type %s',conjuntion_type);
            end
            U={};
            for nic2=1:size(Ic,2),
                a=spm_vol(fullfile(SPM{Ec(nic1,nic2)}.swd,SPM{Ec(nic1,nic2)}.xCon(abs(Ic(nic1,nic2))).Vspm.fname));
                b=spm_read_vols(a);
                idx=find(~isnan(b)&b~=0);
                dof=[SPM{Ec(nic1,nic2)}.xCon(abs(Ic(nic1,nic2))).eidf,SPM{Ec(nic1,nic2)}.xX.erdf];
                STAT=SPM{Ec(nic1,nic2)}.xCon(abs(Ic(nic1,nic2))).STAT;
                R=SPM{Ec(nic1,nic2)}.xVol.R;
                S=SPM{Ec(nic1,nic2)}.xVol.S;
                n=1;
                Y=nan+zeros(size(b)); % p-value of individual contrast
                if isequal(conjunction_type,'omnibus'),
                    Y(idx)=b(idx);
                    if sign(Ic(nic1,nic2))<0, error('omnibus conjunction type does not allow -not flags in contrasts'); end % exclusion mask
                    switch(STAT),
                        case 'T',dfZ=dfZ+dof(2); Z=Z+Y*sqrt(dof(2));
                        case 'F',dfZ=dfZ+dof(2); Z=Z+sqrt(Y)*sqrt(dof(2)); % note: only exact for F(1,df) stats
                        otherwise, error('null');
                    end
                else
                    switch(STAT),
                        case 'Z',Y(idx)=1-spm_Ncdf(b(idx));
                        case 'T',Y(idx)=1-spm_Tcdf(b(idx),dof(2));
                        case 'X',Y(idx)=1-spm_Xcdf(b(idx),dof(2));
                        case 'F',Y(idx)=1-spm_Fcdf(b(idx),dof);
                        otherwise, error('null');
                    end
                    if sign(Ic(nic1,nic2))<0, Y=1-Y; end % exclusion mask
                    switch(conjunction_type)
                        case 'min',     Z=min(Z,Y);
                        case 'max',     Z=max(Z,Y);
                        case 'prod',    Z=Z.*Y;
                        case 'sum',     Z=Z+Y;
                        otherwise,      error('unrecognized conjunction type %s',conjuntion_type);
                    end
                end
                U{nic2}=sprintf('p=%f(%s),STAT=%s,dof=[%f,%f],R=[%f,%f,%f,%f],S=%f,n=%d,I=%d;',thr(nic2),thr_type{nic2},STAT,dof(1),dof(2),R(1),R(2),R(3),R(4),S,n,sign(Ic(nic1,nic2))>0); %p=spm_P_RF(1,0,u,dof,STAT,R,n);
            end                
            if isequal(conjunction_type,'omnibus'), Z=1-spm_Tcdf(Z/sqrt(dfZ), dfZ); end
            if ~isempty(strmatch(thr_type{nic2},{'percentile-whole-brain','percentile-ROI-level','Nvoxels-whole-brain','Nvoxels-ROI-level'},'exact')),
                b=nan+zeros(size(Z));
                idx=find(~isnan(Z));
                b(idx)=exp(-Z(idx)); 
                %switch(STAT),
                %    case 'Z',b(idx)=spm_invNcdf(max(0,min(1, 1-Z(idx) )));
                %    case 'T',b(idx)=spm_invTcdf(max(0,min(1, 1-Z(idx))),dof(2));
                %    case 'X',b(idx)=spm_invXcdf(max(0,min(1, 1-Z(idx))),dof(2));
                %    case 'F',b(idx)=spm_invFcdf(max(0,min(1, 1-Z(idx))),dof);
                %    otherwise, error('null');
                %end
                switch(thr_type{nic2}),
                    %case 'FWE',
                    %    u=spm_uc(thr(nic2),dof,STAT,R,n,S);
                    %    Y=double(b>u);
                    %case 'FDR',
                    %    Y=nan+zeros(size(b));
                    %    switch(STAT),
                    %        case 'Z',Y(idx)=1-spm_Ncdf(b(idx));
                    %        case 'T',Y(idx)=1-spm_Tcdf(b(idx),dof(2));
                    %        case 'X',Y(idx)=1-spm_Xcdf(b(idx),dof(2));
                    %        case 'F',Y(idx)=1-spm_Fcdf(b(idx),dof);
                    %        otherwise, error('null');
                    %    end
                    %    Y(:)=spm_ss_fdr(Y(:));
                    %    Y=double(Y<thr(nic2));
                    %case 'none'
                    %    u=spm_u(thr(nic2),dof,STAT);
                    %    Y=double(b>u);
                    case 'percentile-whole-brain'
                        Y=zeros(size(b));
                        tN=b(idx);
                        tN(tN==0|isnan(tN))=-inf;
                        tNthreshold=sort(tN);
                        tNthreshold=tNthreshold(max(1,min(numel(tN), ceil(numel(tN)*(1-thr(nic2))))));
                        Y(idx)=double(tN>tNthreshold);
                    case 'percentile-ROI-level'
                        if isempty(optionROIs), error('no ROI file provided for percentile-ROI-level option'); end
                        clear XYZ;
                        [XYZ{1},XYZ{2},XYZ{3}]=ndgrid(1:a(1).dim(1),1:a(1).dim(2),1:a(1).dim(3));
                        XYZ=reshape(cat(4,XYZ{:}),[],3)';XYZ=a(1).mat*cat(1,XYZ,ones(1,size(XYZ,2)));
                        frois=round(spm_get_data(optionROIs,pinv(optionROIs.mat)*XYZ));
                        nrois=max(frois(:));
                        Y=zeros(size(b));
                        for nroi=1:nrois,
                            idxroi=find(frois==nroi);
                            if ~isempty(idxroi)
                                tN=b(idxroi);
                                tN(tN==0|isnan(tN))=-inf;
                                tNthreshold=sort(tN);
                                tNthreshold=tNthreshold(max(1,min(numel(tN), ceil(numel(tN)*(1-thr(nic2))))));
                                Y(idxroi)=double(tN>tNthreshold);
                            end
                        end
                    case 'Nvoxels-whole-brain'
                        Y=zeros(size(b));
                        tN=b(idx);
                        tN(tN==0|isnan(tN))=-inf;
                        [tNthreshold,tNidx]=sort(tN);
                        tNidx(tNidx)=1:numel(tNidx);
                        Y(idx)=double(tNidx>numel(tN)-thr(nic2)&~isinf(tN));
                    case 'Nvoxels-ROI-level'
                        if isempty(optionROIs), error('no ROI file provided for Nvoxels-ROI-level option'); end
                        clear XYZ;
                        [XYZ{1},XYZ{2},XYZ{3}]=ndgrid(1:a(1).dim(1),1:a(1).dim(2),1:a(1).dim(3));
                        XYZ=reshape(cat(4,XYZ{:}),[],3)';XYZ=a(1).mat*cat(1,XYZ,ones(1,size(XYZ,2)));
                        frois=round(spm_get_data(optionROIs,pinv(optionROIs.mat)*XYZ));
                        nrois=max(frois(:));
                        Y=zeros(size(b));
                        for nroi=1:nrois,
                            idxroi=find(frois==nroi);
                            if ~isempty(idxroi)
                                tN=b(idxroi);
                                tN(tN==0|isnan(tN))=-inf;
                                [tNthreshold,tNidx]=sort(tN);
                                tNidx(tNidx)=1:numel(tNidx);
                                Y(idxroi)=double(tNidx>numel(tN)-thr(nic2)&~isinf(tN));
                            end
                        end
                end
                Z=Y;
            else
                error(['Incorrect threshold type option',thr_type{nic2}]);
            end
            if numel(Z)>1
                if size(Ic,2)>1,
                    if ~any(Z(:)>0), disp(['Warning! Output localizer volume #',num2str(nic1), ',  in contrast file ',maskfilename{nic1},', contains no supra-threshold voxels']);
                    else disp([num2str(sum(Z(:)>0)),' voxels in output localizer volume #',num2str(nic1), ',  in contrast file ',maskfilename{nic1}]); end
                end
                cwd0=pwd;
                cd(fileparts(maskfilename{nic1}));
                fullfilename=filename; if numel(filename)>64, filename=['locT_conjunction_',char(mlreportgen.utils.hash(filename)),'.nii']; end
                Vo=struct(...
                    'fname',    filename,...
                    'dim',      a.dim,...
                    'dt',       [spm_type('uint8') spm_platform('bigend')],...
                    'mat',      a.mat,...
                    'pinfo',    [1;0;0],...
                    'descrip',  sprintf('SPM_SS LOCALIZER{%s}',cat(2,U{:})));
                Vo=spm_write_vol(Vo,Z);
                spm_jsonwrite([regexprep(filename,'\.nii$|\.img$',''),'.json'],struct('threshold_value',thr,'threshold_type',{thr_type}, 'conjunction_type',conjunction_type, 'threshold_roifile',{parcelfile},'output_nvoxels',nnz(Z>0),'full_filename',fullfilename));
                cd(cwd0);
            end
        end
    end
end
end




