function ss=spm_ss_estimate_ROI(ss)
% SPM_SS_ESTIMATE_ROI subject-specific ROI-based model estimation
% 
% ss=spm_ss_estimate_ROI(ss)
% see SPM_SS_DESIGN, SPM_SS_ESTIMATE

if nargin<1, 
    str='Select spm_ss*.mat analysis file';
    disp(str);
    Pdefault='';objname=findobj('tag','spm_ss');if numel(objname)==1,objdata=get(objname,'userdata');if isfield(objdata,'files_spm_ss'),Pdefault=objdata.files_spm_ss;end;end;
    P=spm_select(1,'^SPM_ss.*\.mat$',str,{Pdefault});
    if numel(objname)==1&&~isempty(P),objdata.files_spm_ss=P;set(objname,'userdata',objdata);end;
    load(P);
    ss.swd=fileparts(P);
end

cwd=pwd;

% explicit mask
if ~isempty(ss.ExplicitMasking),XM=spm_vol(ss.ExplicitMasking);else XM=[];end

% creates transformed measures
neffects=size(ss.EffectOfInterest{1},1);
k=0;
for ne=1:neffects,
    for n=1:ss.n,
        for m=1:numel(ss.Localizer{n}),
            [pth2,nm2,ext2,num2]=spm_fileparts(ss.EffectOfInterest{n}{ne,m});
            Yvolume=[nm2,ext2,num2];
            k=k+1;
            ss.PY{k}=fullfile(pth2,Yvolume);
        end
    end
end
k=0;
ss.PV=[];
ss.PV_subjnames={};
for n=1:ss.n,
    for m=1:numel(ss.Localizer{n}),
        [pth1,nm1,ext1,num1]=spm_fileparts(ss.Localizer{n}{m});
        Nvolume=[nm1,ext1,num1];
        k=k+1;
        ss.PN{k}=fullfile(pth1,Nvolume);
        ss.PV(n,k)=1/numel(ss.Localizer{n});
    end
    ss.PV_subjnames{n}=num2str(n); 
    if numel(ss.Localizer{n})>0
        [pth1a,pth1b]=fileparts(pth1);
        if ~isempty(regexp(pth1b,'^model_|^firstlevel_'))&&~isempty(pth1a), [pth1b,pth1a]=fileparts(pth1a); ss.PV_subjnames{n}=pth1a; end
    end
end

cd(ss.swd); 
ss.VN=spm_vol(char(ss.PN));                             % localizer masks: [sum(cross-validation runs across subjects) , 1]
ss.VY=reshape(spm_vol(char(ss.PY)),numel(ss.VN),[]);    % effects: [sum(cross-validation runs across subjects) , neffects]
ss.PY=reshape(ss.PY,numel(ss.VN),[]);

ss.refspace.masks={};
ss.refspace.effects={};
ss.refspace.effectsinsamespace=false(1,ss.n);            % all effect volumes are in same space
ss.refspace.effectsinsamespaceasmask=false(1,ss.n);      % all effect volumes are in same space as mask volume
for n=1:ss.n,
    idxk=ss.PV(n,:)>0;
    idxk1=find(idxk,1);
    if ~isempty(idxk1)
        [pth1,nm1,ext1,num1]=spm_fileparts(ss.VY(idxk1,1).fname);
        [pth1a,pth1b]=fileparts(pth1);
        if ~isempty(regexp(pth1b,'^model_|^firstlevel_'))&&~isempty(pth1a), [pth1b,pth1a]=fileparts(pth1a); ss.PV_subjnames{n}=pth1a; end
        nvy=ss.VN(idxk1); 
        ss.refspace.masks{n}=struct('mat',nvy.mat,'dim',nvy.dim);
        nvy=ss.VY(idxk,:); 
        ss.refspace.effects{n}=struct('mat',{nvy.mat},'dim',{nvy.dim});
        try, ss.refspace.effectsinsamespace(n)=all(reshape(~diff(cat(4,nvy.mat),1,4),[],1))&all(reshape(~diff(cat(4,nvy.dim),1,4),[],1)); end
        try, ss.refspace.effectsinsamespaceasmask(n)=all(reshape(~diff(cat(4,ss.VN(idxk1).mat,nvy.mat),1,4),[],1))&all(reshape(~diff(cat(4,ss.VN(idxk1).dim,nvy.dim),1,4),[],1)); end
    end
end

% Creates overlap maps
ssPM1=['Overlap',ext1];
p=0;
for n=1:ss.n,
    idx=find(ss.PV(n,:));
    p1=0;
    for k=1:numel(idx),
        a1=ss.VN(idx(k));%spm_vol(ss.PN{idx(k)});
        if n==1, b1=spm_read_vols(a1);  % note: if subjects data is not in same space overlap maps are computed in the (arbitrary) space of first subject
        else b1=reshape(spm_get_data(a1,pinv(a1.mat)*tXYZ),e0(1).dim);
        end
        p1=p1+(b1>0)/numel(idx);
    end
    p=p+(p1>.5);
    if n==1,
        ss.PL=['Localizer',ext1];
        e0=struct('fname',ss.PL,'descrip','spm_ss (localizer mask for each subject)','mat',a1.mat,'dim',a1.dim,'n',[1,1],'pinfo',[1;0;0],'dt',[spm_type('float32'),spm_platform('bigend')]);
        try, spm_unlink(e0.fname); end
        e0=repmat(e0,[ss.n,1]);for nb=1:ss.n,e0(nb).n=[nb,1];end
        e0=spm_create_vol(e0);
        [XYZ{1},XYZ{2},XYZ{3}]=ndgrid(1:a1.dim(1),1:a1.dim(2),1:a1.dim(3));
        tXYZ=reshape(cat(4,XYZ{:}),[],3)';tXYZ=a1.mat*cat(1,tXYZ,ones(1,size(tXYZ,2))); % world-coordinates for first subject
    end
    spm_write_vol(e0(n),p1);
    if 1 % check/debug this
        fname='Localizer.csv';
        fh=fopen(fullfile(ss.swd,fname),'wt');
        fprintf(fh,'subject,index,weight,i,j,k,x(mm),y(mm),z(mm)\n');
        for k=reshape(find(p1>0),1,[]), fprintf(fh,'%d,%d,%.4f,%d,%d,%d,%d,%d,%d\n',n,k,p1(k),XYZ{1}(k),XYZ{2}(k),XYZ{3}(k),round(tXYZ(1,k)),round(tXYZ(2,k)),round(tXYZ(3,k))); end
        fclose(fh);
    end
end
e1=struct('fname',ssPM1,'descrip','spm_ss (inter-subject overlap map)','mat',e0(1).mat,'dim',e0(1).dim,'dt',[spm_type('float32'),spm_platform('bigend')]);
e1=spm_write_vol(e1,p/ss.n);
if ss.typen==2,
    ssPM2=['sOverlap',ext1];
    spm_smooth(e1,ssPM2,ss.smooth*[1,1,1]);
    ss.PM2=ssPM2;
end
fprintf(1,'\n');
ss.PM1=ssPM1;

% defines fROIs
if ss.typen==2,
    ss.PM=['fROIs',ext1];
    disp('GcSS defining ROIs. Please wait...');
    a2=spm_vol(ss.PM2);
    b2=spm_read_vols(a2); 
    if ~isempty(XM),
        [XYZ{1},XYZ{2},XYZ{3}]=ndgrid(1:ss.refspace.masks{1}.dim(1),1:ss.refspace.masks{1}.dim(2),1:ss.refspace.masks{1}.dim(3));
        tXYZ=reshape(cat(4,XYZ{:}),[],3)';tXYZ=cat(1,tXYZ,ones(1,size(tXYZ,2))); % reference: world-coordinates for first localizer mask volume for first subject
        tXM=reshape(double(spm_get_data(XM,pinv(XM.mat)*ss.refspace.masks{1}.mat*tXYZ)>0),ss.refspace.masks{1}.dim(1:3));
    else tXM=1;
    end
    b3=spm_ss_watershed(-b2,find(b2>=ss.overlap_thr_vox&tXM));
    fprintf('Done. Defined %d regions\n',max(b3(:)));
    a3=struct('fname',ss.PM,'mat',a2.mat,'dim',a2.dim,'dt',[spm_type('int16') spm_platform('bigend')],'pinfo',[1;0;0]);
    spm_write_vol(a3,b3);  
else
    ss.PM=ss.ManualROIs;
end
ss.VM=spm_vol(char(ss.PM));

for n=1:ss.n 
    % frois to subject-space
    idxk=ss.PV(n,:)>0;
    idxk1=find(idxk,1);
    [XYZ{1},XYZ{2},XYZ{3}]=ndgrid(1:ss.VN(idxk1).dim(1),1:ss.VN(idxk1).dim(2),1:ss.VN(idxk1).dim(3));
    tXYZ=reshape(cat(4,XYZ{:}),[],3)';tXYZ=cat(1,tXYZ,ones(1,size(tXYZ,2))); % reference: world-coordinates for first localizer mask volume for subject n
    vm=ss.VM(min(numel(ss.VM),n));
    frois{n}=reshape(round(spm_get_data(vm,pinv(vm.mat)*ss.VN(idxk1).mat*tXYZ)),ss.VN(idxk1).dim(1:3));
    if n==1, nrois=max(frois{n}(:));
    else nrois=max(nrois,max(frois{n}(:)));
    %elseif max(frois{n}(:))~=nrois, error('All subject-specific ROI files should contain the same number of ROIs'); 
    end
    if ~isempty(XM),
        tXM=reshape(double(spm_get_data(XM,pinv(XM.mat)*ss.VN(idxk1).mat*tXYZ)>0),ss.VN(idxk1).dim(1:3));
        frois{n}=tXM.*frois{n};
    end
    
    % creates QA images in subject-space
    a4=struct('fname',[sprintf('QA_parcels.%s',ss.PV_subjnames{n}),ext1],'mat',ss.refspace.masks{n}.mat,'dim',ss.refspace.masks{n}.dim,'dt',[spm_type('float32') spm_platform('bigend')],'descrip',sprintf('QA display parcels for subject#%d',n),'pinfo',[1;0;0]);
    spm_write_vol(a4,reshape(frois{n},ss.VN(idxk1).dim(1:3)));
    a4=struct('fname',[sprintf('QA_effects.%s',ss.PV_subjnames{n}),ext1],'mat',ss.refspace.masks{n}.mat,'dim',ss.refspace.masks{n}.dim,'dt',[spm_type('float32') spm_platform('bigend')],'descrip',sprintf('QA display average all effects for subject#%d',n),'pinfo',[1;0;0]);
    tvy=reshape(ss.VY(idxk,:),[],1);
    if isempty(idxk1),                              tY=zeros(numel(tvy),size(tXYZ,2));
    else
        if ss.refspace.effectsinsamespaceasmask(n), tY=spm_get_data(tvy,tXYZ);
        elseif ss.refspace.effectsinsamespace(n),   tY=spm_get_data(tvy,pinv(tvy(1).mat)*ss.refspace.masks{n}.mat*tXYZ);
        else                                        tY=zeros(numel(tvy),size(tXYZ,2)); for nk=1:numel(tvy), tY(nk,:)=spm_get_data(tvy(nk),pinv(tvy(nk).mat)*ss.refspace.masks{n}.mat*tXYZ); end
        end
    end
    tY(isnan(tY))=0;
    spm_write_vol(a4,reshape(mean(tY,1),ss.VN(idxk1).dim(1:3)));
    a4=struct('fname',[sprintf('QA_localizer.%s',ss.PV_subjnames{n}),ext1],'mat',ss.refspace.masks{n}.mat,'dim',ss.refspace.masks{n}.dim,'dt',[spm_type('float32') spm_platform('bigend')],'descrip',sprintf('QA display average localizer mask for subject#%d',n),'pinfo',[1;0;0]);
    tvy=reshape(ss.VN(idxk),[],1);
    if isempty(idxk1),                              tN=zeros(numel(tvy),size(tXYZ,2));
    else                                            tN=spm_get_data(tvy,tXYZ);
    end
    tN(isnan(tN))=0;
    tN=reshape(mean(tN,1),ss.VN(idxk1).dim(1:3));
    spm_write_vol(a4,tN);
end
[ss.VM_roinames,ss.VM_roiids]=spm_ss_roilabels(ss.VM(1).fname);
if isempty(ss.VM_roinames), ss.VM_roinames=arrayfun(@num2str,1:nrois,'uni',0); ss.VM_roiids=1:nrois; 
else nrois=numel(ss.VM_roiids); 
end

% analysis

Nb=[size(ss.X,2),neffects];
extname=['_',ss.type];
VB=struct('fname',['spm_ss',extname,'_beta.img'],...
    'mat',ss.VN(1).mat,...
    'dim',ss.VN(1).dim,...
    'n',[1,1],...
    'pinfo',[1;0;0],...
    'dt',[spm_type('float32'),spm_platform('bigend')],...
    'descrip','spm_ss (effect sizes parameter estimates)');
try, spm_unlink(VB.fname); end
VB=repmat(VB,[prod(Nb),1]);for nb=1:prod(Nb),VB(nb).n=[nb,1];end
VB=spm_create_vol(VB);
% VE=struct('fname',['spm_ss',extname,'_rss.img'],...
%     'mat',ss.VN(1).mat,...
%     'dim',ss.VN(1).dim,...
%     'n',[1,1],...
%     'pinfo',[1;0;0],...
%     'dt',[spm_type('float32'),spm_platform('bigend')],...
%     'descrip','spm_ss (residual sum squares)');
% VE=repmat(VE,[Nb(2)*Nb(2),1]);for nb=1:Nb(2)*Nb(2),VE(nb).n=[nb,1];end
% VE=spm_create_vol(VE);
VO=struct('fname',['spm_ss',extname,'_overlap.img'],...
    'mat',ss.VN(1).mat,...
    'dim',ss.VN(1).dim,...
    'pinfo',[1;0;0],...
    'dt',[spm_type('float32'),spm_platform('bigend')],...
    'descrip','spm_ss (proportion overlap)');
VO=spm_create_vol(VO);


Bplane=nan+zeros([Nb,nrois]);Cplane=zeros(ss.n,nrois);Eplane=nan+zeros(Nb(2),Nb(2),nrois);Oplane=nan+zeros(1,nrois);Zplane=zeros(ss.n,Nb(2),nrois);Nplane=nan+zeros([ss.n,nrois]);Pplane=nan+zeros([1,nrois]);Mplane=nan+zeros([numel(ss.VN),nrois]);

fprintf('Performing model estimation...');
for nroi=1:nrois,
    fprintf(1,'.');
    Y=zeros(numel(ss.VY),1);
    N=zeros(numel(ss.VN),1);
    for n=1:ss.n
        idx=find(frois{n}==ss.VM_roiids(nroi));
        [idx1,idx2,idx3]=ind2sub(size(frois{n}),idx);
        xyz=[idx1,idx2,idx3,ones(numel(idx1),1)]';
        Nplane(n,nroi)=numel(idx);

        idxk=ss.PV(n,:)>0;
        idxk1=find(idxk,1);
        tN=spm_get_data(ss.VN(idxk),xyz);
        tvy=reshape(ss.VY(idxk,:),[],1);
        if isempty(idxk1),                              tY=zeros(numel(tvy),size(xyz,2)); 
        else 
            if ss.refspace.effectsinsamespaceasmask(n), tY=spm_get_data(tvy,xyz);
            elseif ss.refspace.effectsinsamespace(n),   tY=spm_get_data(tvy,pinv(tvy(1).mat)*ss.refspace.masks{n}.mat*xyz);
            else                                        tY=zeros(numel(tvy),size(xyz,2)); for nk=1:numel(tvy), tY(nk,:)=spm_get_data(tvy(nk),pinv(tvy(nk).mat)*ss.refspace.masks{n}.mat*xyz); end
            end
        end
        Z=(tN==0)|isnan(tN);tY(repmat(Z,[Nb(2),1]))=0;tN(Z)=0;tY(isnan(tY))=0;
        Mplane(idxk,nroi)=sum(tN,2);
        tY=mean(tY.*repmat(tN,[Nb(2),1]),2);
        tN=mean(tN,2);
        Y(repmat(idxk',[Nb(2),1]),:)=tY;
        N(idxk)=tN;
    end
%         %data=cell(size(Y,1),1); for n1=1:size(Y,1), data{n1}=Y(n1,N(1+rem(n1-1,size(N,1)),:)); end; save(['data_roi',num2str(nroi),'.mat'],'data');     
%     Y=spm_get_data(ss.VY,xyz);
%     N=spm_get_data(ss.VN,xyz);
%     Z=(N==0)|isnan(N);Y(repmat(Z,[Nb(2),1]))=0;N(Z)=0;Y(isnan(Y))=0;
%     Mplane(:,nroi)=sum(N,2);
%     Y=mean(Y.*repmat(N,[Nb(2),1]),2);
%     N=mean(N,2);
    Y=reshape(Y./max(eps,repmat(N,[Nb(2),1])),[size(N,1),Nb(2)]);
    Y=ss.PV*Y;
    N=1./(ss.PV*(1./max(eps,N)));
    sN=mean(N>1e-4,1);
    if sN>0,
        if strcmpi(ss.estimation,'ols')
            iC=double(N>1e-4);%handles missing-data
        else
            n=N;
            y=Y.*sqrt(n(:,ones(1,Nb(2))));
            x=ss.X.*sqrt(n(:,ones(1,Nb(1))));
            e=Y-ss.X*(pinv(x'*x)*(x'*y));
            [nill,iC]=spm_ss_fls({e,n});%covariance estimation
        end
        y=Y.*iC(:,ones(1,Nb(2)));%whitening
        x=ss.X.*iC(:,ones(1,Nb(1)));
        [b,ee]=spm_ss_glm('estimate',x,y);
        Bplane(:,:,nroi)=b;
        Cplane(:,nroi)=iC;
        Eplane(:,:,nroi)=ee;
    end 
    Pplane(nroi)=mean(N);
    Oplane(nroi)=sN;
    Zplane(:,:,nroi)=Y;
end
fprintf(1,'\n');

% save files
nb=1;for nb1=1:Nb(1),for nb2=1:Nb(2),z=nan+zeros(size(frois{1}));for nroi=1:nrois,z(frois{1}==ss.VM_roiids(nroi))=Bplane(nb1,nb2,nroi);end; spm_write_vol(VB(nb),z);nb=nb+1;end;end
% nb=1;for nb1=1:Nb(2),for nb2=1:Nb(2),z=nan+zeros(size(frois));for nroi=1:nrois,z(frois==nroi)=Eplane(nb1,nb2,nroi);end; spm_write_vol(VE(nb),z);nb=nb+1;end;end
z=nan+zeros(size(frois{1}));for nroi=1:nrois,z(frois{1}==ss.VM_roiids(nroi))=Oplane(nroi);end; spm_write_vol(VO,z);
% disp(['created beta volume       : ',fullfile(ss.swd,VB(1).fname),' - ',num2str(Nb),' volume(s)']); 
% disp(['created rss volume        : ',fullfile(ss.swd,VE(1).fname)]); 
% disp(['created overlap volume    : ',fullfile(ss.swd,VO(1).fname)]); 

ss.estimate=struct('BETA',VB,'OVERLAP',VO,'beta',Bplane,'rss',Eplane,'whitening',Cplane,'overlap',Oplane,'voxels',Nplane,'coverage',Pplane,'qa',Mplane,'y',Zplane); 
% ss.estimate=struct('BETA',VB,'RSS',VE,'OVERLAP',VO,'beta',Bplane,'rss',Eplane,'whitening',Cplane,'overlap',Oplane,'voxels',Nplane,'coverage',Pplane,'qa',Mplane,'y',Zplane); 
save(fullfile(ss.swd,['SPM_ss',extname,'.mat']),'ss');
disp(['Analysis file saved: ',fullfile(ss.swd,['SPM_ss',extname,'.mat'])]);

% csv output
nidxs=zeros(ss.n,1);
nidxs2=zeros(2,size(ss.estimate.qa,1));
for nf=1:size(ss.estimate.qa,1),
    [nill,idxs]=max(ss.PV(:,nf));
    nidxs(idxs)=nidxs(idxs)+1;
    nidxs2(1,nf)=idxs;
    nidxs2(2,nf)=nidxs(idxs);
end

if 0, % note: this will be removed in the future; change to 1 to create old-format _data.csv file
    % Original csv output
    fname=['spm_ss',extname,'_data.csv'];
    fh=fopen(fullfile(ss.swd,fname),'wt');
    fprintf(fh,'Data\n');
    fprintf(fh,'ROI#,average ROI size,average localizer mask size,inter-subject overlap,');
    for ns=1:ss.n,for ne=1:Nb(2),fprintf(fh,'Subject#%d[%d]',ns,ne); if ne<Nb(2)||ns<ss.n, fprintf(fh,','); else fprintf(fh,'\n'); end; end; end
    for nroi=1:nrois,
        fprintf(fh,'%d,%d,%d,%f,',nroi,round(mean(ss.estimate.voxels(:,nroi))),round(mean(ss.estimate.voxels(:,nroi))*ss.estimate.coverage(nroi)),ss.estimate.overlap(nroi));
        for ns=1:ss.n, for ne=1:Nb(2),fprintf(fh,'%f',Zplane(ns,ne,nroi)); if ne<Nb(2)||ns<ss.n, fprintf(fh,','); else fprintf(fh,'\n'); end; end; end
    end
    fprintf(fh,'\nWeights\n');
    fprintf(fh,'ROI#,');
    for ns=1:ss.n,fprintf(fh,'Subject#%d',ns); if ns<ss.n, fprintf(fh,','); else fprintf(fh,'\n'); end; end;
    for nroi=1:nrois,
        fprintf(fh,'%d,',nroi);
        for ns=1:ss.n, fprintf(fh,'%f',Cplane(ns,nroi)); if ns<ss.n, fprintf(fh,','); else fprintf(fh,'\n'); end; end;
    end
    fprintf(fh,'\nquality control (localizer mask sizes)\n');
    fprintf(fh,'Subject#,Session/Partition#,filename'); for nroi=1:nrois,fprintf(fh,',ROI#%d',nroi);end;fprintf(fh,'\n');
    for nf=1:size(ss.estimate.qa,1),
        fprintf(fh,'%d,%d,%s',nidxs2(1,nf),nidxs2(2,nf),ss.PN{nf});
        for nroi=1:nrois,fprintf(fh,',%d',round(ss.estimate.qa(nf,nroi)));end
        fprintf(fh,'\n');
    end
    fclose(fh);
end
% New format output
if size(ss.EffectOfInterest_contrasts,2)==Nb(2), effect_names=ss.EffectOfInterest_contrasts; 
else effect_names=arrayfun(@num2str,1:Nb(2),'uni',0); 
end
fname=['spm_ss',extname,'_data.details.Weight.csv'];
fh=fopen(fullfile(ss.swd,fname),'wt');
fprintf(fh,'ROI,Subject,Weight\n');
for nroi=1:nrois,for ns=1:ss.n, fprintf(fh,'%s,%s,%f\n',ss.VM_roinames{nroi},ss.PV_subjnames{ns},Cplane(ns,nroi)); end; end
fclose(fh);
fname=['spm_ss',extname,'_data.summaries.EffectSize.csv'];
fh=fopen(fullfile(ss.swd,fname),'wt');
fprintf(fh,'ROI,Effect,MeanEffect,StdEffect,StderrEffect\n');
for nroi=1:nrois,for ne=1:Nb(2), fprintf(fh,'%s,%s,%f,%f,%f\n',ss.VM_roinames{nroi},effect_names{1,ne},mean(Zplane(:,ne,nroi)),std(Zplane(:,ne,nroi)),std(Zplane(:,ne,nroi))/sqrt(size(Zplane,1))); end; end
fclose(fh);
fname=['spm_ss',extname,'_data.details.ROIsize.csv'];
fh=fopen(fullfile(ss.swd,fname),'wt');
fprintf(fh,'ROI,Subject,Session/Partition,LocalizerSize\n');
for nroi=1:nrois, for nf=1:size(ss.estimate.qa,1), fprintf(fh,'%s,%s,%d,%d\n',ss.VM_roinames{nroi},ss.PV_subjnames{nidxs2(1,nf)},nidxs2(2,nf),round(ss.estimate.qa(nf,nroi)));end; end
fclose(fh);
fname=['spm_ss',extname,'_data.summaries.ROIsize.csv'];
fh=fopen(fullfile(ss.swd,fname),'wt');
fprintf(fh,'ROI,ROISize,LocalizerSize,LocalizerIntersubjectOverlap\n');
for nroi=1:nrois, fprintf(fh,'%d,%d,%d,%f\n',ss.VM_roinames{nroi},round(mean(ss.estimate.voxels(:,nroi))),round(mean(ss.estimate.voxels(:,nroi))*ss.estimate.coverage(nroi)),ss.estimate.overlap(nroi)); end
fname=['spm_ss',extname,'_data.details.SourceFiles.csv'];
fh=fopen(fullfile(ss.swd,fname),'wt');
fprintf(fh,'Subject,Session/Partition,Localizer');
for ne=1:Nb(2), fprintf(fh,',%s',effect_names{1,ne}); end; fprintf(fh,'\n'); 
for nf=1:size(ss.estimate.qa,1),
    fprintf(fh,'%s,%d,%s',ss.PV_subjnames{nidxs2(1,nf)},nidxs2(2,nf),ss.PN{nf}); 
    for ne=1:size(ss.PY,2), fprintf(fh,',%s',ss.PY{nf,ne}); end; fprintf(fh,'\n');
end
fclose(fh);
fname=['spm_ss',extname,'_data.csv'];
fh=fopen(fullfile(ss.swd,fname),'wt');
fprintf(fh,'ROI,Subject,Effect, LocalizerSize,EffectSize\n');
for nroi=1:nrois,for ns=1:ss.n, for ne=1:Nb(2), fprintf(fh,'%s,%s,%s, %d,%f\n',ss.VM_roinames{nroi},ss.PV_subjnames{ns},effect_names{1,ne}, round(ss.PV(ns,:)*ss.estimate.qa(:,nroi)),Zplane(ns,ne,nroi)); end; end; end
fclose(fh);

% creates QA plots
try, spm_ss qacreate;
catch, fprintf('warning: unable to generate QA plots (possibly missing graphic display capabilities). Please use "spm_ss qa ''%s''" syntax to create these plots at a later time\n',ss.swd);
end

% estimates defined contrasts
ss=spm_ss_contrast_ROI(ss);
cd(cwd);

end


function conn_savetextfile(tfilename,data,names,descrip)
% conn_savetextfile saves numeric data data to text file
% conn_savetextfile(tfilename,data [,names,descrip])
%

if nargin<4||isempty(descrip), descrip={}; end
if nargin<3||isempty(names), names={}; end
[nill,nill,tfileext]=fileparts(tfilename);
switch(tfileext)
    case '.mat'
        if ~isempty(names)&&~isempty(descrip), save(tfilename,'data','names','descrip');
        elseif ~isempty(names), save(tfilename,'data','names');
        else save(tfilename,'data');
        end
    otherwise,
        if strcmp(tfileext,'.txt'), names=regexprep(names,'\s','');
        else                        names=regexprep(names,'\,','');
        end
        fh=fopen(tfilename,'wt');
        for n1=1:numel(names),
            if isempty(names{n1}), names{n1}='-'; end
            fprintf(fh,'%s',names{n1});
            if n1<numel(names)&&strcmp(tfileext,'.csv'), fprintf(fh,','); 
            elseif n1<numel(names), fprintf(fh,' '); 
            else fprintf(fh,'\n'); 
            end
        end
        for n2=1:size(data,1),
            for n1=1:size(data,2),
                if iscell(data(n2,n1))&&ischar(data{n2,n1}), fprintf(fh,'%s',data{n2,n1});
                else fprintf(fh,'%s',mat2str(data(n2,n1)));
                end
                if n1<size(data,2)&&strcmp(tfileext,'.csv'), fprintf(fh,','); 
                elseif n1<size(data,2), fprintf(fh,' '); 
                else fprintf(fh,'\n'); 
                end
            end
        end
        fclose(fh);
end
end


