function bcc=spm_bcc_estimate(bcc)
% SPM_BCC_ESTIMATE between-condition-correlations ROI-based model estimation
% 
% bcc=spm_bcc_estimate(bcc)
% see SPM_BCC_DESIGN

% spm_ss initialization
spm_ss init silent;

if nargin<1, 
    str='Select spm_bcc*.mat analysis file';
    disp(str);
    Pdefault='';objname=findobj('tag','spm_bcc');if numel(objname)==1,objdata=get(objname,'userdata');if isfield(objdata,'files_spm_bcc'),Pdefault=objdata.files_spm_bcc;end;end;
    P=spm_select(1,'^SPM_bcc.*\.mat$',str,{Pdefault});
    if numel(objname)==1&&~isempty(P),objdata.files_spm_bcc=P;set(objname,'userdata',objdata);end;
    load(P);
    bcc.swd=fileparts(P);
end

bcc=spm_bcc_design(bcc); % checks that information is complete

cwd=pwd;

% explicit mask
%if ~isempty(bcc.ExplicitMasking),XM=spm_read_vols(spm_vol(bcc.ExplicitMasking))>0;else XM=1;end
if ~isempty(bcc.ExplicitMasking),XM=spm_vol(bcc.ExplicitMasking);else XM=[];end

neffects1=size(bcc.EffectOfInterest1{1},1);
neffects2=size(bcc.EffectOfInterest2{1},1);
neffects=neffects1*neffects2;
try
    effect_names=cellfun(@(a,b)[a,'/',b],repmat(bcc.EffectOfInterest1_contrasts(:),1,numel(bcc.EffectOfInterest2_contrasts)),repmat(bcc.EffectOfInterest2_contrasts(:)',numel(bcc.EffectOfInterest1_contrasts),1),'uni',0);
catch
    effect_names={};
end
k=0;
bcc.PV=[];
for n=1:bcc.n,
    ncross=size(bcc.EffectOfInterest1{n},2);
    for m=1:ncross
        k=k+1;
        for ne=1:neffects1,
            [pth2,nm2,ext2,num2]=spm_fileparts(bcc.EffectOfInterest1{n}{ne,m});
            Yvolume=[nm2,ext2,num2];
            bcc.PY1{k,ne}=fullfile(pth2,Yvolume);
        end
        for ne=1:neffects2,
            [pth2,nm2,ext2,num2]=spm_fileparts(bcc.EffectOfInterest2{n}{ne,m});
            Yvolume=[nm2,ext2,num2];
            bcc.PY2{k,ne}=fullfile(pth2,Yvolume);
        end
        bcc.PV(n,k)=1/ncross;
    end
end

cd(bcc.swd);
bcc.VY1=reshape(spm_vol(char(bcc.PY1)),size(bcc.PY1));
bcc.VY2=reshape(spm_vol(char(bcc.PY2)),size(bcc.PY2));

% defines fROIs
bcc.PM=bcc.ManualROIs;
if isempty(bcc.PM), bcc.VM=[];
else bcc.VM=spm_vol(char(bcc.PM));
end
% if numel(bcc.VM)>1,error('multiple ROI image files not supported'); end
% [XYZ{1},XYZ{2},XYZ{3}]=ndgrid(1:bcc.VY1(1).dim(1),1:bcc.VY1(1).dim(2),1:bcc.VY1(1).dim(3));
% XYZ=reshape(cat(4,XYZ{:}),[],3)';XYZ=bcc.VY1(1).mat*cat(1,XYZ,ones(1,size(XYZ,2)));
% if isempty(bcc.VM), frois=ones(bcc.VY1(1).dim); 
% else                frois=reshape(round(spm_get_data(bcc.VM,pinv(bcc.VM.mat)*XYZ))',[bcc.VY1(1).dim(1:3),numel(bcc.VM)]);
% end
% if numel(XM)>1&&any(size(frois)~=size(XM)),error('mismatched volume dimensions between functional volumes and explicit mask images'); end
% frois=XM.*frois;
% nrois=max(frois(:));
for n=1:bcc.n % frois in subject-space
    idxk=bcc.PV(n,:)>0;
    idxk1=find(idxk,1);
    [XYZ{1},XYZ{2},XYZ{3}]=ndgrid(1:bcc.VY1(idxk1).dim(1),1:bcc.VY1(idxk1).dim(2),1:bcc.VY1(idxk1).dim(3));
    tXYZ=reshape(cat(4,XYZ{:}),[],3)';tXYZ=bcc.VY1(idxk1).mat*cat(1,tXYZ,ones(1,size(tXYZ,2))); % voxel-coordinates for first functional volume for subject n
    if isempty(bcc.VM), frois{n}=ones(bcc.VY1(idxk1).dim); 
    else vm=bcc.VM(min(numel(bcc.VM),n)); frois{n}=reshape(round(spm_get_data(vm,pinv(vm.mat)*tXYZ)),bcc.VY1(idxk1).dim(1:3));
    end
    if n==1, nrois=max(frois{n}(:));
    elseif max(frois{n}(:))~=nrois, error('All subject-specific ROI files should contain the same number of ROIs');
    end
    if ~isempty(XM),
        tXM=reshape(double(spm_get_data(XM,pinv(XM.mat)*tXYZ)>0),bcc.VY1(idxk1).dim(1:3));
        frois{n}=tXM.*frois{n};
    end
end

% analysis

Nb=[size(bcc.X,2),neffects];
extname='';

Bplane=nan+zeros([Nb,nrois]);Eplane=nan+zeros(Nb(2),Nb(2),nrois);Zplane=zeros(bcc.n,Nb(2),nrois);Nplane=nan+zeros([bcc.n,nrois]);
fprintf('Performing model estimation...');
for nroi=1:nrois,
    fprintf(1,'.');
    R=zeros([size(bcc.VY1,1),neffects1,neffects2]);
    for n=1:bcc.n
        idx=find(frois{n}==nroi);
        [idx1,idx2,idx3]=ind2sub(size(frois{n}),idx);
        xyz=[idx1,idx2,idx3,ones(numel(idx1),1)]';
        Nplane(n,nroi)=numel(idx);

        idxk=bcc.PV(n,:)>0;
        Y1=spm_get_data(reshape(bcc.VY1(idxk,:),[],1),xyz); % (crossv x neffects1) x nvoxels
        Y2=spm_get_data(reshape(bcc.VY2(idxk,:),[],1),xyz); % (crossv x neffects2) x nvoxels
        Y1(isnan(Y1))=0;
        Y2(isnan(Y2))=0;
        nY1=Y1~=0; % masks out-of-brain voxels
        nY2=Y2~=0;
        mY1=sum(Y1,2)./max(eps,sum(nY1,2));
        mY2=sum(Y2,2)./max(eps,sum(nY2,2));
        Y1=bsxfun(@minus,Y1,mY1);
        Y2=bsxfun(@minus,Y2,mY2);
        Y1(~nY1)=0;
        Y2(~nY2)=0;
        Y1=bsxfun(@rdivide,Y1,max(eps,sqrt(sum(abs(Y1).^2,2))));
        Y2=bsxfun(@rdivide,Y2,max(eps,sqrt(sum(abs(Y2).^2,2))));
        Y1=reshape(Y1,[],neffects1,1,size(Y1,2));
        Y2=reshape(Y2,[],1,neffects2,size(Y2,2));
        R(idxk,:,:)=sum(bsxfun(@times,Y1,Y2),4); % corrs ( crossv x neffects1 x neffects2 )
    end
%     idx=find(frois==nroi);
%     [idx1,idx2,idx3]=ind2sub(size(frois),idx);
%     xyz=[idx1,idx2,idx3,ones(numel(idx1),1)]';
%     Y1=spm_get_data(bcc.VY1,xyz); % (subjects x crossv x neffects1) x nvoxels
%     Y2=spm_get_data(bcc.VY2,xyz); % (subjects x crossv x neffects2) x nvoxels
%     Y1(isnan(Y1))=0;
%     Y2(isnan(Y2))=0;
%     nY1=Y1~=0; % masks out-of-brain voxels
%     nY2=Y2~=0;
%     mY1=sum(Y1,2)./max(eps,sum(nY1,2));
%     mY2=sum(Y2,2)./max(eps,sum(nY2,2));
%     Y1=bsxfun(@minus,Y1,mY1);
%     Y2=bsxfun(@minus,Y2,mY2);
%     Y1(~nY1)=0;
%     Y2(~nY2)=0;
%     Y1=bsxfun(@rdivide,Y1,max(eps,sqrt(sum(abs(Y1).^2,2))));
%     Y2=bsxfun(@rdivide,Y2,max(eps,sqrt(sum(abs(Y2).^2,2))));
%     Y1=reshape(Y1,[],neffects1,1,size(Y1,2));
%     Y2=reshape(Y2,[],1,neffects2,size(Y2,2));
%     R=sum(bsxfun(@times,Y1,Y2),4); % corr
    
    y=bcc.PV*atanh(R(:,:));        % fisher-transformed corrs averaged across crossvalidation partitions

    x=bcc.X;
    [b,ee]=spm_ss_glm('estimate',x,y);
    Bplane(:,:,nroi)=b;
    Eplane(:,:,nroi)=ee;
%     Nplane(nroi)=numel(idx);
    Zplane(:,:,nroi)=y;
end
fprintf(1,'\n');

% save files
bcc.estimate=struct('beta',Bplane,'rss',Eplane,'voxels',Nplane,'y',Zplane,'effect_names',{effect_names}); 
save(fullfile(bcc.swd,['SPM_bcc',extname,'.mat']),'bcc');
disp(['Analysis file saved: ',fullfile(bcc.swd,['SPM_bcc',extname,'.mat'])]);

if 0,% note: this will be removed in the future; change to 1 to create old-format _data.csv file
    fname=['spm_bcc',extname,'_data.csv'];
    fh=fopen(fullfile(bcc.swd,fname),'wt');
    fprintf(fh,'Data\n');
    fprintf(fh,'ROI#,average ROI size,');
    if numel(effect_names)==Nb(2)
        for ns=1:bcc.n,for ne=1:Nb(2),fprintf(fh,'Subject#%d[%s]',ns,effect_names{ne}); if ne<Nb(2)||ns<bcc.n, fprintf(fh,','); else fprintf(fh,'\n'); end; end; end
    else
        for ns=1:bcc.n,for ne=1:Nb(2),fprintf(fh,'Subject#%d[%d]',ns,ne); if ne<Nb(2)||ns<bcc.n, fprintf(fh,','); else fprintf(fh,'\n'); end; end; end
    end
    for nroi=1:nrois, for ns=1:bcc.n, for ne=1:Nb(2), 
        fprintf(fh,'%d,%d,',nroi,round(mean(bcc.estimate.voxels(:,nroi))));
        fprintf(fh,'%f',Zplane(ns,ne,nroi)); if ne<Nb(2)||ns<bcc.n, fprintf(fh,','); else fprintf(fh,'\n'); end; end; end
    end
    fprintf(fh,'\n');
    fclose(fh);
end
fname=['spm_bcc',extname,'_data.EffectSize.csv'];
fh=fopen(fullfile(bcc.swd,fname),'wt');
fprintf(fh,'ROI,Subject,Effects,Fisher transformed correlation coefficients\n');
if numel(effect_names)==Nb(2), for nroi=1:nrois,for ns=1:bcc.n, for ne=1:Nb(2), fprintf(fh,'%d,%d,%s,%f\n',nroi,ns,effect_names{ne},Zplane(ns,ne,nroi)); end; end; end
else for nroi=1:nrois,for ns=1:bcc.n, for ne=1:Nb(2), fprintf(fh,'%d,%d,%d,%f\n',nroi,ns,ne,Zplane(ns,ne,nroi)); end; end; end
end
fclose(fh);

% estimates defined contrasts
bcc=spm_bcc_contrast(bcc);
cd(cwd);

end



