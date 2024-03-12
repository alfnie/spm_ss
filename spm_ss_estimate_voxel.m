function ss=spm_ss_estimate_voxel(ss)
% SPM_SS_ESTIMATE_VOXEL subject-specific voxel-based model estimation
% 
% ss=spm_ss_estimate_voxel(ss)
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
if ~isempty(ss.ExplicitMasking),XM=spm_read_vols(spm_vol(ss.ExplicitMasking))>0;else XM=1;end

% creates transformed measures
fprintf('Creating transformed files...');
neffects=size(ss.EffectOfInterest{1},1);
Mplane=[];nm=1;
for n=1:ss.n,
    p1=0;
    for m=1:numel(ss.Localizer{n}),
        [pth1,nm1,ext1,num1]=spm_fileparts(ss.Localizer{n}{m});
        Nvolume=['t',nm1,ext1,num1];
        fprintf(1,'.');
        a1=spm_vol(ss.Localizer{n}{m});
        b1=spm_read_vols(a1);
        if numel(XM)>1&&any(size(XM)~=size(b1)),error('mismatched volume dimensions between explicit mask and localizer image file'); end
        b1=XM.*b1;
        Z1=isnan(b1)|(b1==0);
        b1(Z1)=0;
        Mplane(nm)=sum(b1(:));nm=nm+1;
        p1=p1+b1/numel(ss.Localizer{n});
        c1=struct('fname',Nvolume,'descrip','Masked','mat',a1.mat,'dim',a1.dim,'dt',[spm_type('float32'),spm_platform('bigend')]);
        cd(pth1);   c1=spm_write_vol(c1,b1);
        spm_smooth(a1,['s',Nvolume],ss.smooth*[1,1,1]);
        for ne=1:neffects,
            [pth2,nm2,ext2,num2]=spm_fileparts(ss.EffectOfInterest{n}{ne,m});
            Yvolume=['t',nm2,ext2,num2];
            a2=spm_vol(ss.EffectOfInterest{n}{ne,m});
            b2=spm_read_vols(a2);
            if any(any(a1.mat~=a2.mat))||any(any(a1.dim~=a2.dim)), error(['Mismatched orientations/dimensions between localizer volume and contrast of interest volume for subject ',num2str(n)]); end
            Z2=isnan(b2);
            b2(Z1|Z2)=0;
            d1=struct('fname',Yvolume,'descrip','Masked','mat',a2.mat,'dim',a2.dim,'dt',[spm_type('float32'),spm_platform('bigend')]);
            cd(pth2);   d1=spm_write_vol(d1,b1.*b2);
            spm_smooth(d1,['s',Yvolume],ss.smooth*[1,1,1]);
            e1=spm_vol(['s',Yvolume]);c1=spm_read_vols(e1);c1(Z2&~Z1)=nan;spm_write_vol(e1,c1); % keeps nan
        end
    end
    if n==1,
        ss.PL=['Localizer',ext1];
        e0=struct('fname',ss.PL,'descrip','spm_ss (localizer mask for each subject)','mat',a1.mat,'dim',a1.dim,'n',[1,1],'pinfo',[1;0;0],'dt',[spm_type('float32'),spm_platform('bigend')]);
        e0=repmat(e0,[ss.n,1]);for nb=1:ss.n,e0(nb).n=[nb,1];end
        cd(ss.swd); e0=spm_create_vol(e0);
    end
    cd(ss.swd); spm_write_vol(e0(n),p1);
end

neffects=size(ss.EffectOfInterest{1},1);
k=0;
for ne=1:neffects,
    for n=1:ss.n,
        for m=1:numel(ss.Localizer{n}),
            [pth2,nm2,ext2,num2]=spm_fileparts(ss.EffectOfInterest{n}{ne,m});
            Yvolume=['t',nm2,ext2,num2];
            k=k+1;
            ss.PY{k}=fullfile(pth2,['s',Yvolume]);
        end
    end
end
k=0;
ss.PV=[];
for n=1:ss.n,
    for m=1:numel(ss.Localizer{n}),
        [pth1,nm1,ext1,num1]=spm_fileparts(ss.Localizer{n}{m});
        Nvolume=['t',nm1,ext1,num1];
        k=k+1;
        ss.PN{k}=fullfile(pth1,['s',Nvolume]);
        ss.PV(n,k)=1/numel(ss.Localizer{n});
    end
end
fprintf(1,'\n');
cd(ss.swd);

% analysis
ss.VN=spm_vol(char(ss.PN));
ss.VY=spm_vol(char(ss.PY));

Nb=[size(ss.X,2),neffects];
extname=['_',ss.type];
VB=struct('fname',['spm_ss',extname,'_beta.img'],...
    'mat',ss.VN(1).mat,...
    'dim',ss.VN(1).dim,...
    'n',[1,1],...
    'pinfo',[1;0;0],...
    'dt',[spm_type('float32'),spm_platform('bigend')],...
    'descrip','spm_ss (effect sizes parameter estimates)');
VB=repmat(VB,[prod(Nb),1]);for nb=1:prod(Nb),VB(nb).n=[nb,1];end
VB=spm_create_vol(VB);
VC=struct('fname',['spm_ss',extname,'_whitening.img'],...
    'mat',ss.VN(1).mat,...
    'dim',ss.VN(1).dim,...
    'n',[1,1],...
    'pinfo',[1;0;0],...
    'dt',[spm_type('float32'),spm_platform('bigend')],...
    'descrip','spm_ss (whitening vector)');
VC=repmat(VC,[ss.n,1]);for nc=1:ss.n,VC(nc).n=[nc,1];end
VC=spm_create_vol(VC);
VE=struct('fname',['spm_ss',extname,'_rss.img'],...
    'mat',ss.VN(1).mat,...
    'dim',ss.VN(1).dim,...
    'pinfo',[1;0;0],...
    'dt',[spm_type('float32'),spm_platform('bigend')],...
    'descrip','spm_ss (residual sum squares)');
VE=repmat(VE,[Nb(2)*Nb(2),1]);for nb=1:Nb(2)*Nb(2),VE(nb).n=[nb,1];end
VE=spm_create_vol(VE);
VO=struct('fname',['spm_ss',extname,'_overlap.img'],...
    'mat',ss.VN(1).mat,...
    'dim',ss.VN(1).dim,...
    'pinfo',[1;0;0],...
    'dt',[spm_type('float32'),spm_platform('bigend')],...
    'descrip','spm_ss (proportion overlap)');
VO=spm_create_vol(VO);
VS=struct('fname',['spm_ss',extname,'_bw.img'],...
    'mat',ss.VN(1).mat,...
    'dim',ss.VN(1).dim,...
    'pinfo',[1;0;0],...
    'dt',[spm_type('float32'),spm_platform('bigend')],...
    'descrip','spm_ss (sigma ratio)');
VS=spm_create_vol(VS);

nsize=ss.VN(1).dim;
xyz=[1:nsize(1);ones(2,nsize(1))];
fprintf('Performing model estimation...');
for nz=1:nsize(3),
    fprintf(1,'.');
    Bplane=nan+zeros([Nb,nsize(1:2)]);Cplane=zeros([ss.n,nsize(1:2)]);Eplane=nan+zeros([Nb(2),Nb(2),nsize(1:2)]);Oplane=nan+zeros(nsize(1:2));Splane=nan+zeros(nsize(1:2));
    xyz(3,:)=nz;
    for ny=1:nsize(2),
        xyz(2,:)=ny;        
        Y=spm_get_data(ss.VY(:),xyz);
        N=spm_get_data(ss.VN(:),xyz);
        
        Y=reshape(Y./max(eps,repmat(N,[Nb(2),1])),[size(N,1),Nb(2),size(xyz,2)]);
        Y=reshape(ss.PV*Y(:,:),[ss.n,Nb(2),size(xyz,2)]);
        N=1./(ss.PV*(1./max(eps,N)));
        sN=mean(N>1e-4,1).*(~any(isnan(N),1));
        idx=find(sN>0);
        for nidx=1:numel(idx),
            nx=idx(nidx);
            if strcmpi(ss.estimation,'ols')
                iC=double(N(:,nx)>1e-4);%handles missing-data
                bw=nan;%(not estimated; equivalent to bw=inf)
            else
                n=N(:,nx);
                y=Y(:,:,nx).*sqrt(n(:,ones(1,Nb(2))));
                x=ss.X.*sqrt(n(:,ones(1,Nb(1))));
                e=Y(:,:,nx)-ss.X*(pinv(x'*x)*(x'*y));
                [nill,iC,nill,bw]=spm_ss_fls({e,n});%covariance estimation
            end
            y=Y(:,:,nx).*iC(:,ones(1,Nb(2)));%whitening
            x=ss.X.*iC(:,ones(1,Nb(1)));
            [b,ee]=spm_ss_glm('estimate',x,y);
            Bplane(:,:,nx,ny)=b;
            Cplane(:,nx,ny)=iC;
            Eplane(:,:,nx,ny)=ee;
            Splane(nx,ny)=bw;
        end
        Oplane(:,ny)=sN';
    end
    nb=1;for nb1=1:Nb(1),for nb2=1:Nb(2),VB(nb)=spm_write_plane(VB(nb),shiftdim(Bplane(nb1,nb2,:,:),2),nz);nb=nb+1;end;end
    for nn=1:ss.n,VC(nn)=spm_write_plane(VC(nn),shiftdim(Cplane(nn,:,:),1),nz);end
    nb=1;for nb1=1:Nb(2),for nb2=1:Nb(2),VE(nb)=spm_write_plane(VE(nb),shiftdim(Eplane(nb1,nb2,:,:),2),nz);nb=nb+1;end;end
    VO=spm_write_plane(VO,Oplane,nz);
    VS=spm_write_plane(VS,Splane,nz);
end
fprintf(1,'\n');

disp(['created beta volume       : ',fullfile(ss.swd,VB(1).fname),' - ',num2str(prod(Nb)),' volume(s)']); 
disp(['created whitening volume  : ',fullfile(ss.swd,VC(1).fname),' - ',num2str(ss.n),' volume(s)']); 
disp(['created rss volume        : ',fullfile(ss.swd,VE(1).fname),' - ',num2str(Nb(2)*Nb(2)),' volume(s)']); 
disp(['created overlap volume    : ',fullfile(ss.swd,VO(1).fname)]); 

ss.estimate=struct('BETA',VB,'RSS',VE,'WHITENING',VC,'OVERLAP',VO,'qa',Mplane);
save(fullfile(ss.swd,['SPM_ss',extname,'.mat']),'ss');
disp(['Analysis file saved: ',fullfile(ss.swd,['SPM_ss',extname,'.mat'])]);

fname=['spm_ss',extname,'_data.csv'];
fh=fopen(fullfile(ss.swd,fname),'wt');
fprintf(fh,'\nquality control (localizer mask sizes)\n');
fprintf(fh,'Subject#,Session/Partition#,filename,localizer mask size\n');
nidxs=zeros(ss.n,1);
for nf=1:numel(ss.VN,1),
    [nill,idxs]=max(ss.PV(:,nf));
    nidxs(idxs)=nidxs(idxs)+1;
    fprintf(fh,'%d,%d,%s,%d\n',idxs,nidxs(idxs),ss.PN{nf},round(ss.estimate.qa(nf)));
end
fclose(fh);


% estimates defined contrasts
spm_ss_contrast_voxel(ss);
cd(cwd);

end
