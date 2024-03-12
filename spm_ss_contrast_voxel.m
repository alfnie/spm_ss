function ss=spm_ss_contrast_voxel(ss,ssIc,overwrite)
% SPM_SS_CONTRAST_VOXEL evaluates contrast in subject-specific voxel-based analyses
% 
% ss=spm_ss_contrast_voxel(ss,Ic)
% see SPM_SS_DESIGN, SPM_SS_CONTRAST

if nargin<1, 
    str='Select spm_ss*.mat analysis file';
    disp(str);
    Pdefault='';objname=findobj('tag','spm_ss');if numel(objname)==1,objdata=get(objname,'userdata');if isfield(objdata,'files_spm_ss'),Pdefault=objdata.files_spm_ss;end;end;
    P=spm_select(1,'^SPM_ss.*\.mat$',str,{Pdefault});
    if numel(objname)==1&&~isempty(P),objdata.files_spm_ss=P;set(objname,'userdata',objdata);end;
    load(P);
    ss.swd=fileparts(P);
end
if nargin<2,
    ssIc=1:numel(ss.C);
end
if nargin<3,
    overwrite=1;
end

cwd=pwd;
cd(ss.swd);

for nIc=1:numel(ssIc),
    Ic=ssIc(nIc);
    if overwrite||~isfield(ss,'evaluate')||numel(ss.evaluate)<Ic||isempty(ss.evaluate{Ic}),
        Nh=[size(ss.C(Ic).between,1),size(ss.C(Ic).within,1)];
        Nb=size(ss.X,2);
        extname=['_',ss.type];
        if prod(Nh)==1,fname=['spm_ss',extname,'_T'];else fname=['spm_ss',extname,'_F'];end
        VF=struct('fname',[fname,'_',num2str(Ic,'%04d'),'.img'],...
            'mat',ss.VN(1).mat,...
            'dim',ss.VN(1).dim,...
            'pinfo',[1;0;0],...
            'dt',[spm_type('float32'),spm_platform('bigend')],...
            'descrip','spm_ss (T/F statistics)');
        VF=spm_create_vol(VF);
        VP=struct('fname',['spm_ss',extname,'_P_',num2str(Ic,'%04d'),'.img'],...
            'mat',ss.VN(1).mat,...
            'dim',ss.VN(1).dim,...
            'pinfo',[1;0;0],...
            'dt',[spm_type('float32'),spm_platform('bigend')],...
            'descrip','spm_ss (uncorrected voxel-wise p-value)');
        VP=spm_create_vol(VP);
        VDOF=struct('fname',['spm_ss',extname,'_dof_',num2str(Ic,'%04d'),'.img'],...
            'mat',ss.VN(1).mat,...
            'dim',ss.VN(1).dim,...
            'pinfo',[1;0;0],...
            'dt',[spm_type('float32'),spm_platform('bigend')],...
            'descrip','spm_ss (degrees of freedom)');
        VDOF=spm_create_vol(VDOF);
        VH=struct('fname',['spm_ss',extname,'_con_',num2str(Ic,'%04d'),'.img'],...
            'mat',ss.VN(1).mat,...
            'dim',ss.VN(1).dim,...
            'n',[1,1],...
            'pinfo',[1;0;0],...
            'dt',[spm_type('float32'),spm_platform('bigend')],...
            'descrip','spm_ss (effect sizes contrast estimates)');
        VH=repmat(VH,[prod(Nh),1]);for nh=1:prod(Nh),VH(nh).n=[nh,1];end
        VH=spm_create_vol(VH);
        VPfdr=struct('fname',['spm_ss',extname,'_Pfdr_',num2str(Ic,'%04d'),'.img'],'mat',VP.mat,'dim',VP.dim,'dt',VP.dt);
        
        nsize=ss.VN(1).dim;
        Nx=rank(ss.X);
        Nc0=rank(ss.X*ss.C(Ic).between');
        Ns=rank(ss.C(Ic).within);
        xyz=[1:nsize(1);ones(2,nsize(1))];
        fprintf(['Estimating contrast #',num2str(Ic)]);
        for nz=1:nsize(3),
            fprintf('.');
            Fplane=nan+zeros(nsize(1:2));Pplane=nan+zeros(nsize(1:2));Hplane=nan+zeros([Nh,nsize(1:2)]);DOFplane=nan+zeros(nsize(1:2));
            xyz(3,:)=nz;
            for ny=1:nsize(2),
                xyz(2,:)=ny;
                B=spm_get_data(ss.estimate.BETA,xyz);%,0);
                idx0=find(~any(isnan(B),1));
                if ~isempty(idx0),
                    B=B(:,idx0);
                    W=spm_get_data(ss.estimate.WHITENING,xyz(:,idx0));%,0);
                    E=spm_get_data(ss.estimate.RSS,xyz(:,idx0));%,0);
                    dof=(sum(W,1).^2)./max(eps,sum(W.^2,1));%Welch–Satterthwaite
                    idx=find(~any(isnan(E),1) & dof-Nx>2);
                    for nidx=1:numel(idx),
                        nx=idx(nidx);nxx=idx0(nx);
                        x=ss.X.*W(:,nx+zeros(1,Nb));
                        b=reshape(B(:,nx),Nb,[]);
                        rss=reshape(E(:,nx),size(b,2),size(b,2));
                        [h,f,p,edof]=spm_ss_glm('evaluate',x,b,rss,ss.C(Ic).between,ss.C(Ic).within,dof(nx),Nx,Nc0,Nh,Ns); % general linear model
                        dof(nx)=edof(end);
                        Hplane(:,:,nxx,ny)=h;
                        Fplane(nxx,ny)=f;
                        Pplane(nxx,ny)=p;
                    end
                    DOFplane(idx0,ny)=dof';
                end
            end
            VF=spm_write_plane(VF,Fplane,nz);
            VP=spm_write_plane(VP,Pplane,nz);
            VDOF=spm_write_plane(VDOF,DOFplane,nz);
            nh=1;for nh1=1:Nh(1),for nh2=1:Nh(2),VH(nh)=spm_write_plane(VH(nh),shiftdim(Hplane(nh1,nh2,:,:),2),nz);nh=nh+1;end;end
        end
        p=spm_read_vols(VP);
        Pfdr=p;Pfdr(:)=spm_ss_fdr(p(:));
        spm_write_vol(VPfdr,Pfdr);
        fprintf('\nDone\n');
        
        disp(['contrast #: ',num2str(Ic)]);
        disp(['contrast volume : ',fullfile(ss.swd,VH(1).fname),' - ',num2str(prod(Nh)),' volume(s)']);
        if prod(Nh)==1,   disp(['T statistics    : ',fullfile(ss.swd,VF(1).fname)]);
        else        disp(['F statistics    : ',fullfile(ss.swd,VF(1).fname)]);end
        disp(['dof             : ',fullfile(ss.swd,VDOF(1).fname)]);
        disp(['p values        : ',fullfile(ss.swd,VP(1).fname)]);
        disp(['p-fdr values    : ',fullfile(ss.swd,VPfdr(1).fname)]);
        
        ss.evaluate{Ic}=struct('name',ss.C(Ic).name,'CON',VH,'T',VF,'DOF',VDOF,'P',VP,'PFDR',VPfdr);
        save(fullfile(ss.swd,['SPM_ss',extname,'.mat']),'ss');
    end
end

cd(cwd);
    
end

