function ss=spm_ss_contrast_ROI(ss,ssIc,overwrite)
% SPM_SS_CONTRAST_ROI evaluates contrast in subject-specific ROI-based analyses
% 
% ss=spm_ss_contrast_ROI(ss,Ic)
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
ss.VM=spm_vol(char(ss.PM));
% [XYZ{1},XYZ{2},XYZ{3}]=ndgrid(1:ss.VN(1).dim(1),1:ss.VN(1).dim(2),1:ss.VN(1).dim(3));
% XYZ=reshape(cat(4,XYZ{:}),[],3)';XYZ=ss.VN(1).mat*cat(1,XYZ,ones(1,size(XYZ,2)));
% frois=reshape(round(spm_get_data(ss.VM,pinv(ss.VM.mat)*XYZ))',[ss.VN(1).dim(1:3),numel(ss.VM)]);
% nrois=max(frois(:));
for n=1 % frois in subject-space (first-subject)
    idxk=ss.PV(n,:)>0;
    idxk1=find(idxk,1);
    [XYZ{1},XYZ{2},XYZ{3}]=ndgrid(1:ss.VN(idxk1).dim(1),1:ss.VN(idxk1).dim(2),1:ss.VN(idxk1).dim(3));
    tXYZ=reshape(cat(4,XYZ{:}),[],3)';tXYZ=ss.VN(idxk1).mat*cat(1,tXYZ,ones(1,size(tXYZ,2))); % voxel-coordinates for first functional volume for subject n
    vm=ss.VM(min(numel(ss.VM),n));
    frois=reshape(round(spm_get_data(vm,pinv(vm.mat)*tXYZ)),ss.VN(idxk1).dim(1:3));
    nrois=max(frois(:));
end
if ~isfield(ss,'VM_roinames')
    [ss.VM_roinames,ss.VM_roiids]=spm_ss_roilabels(ss.VM(1).fname);
    if isempty(ss.VM_roinames), ss.VM_roinames=arrayfun(@num2str,1:nrois,'uni',0); ss.VM_roiids=1:nrois;
    else nrois=numel(ss.VM_roiids);
    end
end
        
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
        try, spm_unlink(VH.fname); end
        VH=repmat(VH,[prod(Nh),1]);for nh=1:prod(Nh),VH(nh).n=[nh,1];end
        VH=spm_create_vol(VH);
        VPfdr=struct('fname',['spm_ss',extname,'_Pfdr_',num2str(Ic,'%04d'),'.img'],'mat',VP.mat,'dim',VP.dim,'dt',VP.dt);
        
        Nx=rank(ss.X);
        Nc0=rank(ss.X*ss.C(Ic).between');
        Ns=rank(ss.C(Ic).within);
        Fplane=nan+zeros(1,nrois);Pplane=nan+zeros(1,nrois);Hplane=nan+zeros([Nh,nrois]);DOFplane=zeros(1,nrois);Eplane=nan+zeros([Nh,nrois]);
        fprintf(['Estimating contrast #',num2str(Ic)]);
        for nroi=1:nrois,
            fprintf('.');
            B=ss.estimate.beta(:,:,nroi);
            W=ss.estimate.whitening(:,nroi);
            E=ss.estimate.rss(:,:,nroi);
            dof=(sum(W,1).^2)./max(eps,sum(W.^2,1));%Welch–Satterthwaite
            if ~any(any(isnan(B),1),2) && ~any(any(isnan(E))),% && dof-Nx>2,
                x=ss.X.*W(:,ones(1,Nb));
                b=B;
                rss=E;
                [h,f,p,edof,stderr]=spm_ss_glm('evaluate',x,b,rss,ss.C(Ic).between,ss.C(Ic).within,dof,Nx,Nc0,Nh,Ns); % general linear model
                dof=edof(end);
                Hplane(:,:,nroi)=h;
                Fplane(nroi)=f;
                Pplane(nroi)=p;
                Eplane(:,:,nroi)=stderr;
            end
            DOFplane(nroi)=dof;
        end
        PFDRplane=Pplane;
        PFDRplane(ss.estimate.overlap<ss.overlap_thr_roi)=nan;
        PFDRplane(:)=spm_ss_fdr(PFDRplane(:));
        fprintf(1,'\nDone\n');
        
        nh=1;for nh1=1:Nh(1),for nh2=1:Nh(2),z=nan+zeros(size(frois));for nroi=1:nrois,z(frois==ss.VM_roiids(nroi))=Hplane(nh1,nh2,nroi);end; spm_write_vol(VH(nh),z);nh=nh+1;end;end
        z=nan+zeros(size(frois));for nroi=1:nrois,z(frois==ss.VM_roiids(nroi))=Fplane(nroi);end; spm_write_vol(VF,z);
        z=nan+zeros(size(frois));for nroi=1:nrois,z(frois==ss.VM_roiids(nroi))=Pplane(nroi);end; spm_write_vol(VP,z);
        z=nan+zeros(size(frois));for nroi=1:nrois,z(frois==ss.VM_roiids(nroi))=DOFplane(nroi);end; spm_write_vol(VDOF,z);
        z=nan+zeros(size(frois));for nroi=1:nrois,z(frois==ss.VM_roiids(nroi))=PFDRplane(nroi);end; spm_write_vol(VPfdr,z);
        
        %     disp(['contrast #: ',num2str(Ic)]);
        %     disp(['contrast volume : ',fullfile(ss.swd,VH(1).fname),' - ',num2str(Nh),' volume(s)']);
        %     if Nh==1,   disp(['T statistics    : ',fullfile(ss.swd,VF(1).fname)]);
        %     else        disp(['F statistics    : ',fullfile(ss.swd,VF(1).fname)]);end
        %     disp(['dof             : ',fullfile(ss.swd,VDOF(1).fname)]);
        %     disp(['p values        : ',fullfile(ss.swd,VP(1).fname)]);
        %     disp(['p-fdr values    : ',fullfile(ss.swd,VPfdr(1).fname)]);
        
        ss.evaluate{Ic}=struct('name',ss.C(Ic).name,'CON',VH,'T',VF,'DOF',VDOF,'P',VP,'PFDR',VPfdr,'con',Hplane,'stderr',Eplane,'t',Fplane,'dof',DOFplane,'p',Pplane,'pfdr',PFDRplane);
        save(fullfile(ss.swd,['SPM_ss',extname,'.mat']),'ss');
        
        if 0,% note: this will be removed in the future; change to 1 to create old-format _data.csv file
            fname=['spm_ss',extname,'_results_',num2str(Ic,'%04d'),'.csv'];
            fh=fopen(fullfile(ss.swd,fname),'wt');
            fprintf(fh,'Contrast name:,%s\n',ss.C(Ic).name);
            fprintf(fh,'Within-subjects contrast matrix:,');for nc1=1:size(ss.C(Ic).within,1),for nc2=1:size(ss.C(Ic).within,2),fprintf(fh,'c[%d;%d],',nc1,nc2);end;end;fprintf(fh,'\n');
            fprintf(fh,',');for nc1=1:size(ss.C(Ic).within,1),for nc2=1:size(ss.C(Ic).within,2),fprintf(fh,'%f,',ss.C(Ic).within(nc1,nc2));end;end;fprintf(fh,'\n');
            
            fprintf(fh,'Between-subjects contrast matrix:,');for nc1=1:size(ss.C(Ic).between,1),for nc2=1:size(ss.C(Ic).between,2),fprintf(fh,'c[%d;%d],',nc1,nc2);end;end;fprintf(fh,'\n');
            fprintf(fh,',');for nc1=1:size(ss.C(Ic).between,1),for nc2=1:size(ss.C(Ic).between,2),fprintf(fh,'%f,',ss.C(Ic).between(nc1,nc2));end;end;fprintf(fh,'\n');
            
            fprintf(fh,'Overlap threshold:,%f\n',ss.overlap_thr_roi);
            fprintf(fh,'\nResults\n');
            fprintf(fh,'ROI#,average ROI size,average localizer mask size,inter-subject overlap,');
            for nh1=1:Nh(1),for nh2=1:Nh(2),fprintf(fh,'h[%d;%d],',nh1,nh2);end; end;
            for nh1=1:Nh(1),for nh2=1:Nh(2),fprintf(fh,'stderr[%d;%d],',nh1,nh2);end; end;
            %if prod(Nh)==1,fprintf(fh,'stderr,'); end;
            fprintf(fh,'T/F,dof,p,p-fdr\n');
            for nroi=1:nrois,
                if ss.estimate.overlap(nroi)>=ss.overlap_thr_roi,
                    fprintf(fh,'%d,%d,%d,%f,',nroi,round(mean(ss.estimate.voxels(:,nroi))),round(mean(ss.estimate.voxels(:,nroi))*ss.estimate.coverage(nroi)),ss.estimate.overlap(nroi));
                    for nh1=1:Nh(1),for nh2=1:Nh(2),fprintf(fh,'%f,',ss.evaluate{Ic}.con(nh1,nh2,nroi));end; end;
                    for nh1=1:Nh(1),for nh2=1:Nh(2),fprintf(fh,'%f,',ss.evaluate{Ic}.stderr(nh1,nh2,nroi));end; end;
                    %if prod(Nh)==1,fprintf(fh,'%f,',ss.evaluate{Ic}.con(1,1,nroi)./ss.evaluate{Ic}.t(nroi));end;
                    fprintf(fh,'%f,%f,%f,%f\n',ss.evaluate{Ic}.t(nroi),ss.evaluate{Ic}.dof(nroi),ss.evaluate{Ic}.p(nroi),ss.evaluate{Ic}.pfdr(nroi));
                end
            end
            fclose(fh);
        end
        if size(ss.EffectOfInterest_contrasts,2)==Nh(2)&&isequal(ss.C(Ic).within,eye(Nh(2))), effect_names=ss.EffectOfInterest_contrasts;
        else effect_names=arrayfun(@num2str,1:Nh(2),'uni',0);
        end
        fname=['spm_ss',extname,'_results_',num2str(Ic,'%04d'),'.Descrip.csv'];
        fh=fopen(fullfile(ss.swd,fname),'wt');
        fprintf(fh,'Contrast name:,%s\n',ss.C(Ic).name);
        fprintf(fh,'Within-subjects contrast matrix:,');for nc1=1:size(ss.C(Ic).within,1),for nc2=1:size(ss.C(Ic).within,2),fprintf(fh,'c[%d;%d],',nc1,nc2);end;end;fprintf(fh,'\n');
        fprintf(fh,',');for nc1=1:size(ss.C(Ic).within,1),for nc2=1:size(ss.C(Ic).within,2),fprintf(fh,'%f,',ss.C(Ic).within(nc1,nc2));end;end;fprintf(fh,'\n');
        fprintf(fh,'Between-subjects contrast matrix:,');for nc1=1:size(ss.C(Ic).between,1),for nc2=1:size(ss.C(Ic).between,2),fprintf(fh,'c[%d;%d],',nc1,nc2);end;end;fprintf(fh,'\n');
        fprintf(fh,',');for nc1=1:size(ss.C(Ic).between,1),for nc2=1:size(ss.C(Ic).between,2),fprintf(fh,'%f,',ss.C(Ic).between(nc1,nc2));end;end;fprintf(fh,'\n');
        fprintf(fh,'Overlap threshold:,%f\n',ss.overlap_thr_roi);
        fclose(fh);
        fname=['spm_ss',extname,'_results_',num2str(Ic,'%04d'),'.Stats.csv'];
        fh=fopen(fullfile(ss.swd,fname),'wt');
        fprintf(fh,'ROI,average ROI size,average localizer mask size,inter-subject overlap,');
        for nh1=1:Nh(1),for nh2=1:Nh(2),fprintf(fh,'mean[%d;%s],',nh1,effect_names{nh2});end; end;
        for nh1=1:Nh(1),for nh2=1:Nh(2),fprintf(fh,'stderr[%d;%s],',nh1,effect_names{nh2});end; end;
        %if prod(Nh)==1,fprintf(fh,'stderr,'); end;
        fprintf(fh,'T/F,dof,p,p-fdr\n');
        for nroi=1:nrois,
            if ss.estimate.overlap(nroi)>=ss.overlap_thr_roi,
                fprintf(fh,'%s,%d,%d,%f,',ss.VM_roinames{nroi},round(mean(ss.estimate.voxels(:,nroi))),round(mean(ss.estimate.voxels(:,nroi))*ss.estimate.coverage(nroi)),ss.estimate.overlap(nroi));
                for nh1=1:Nh(1),for nh2=1:Nh(2),fprintf(fh,'%f,',ss.evaluate{Ic}.con(nh1,nh2,nroi));end; end;
                for nh1=1:Nh(1),for nh2=1:Nh(2),fprintf(fh,'%f,',ss.evaluate{Ic}.stderr(nh1,nh2,nroi));end; end;
                %if prod(Nh)==1,fprintf(fh,'%f,',ss.evaluate{Ic}.con(1,1,nroi)./ss.evaluate{Ic}.t(nroi));end;
                fprintf(fh,'%f,%f,%f,%f\n',ss.evaluate{Ic}.t(nroi),ss.evaluate{Ic}.dof(nroi),ss.evaluate{Ic}.p(nroi),ss.evaluate{Ic}.pfdr(nroi));
            end
        end
        fclose(fh);
    end
end

cd(cwd);
    
end

