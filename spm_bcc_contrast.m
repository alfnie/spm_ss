function bcc=spm_bcc_contrast(bcc,ssIc,overwrite)
% SPM_BCC_CONTRAST evaluates contrast in between-conditions-correlation analyses
% 
% bcc=spm_bcc_contrast(bcc,Ic)
% see SPM_BCC_DESIGN, SPM_BCC_ESTIMATE

if nargin<1, 
    str='Select spm_bcc*.mat analysis file';
    disp(str);
    Pdefault='';objname=findobj('tag','spm_bcc');if numel(objname)==1,objdata=get(objname,'userdata');if isfield(objdata,'files_spm_bcc'),Pdefault=objdata.files_spm_bcc;end;end;
    P=spm_select(1,'^SPM_bcc.*\.mat$',str,{Pdefault});
    if numel(objname)==1&&~isempty(P),objdata.files_spm_bcc=P;set(objname,'userdata',objdata);end;
    load(P);
    bcc.swd=fileparts(P);
end
if nargin<2,
    ssIc=1:numel(bcc.C);
end
if nargin<3,
    overwrite=1;
end

cwd=pwd;
cd(bcc.swd);
if isempty(bcc.PM), bcc.VM=[];
else bcc.VM=spm_vol(char(bcc.PM));
end
% bcc.VM=spm_vol(bcc.PM);
% [XYZ{1},XYZ{2},XYZ{3}]=ndgrid(1:bcc.VY1(1).dim(1),1:bcc.VY1(1).dim(2),1:bcc.VY1(1).dim(3));
% XYZ=reshape(cat(4,XYZ{:}),[],3)';XYZ=bcc.VY1(1).mat*cat(1,XYZ,ones(1,size(XYZ,2)));
% if isempty(bcc.VM), frois=ones(bcc.VY1(1).dim); 
% else                frois=reshape(round(spm_get_data(bcc.VM,pinv(bcc.VM.mat)*XYZ))',[bcc.VY1(1).dim(1:3),numel(bcc.VM)]);
% end
% nrois=max(frois(:));
for n=1 % frois in subject-space (first-subject only)
    idxk=bcc.PV(n,:)>0;
    idxk1=find(idxk,1);
    [XYZ{1},XYZ{2},XYZ{3}]=ndgrid(1:bcc.VY1(idxk1).dim(1),1:bcc.VY1(idxk1).dim(2),1:bcc.VY1(idxk1).dim(3));
    tXYZ=reshape(cat(4,XYZ{:}),[],3)';tXYZ=bcc.VY1(idxk1).mat*cat(1,tXYZ,ones(1,size(tXYZ,2))); % voxel-coordinates for first functional volume for subject n
    if isempty(bcc.VM), frois{n}=ones(bcc.VY1(idxk1).dim); 
    else vm=bcc.VM(min(numel(bcc.VM),n)); frois{n}=reshape(round(spm_get_data(vm,pinv(vm.mat)*tXYZ)),bcc.VY1(idxk1).dim(1:3));
    end
    nrois=max(frois{n}(:));
end

for nIc=1:numel(ssIc),
    Ic=ssIc(nIc);
    if overwrite||~isfield(bcc,'evaluate')||numel(bcc.evaluate)<Ic||isempty(bcc.evaluate{Ic}),
        Nh=[size(bcc.C(Ic).between,1),size(bcc.C(Ic).within,1)];
        Nb=size(bcc.X,2);
        extname='';%['_',bcc.type];
        
        Nx=rank(bcc.X);
        Nc0=rank(bcc.X*bcc.C(Ic).between');
        Ns=rank(bcc.C(Ic).within);
        Fplane=nan+zeros(1,nrois);Pplane=nan+zeros(1,nrois);Hplane=nan+zeros([Nh,nrois]);DOFplane=zeros(1,nrois);Eplane=nan+zeros([Nh,nrois]);
        fprintf(['Estimating contrast #',num2str(Ic)]);
        for nroi=1:nrois,
            fprintf('.');
            B=bcc.estimate.beta(:,:,nroi);
            %W=bcc.estimate.whitening(:,nroi);
            E=bcc.estimate.rss(:,:,nroi);
            dof=size(bcc.X,1);%(sum(W,1).^2)./max(eps,sum(W.^2,1));%Welch–Satterthwaite
            if ~any(any(isnan(B),1),2) && ~any(any(isnan(E))),% && dof-Nx>2,
                x=bcc.X;%.*W(:,ones(1,Nb));
                b=B;
                rss=E;
                [h,f,p,edof,stderr]=spm_ss_glm('evaluate',x,b,rss,bcc.C(Ic).between,bcc.C(Ic).within,dof,Nx,Nc0,Nh,Ns); % general linear model
                dof=edof(end);
                Hplane(:,:,nroi)=h;
                Fplane(nroi)=f;
                Pplane(nroi)=p;
                Eplane(:,:,nroi)=stderr;
            end
            DOFplane(nroi)=dof;
        end
        PFDRplane=Pplane;
        %PFDRplane(bcc.estimate.overlap<bcc.overlap_thr_roi)=nan;
        PFDRplane(:)=spm_ss_fdr(PFDRplane(:));
        fprintf(1,'\nDone\n');
        
        bcc.evaluate{Ic}=struct('name',bcc.C(Ic).name,'con',Hplane,'stderr',Eplane,'t',Fplane,'dof',DOFplane,'p',Pplane,'pfdr',PFDRplane);
        save(fullfile(bcc.swd,['SPM_bcc',extname,'.mat']),'bcc');
        
        if 0,% note: this will be removed in the future; change to 1 to create old-format _data.csv file
            fname=['spm_bcc',extname,'_results_',num2str(Ic,'%04d'),'.csv'];
            fh=fopen(fullfile(bcc.swd,fname),'wt');
            fprintf(fh,'Contrast name:,%s\n',bcc.C(Ic).name);
            fprintf(fh,'Within-subjects contrast matrix:,');for nc1=1:size(bcc.C(Ic).within,1),for nc2=1:size(bcc.C(Ic).within,2),fprintf(fh,'c[%d;%d],',nc1,nc2);end;end;fprintf(fh,'\n');
            fprintf(fh,',');for nc1=1:size(bcc.C(Ic).within,1),for nc2=1:size(bcc.C(Ic).within,2),fprintf(fh,'%f,',bcc.C(Ic).within(nc1,nc2));end;end;fprintf(fh,'\n');
            fprintf(fh,'Between-subjects contrast matrix:,');for nc1=1:size(bcc.C(Ic).between,1),for nc2=1:size(bcc.C(Ic).between,2),fprintf(fh,'c[%d;%d],',nc1,nc2);end;end;fprintf(fh,'\n');
            fprintf(fh,',');for nc1=1:size(bcc.C(Ic).between,1),for nc2=1:size(bcc.C(Ic).between,2),fprintf(fh,'%f,',bcc.C(Ic).between(nc1,nc2));end;end;fprintf(fh,'\n');
            %fprintf(fh,'Overlap threshold:,%f\n',bcc.overlap_thr_roi);
            fprintf(fh,'\nResults\n');
            fprintf(fh,'ROI#,average ROI size,');
            for nh1=1:Nh(1),for nh2=1:Nh(2),fprintf(fh,'h[%d;%d],',nh1,nh2);end; end;
            for nh1=1:Nh(1),for nh2=1:Nh(2),fprintf(fh,'stderr[%d;%d],',nh1,nh2);end; end;
            %if prod(Nh)==1,fprintf(fh,'stderr,'); end;
            fprintf(fh,'T/F,dof,p,p-fdr\n');
            for nroi=1:nrois,
                fprintf(fh,'%d,%d,',nroi,round(mean(bcc.estimate.voxels(:,nroi))));
                for nh1=1:Nh(1),for nh2=1:Nh(2),fprintf(fh,'%f,',bcc.evaluate{Ic}.con(nh1,nh2,nroi));end; end;
                for nh1=1:Nh(1),for nh2=1:Nh(2),fprintf(fh,'%f,',bcc.evaluate{Ic}.stderr(nh1,nh2,nroi));end; end;
                %if prod(Nh)==1,fprintf(fh,'%f,',bcc.evaluate{Ic}.con(1,1,nroi)./bcc.evaluate{Ic}.t(nroi));end;
                fprintf(fh,'%f,%f,%f,%f\n',bcc.evaluate{Ic}.t(nroi),bcc.evaluate{Ic}.dof(nroi),bcc.evaluate{Ic}.p(nroi),bcc.evaluate{Ic}.pfdr(nroi));
            end
            fclose(fh);
            fprintf('Contrast file saved: %s\n',fullfile(bcc.swd,fname));
        end
        fname=['spm_bcc',extname,'_results_',num2str(Ic,'%04d'),'.descrip.csv'];
        fh=fopen(fullfile(bcc.swd,fname),'wt');
        fprintf(fh,'Contrast name:,%s\n',bcc.C(Ic).name);
        fprintf(fh,'Within-subjects contrast matrix:,');for nc1=1:size(bcc.C(Ic).within,1),for nc2=1:size(bcc.C(Ic).within,2),fprintf(fh,'c[%d;%d],',nc1,nc2);end;end;fprintf(fh,'\n');
        fprintf(fh,',');for nc1=1:size(bcc.C(Ic).within,1),for nc2=1:size(bcc.C(Ic).within,2),fprintf(fh,'%f,',bcc.C(Ic).within(nc1,nc2));end;end;fprintf(fh,'\n');
        fprintf(fh,'Between-subjects contrast matrix:,');for nc1=1:size(bcc.C(Ic).between,1),for nc2=1:size(bcc.C(Ic).between,2),fprintf(fh,'c[%d;%d],',nc1,nc2);end;end;fprintf(fh,'\n');
        fprintf(fh,',');for nc1=1:size(bcc.C(Ic).between,1),for nc2=1:size(bcc.C(Ic).between,2),fprintf(fh,'%f,',bcc.C(Ic).between(nc1,nc2));end;end;fprintf(fh,'\n');
        %fprintf(fh,'Overlap threshold:,%f\n',bcc.overlap_thr_roi);
        fname=['spm_bcc',extname,'_results_',num2str(Ic,'%04d'),'.Stats.csv'];
        fh=fopen(fullfile(bcc.swd,fname),'wt');
        fprintf(fh,'ROI,average ROI size,');
        for nh1=1:Nh(1),for nh2=1:Nh(2),fprintf(fh,'h[%d;%d],',nh1,nh2);end; end;
        for nh1=1:Nh(1),for nh2=1:Nh(2),fprintf(fh,'stderr[%d;%d],',nh1,nh2);end; end;
        %if prod(Nh)==1,fprintf(fh,'stderr,'); end;
        fprintf(fh,'T/F,dof,p,p-fdr\n');
        for nroi=1:nrois,
            fprintf(fh,'%d,%d,',nroi,round(mean(bcc.estimate.voxels(:,nroi))));
            for nh1=1:Nh(1),for nh2=1:Nh(2),fprintf(fh,'%f,',bcc.evaluate{Ic}.con(nh1,nh2,nroi));end; end;
            for nh1=1:Nh(1),for nh2=1:Nh(2),fprintf(fh,'%f,',bcc.evaluate{Ic}.stderr(nh1,nh2,nroi));end; end;
            %if prod(Nh)==1,fprintf(fh,'%f,',bcc.evaluate{Ic}.con(1,1,nroi)./bcc.evaluate{Ic}.t(nroi));end;
            fprintf(fh,'%f,%f,%f,%f\n',bcc.evaluate{Ic}.t(nroi),bcc.evaluate{Ic}.dof(nroi),bcc.evaluate{Ic}.p(nroi),bcc.evaluate{Ic}.pfdr(nroi));
        end
        fclose(fh);
        fprintf('Contrast file saved: %s\n',fullfile(bcc.swd,fname));
    end
end

cd(cwd);
    
end

