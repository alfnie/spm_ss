function spm_ss_results(varargin)
% SPM_SS_RESULTS Displays results of ROI- or voxel- based subject-specific analyses
%
% spm_ss_results;
% spm_ss_results(ss,Ic);
%

if nargin<1,        [ss,Ic]=spm_ss_selectcontrast; ss=spm_ss_contrast(ss,Ic,0); option='init';
elseif nargin==1,   ss=varargin{1};[ss,Ic]=spm_ss_selectcontrast(ss); ss=spm_ss_contrast(ss,Ic,0); option='init'; 
elseif nargin==2,   ss=varargin{1};Ic=varargin{2};option='init'; 
elseif nargin==3,   if isstruct(varargin{1}), data=varargin{1}; option=varargin{3}; 
                    else data=get(gcbf,'userdata'); option=varargin{3}; end; 
end

cwd=pwd;
switch(option),
    case 'init',
        data.thr.type=1;
        data.thr.value=.001;
        data.thr.ovr=.5;
        data.ss=ss;
        data.Ic=Ic;
        data.handles.fig=figure('units','norm','position',[.2,.1,.6,.8],'color','w','name',[mfilename,' : ',ss.evaluate{Ic}.name],'numbertitle','off','colormap',gray);
        if ss.typen==1,
            data.handles.ax1=axes('units','norm','position',[.2,.3,.6,.6]);
            %cdata=spm_ss_mip(ones(1,0),ones(3,0));data.handles.im1=imagesc(cdata);axis equal;axis tight;data.handles.ax1=gca;
            set(data.handles.ax1,'xtick',[],'ytick',[],'xcolor','w','ycolor','w');
            data.handles.cb1=colorbar; set(data.handles.cb1,'units','norm','position',[.8,.4,.02,.4],'ydir','reverse','ytick',[0,64*8/9],'yticklabel',[]);         
            str=sprintf('%8s%7s%7s%6s%7s%7s%18s','k','h','F','dof','p','p-fdr','(x,y,z)   ');
            uicontrol('units','norm','position',[.05,.2,.9,.035],'style','text','backgroundcolor','w','foregroundcolor','k','string',str,'horizontalalignment','left','fontname','monospaced');
            data.opts.displays={'Effect sizes (whole brain)','T/F statistics (whole brain)','Proportion subjects (whole brain)'};
            data.opts.displayn=[1:3;ones(1,3)];
            data.opts.display=1;
        else
            data.handles.ax1=axes('units','norm','position',[.2,.3,.6,.6]);
            %cdata=spm_ss_mip(ones(1,0),ones(3,0));data.handles.im1=imagesc(cdata);axis equal;axis tight;data.handles.ax1=gca;
            set(data.handles.ax1,'xtick',[],'ytick',[],'xcolor','w','ycolor','w','visible','off');
            data.handles.cb1=colorbar; set(data.handles.cb1,'units','norm','position',[.8,.4,.02,.4],'ydir','reverse','ytick',[0,64*8/9],'yticklabel',[]);
            data.handles.ax2=axes('units','norm','position',[.3,.4,.4,.4],'xtick',[],'ytick',[],'xcolor','w','ycolor','w');
            str=sprintf('%3s%8s%7s%7s%6s%7s%7s%18s','ROI','k','h','F','dof','p','p-fdr','(x,y,z)   ');
            uicontrol('units','norm','position',[.05,.2,.9,.035],'style','text','backgroundcolor','w','foregroundcolor','k','string',str,'horizontalalignment','left','fontname','monospaced');
            data.opts.displays={'Effect sizes (barplot)','T/F statistics (barplot)','Proportion subjects (barplot)','Effect sizes (whole brain)','T/F statistics (whole brain)','Proportion subjects (whole brain)'};
            data.opts.displayn=[[1:3,1:3];[2+zeros(1,3),1+zeros(1,3)]];
            data.opts.display=1;
        end
        data.handles.frame=uicontrol('units','norm','position',[0,.9,1,.1],'style','frame','backgroundcolor',0*[1,1,1],'foregroundcolor','w');
        data.handles.list=uicontrol('units','norm','position',[.05,.05,.9,.15],'style','listbox','string','');
        if ss.typen==1, uicontrol('units','norm','position',[.01,.94,.29,.035],'style','text','backgroundcolor',0*[1,1,1],'foregroundcolor','w','string','voxel-level p-value < ','horizontalalignment','right','fontsize',12);
        else uicontrol('units','norm','position',[.01,.94,.29,.035],'style','text','backgroundcolor',0*[1,1,1],'foregroundcolor','w','string','ROI-level p-value < ','horizontalalignment','right','fontsize',12); end
        data.handles.thrvalue=uicontrol('units','norm','position',[.3,.94,.1,.035],'style','edit','backgroundcolor',0*[1,1,1],'foregroundcolor','w','string',num2str(data.thr.value),'callback',{@spm_ss_results,'draw'});
        data.handles.thrtype=uicontrol('units','norm','position',[.4,.925,.1,.05],'style','popupmenu','backgroundcolor',0*[1,1,1],'foregroundcolor','w','string',{'uncorrected','FDR-corrected'},'value',data.thr.type,'callback',{@spm_ss_results,'draw'});
        uicontrol('units','norm','position',[.5,.94,.3,.035],'style','text','backgroundcolor',0*[1,1,1],'foregroundcolor','w','string','proportion subjects > ','horizontalalignment','right','fontsize',12);
        data.handles.throvr=uicontrol('units','norm','position',[.8,.94,.05,.035],'style','edit','backgroundcolor',0*[1,1,1],'foregroundcolor','w','string',num2str(data.thr.ovr),'callback',{@spm_ss_results,'draw'});
        data.handles.display=uicontrol('units','norm','position',[.3,.8,.4,.1],'style','popupmenu','backgroundcolor','w','foregroundcolor','k','string',data.opts.displays,'value',data.opts.display,'callback',{@spm_ss_results,'draw'});
        
        cd(ss.swd);
        data.matdim.mat=ss.evaluate{Ic}.CON(1).mat;
        data.matdim.dim=ss.evaluate{Ic}.CON(1).dim;
        data.CON=spm_read_vols(ss.evaluate{Ic}.CON);
        if size(data.CON,4)>1,data.CON=sqrt(sum(abs(data.CON).^2,4));end
        data.T=spm_read_vols(ss.evaluate{Ic}.T);
        data.P=spm_read_vols(ss.evaluate{Ic}.P);
        data.Pfdr=spm_read_vols(ss.evaluate{Ic}.PFDR);
        data.OVR=spm_read_vols(ss.estimate.OVERLAP);
        data.DOF=spm_read_vols(ss.evaluate{Ic}.DOF);
        if data.ss.typen>1,
            data.con=ss.evaluate{Ic}.con;
            data.con=reshape(data.con,[],size(data.con,3));
            data.stderr=ss.evaluate{Ic}.stderr;
            data.stderr=reshape(data.stderr,[],size(data.stderr,3));
            data.t=ss.evaluate{Ic}.t;
            data.p=ss.evaluate{Ic}.p;
            data.pfdr=ss.evaluate{Ic}.pfdr;
            data.ovr=ss.estimate.overlap;
            data.dof=ss.evaluate{Ic}.dof;
            [XYZ{1},XYZ{2},XYZ{3}]=ndgrid(1:ss.VN(1).dim(1),1:ss.VN(1).dim(2),1:ss.VN(1).dim(3));
            XYZ=reshape(cat(4,XYZ{:}),[],3)';XYZ=ss.VN(1).mat*cat(1,XYZ,ones(1,size(XYZ,2)));
            data.froi=reshape(round(spm_get_data(ss.VM,pinv(ss.VM(1).mat)*XYZ))',[ss.VN(1).dim(1:3),numel(ss.VM)]);
            %data.froi=round(spm_read_vols(ss.VM));
            nrois=max(data.froi(:));
            for nr=1:nrois,idx=find(data.froi==nr);data.size(nr)=numel(idx);[idx1,idx2,idx3]=ind2sub(size(data.froi),idx(:));data.xyz{nr}=mean([idx1,idx2,idx3],1);end
        end
        cd(cwd);
        set(data.handles.fig,'userdata',data);
        spm_ss_results(data,[],'draw');

    case 'list'
        nr=get(data.handles.list,'value');
        switch(data.opts.display(2)),
            case 1,%volume
                switch(data.ss.typen)
                    case 1,%voxel
                        ncluster=data.clusters.idxsort(nr);
                        npeak=data.clusters.idxsortpeak(nr);
                        idxnr=find(data.clusters.C==ncluster);
                        [xt,yt,zt]=ind2sub(size(data.X),data.IDX(idxnr));xyz=data.matdim.mat(1:3,:)*[xt,yt,zt,ones(size(xt))]';
                        [nill,cdata]=spm_ss_mip(ones(1,numel(idxnr)),xyz,data.matdim.mat);
                        xyz=data.matdim.mat(1:3,:)*[data.clusters.peakxyz{ncluster}(:,npeak);1];
                        [nill,cdata0]=spm_ss_mip(1,xyz,data.matdim.mat);
                        if isfield(data.handles,'cont')&&ishandle(data.handles.cont),delete(data.handles.cont);end
                        axes(data.handles.ax1);hold on; [nill,data.handles.cont]=contour(double(cdata>0),[.5,.5],'r'); hold off;
                        set(data.handles.cont,'linewidth',2);
                        if isfield(data.handles,'lines')&&any(ishandle(data.handles.lines(:))),delete(data.handles.lines(ishandle(data.handles.lines)));end
                        hold on; idx1=find(any(cdata0,1)); idx2=find(any(cdata0,2)); idx1(diff(idx1)==1)=[];idx2(diff(idx2)==1)=[];data.handles.lines=[plot([1;1]*idx1(:)',repmat(get(data.handles.ax1,'ylim')',[1,numel(idx1)]),'k--'); plot(repmat(get(data.handles.ax1,'xlim')',[1,numel(idx2)]),[1;1]*idx2(:)','k--')]; hold off;
                        set(data.handles.lines,'linewidth',1,'color',.75*[1,1,1]);
                    case {2,3},%roi
                        nroi=data.idx(nr);
                        idxnr=find(data.froi==nroi);
                        [xt,yt,zt]=ind2sub(size(data.X),(idxnr));xyz=data.matdim.mat(1:3,:)*[xt,yt,zt,ones(size(xt))]';
                        [nill,cdata]=spm_ss_mip(ones(1,numel(idxnr)),xyz,data.matdim.mat);
                        if isfield(data.handles,'cont')&&ishandle(data.handles.cont),delete(data.handles.cont);end
                        axes(data.handles.ax1);hold on; [nill,data.handles.cont]=contour(double(cdata>0),[.5,.5],'r'); hold off;
                        set(data.handles.cont,'linewidth',2);
                end
            case 2,%bar
        end
        set(data.handles.fig,'userdata',data);
        
    case 'draw'
        data.thr.type=get(data.handles.thrtype,'value');
        data.thr.value=str2double(get(data.handles.thrvalue,'string'));
        data.thr.ovr=str2double(get(data.handles.throvr,'string'));
        data.opts.display=data.opts.displayn(:,get(data.handles.display,'value'));
        switch(data.ss.typen),
            case 1,
                CON=data.CON;P=data.P;OVR=data.OVR;T=data.T;DOF=data.DOF;
            case {2,3},
                CON=data.con;P=data.p;OVR=data.ovr;T=data.t;dof=data.dof;
        end
        switch(data.opts.display(2)),
            case 1,
                switch(data.opts.display(1)),
                    case 1, data.X=data.CON;
                    case 2, data.X=data.T;
                    case 3, data.X=data.OVR;
                end
            case 2,
                switch(data.opts.display(1)),
                    case 1, data.X=data.con;
                    case 2, data.X=data.t; 
                    case 3, data.X=data.ovr;
                end
        end
        data.idx=find(OVR>data.thr.ovr); 
        Pfdr=P;Pfdr(data.idx)=spm_ss_fdr(Pfdr(data.idx)); % FDR correction over subset above proportion-overlap threshold
        %Pfdr=P;Pfdr(:)=spm_ss_fdr(Pfdr(:)); % FDR correction over all voxels/ROIs
        if data.thr.type==2,P=Pfdr; end
        data.idx=data.idx(find(P(data.idx)<data.thr.value|data.thr.value==1));
        if data.ss.typen==1, data.IDX=data.idx; else data.IDX=[];for nr=1:numel(data.idx),data.IDX=cat(1,data.IDX,find(data.froi==data.idx(nr)));end; end
        
        switch(data.opts.display(2)),
            case 1,
                if ~isempty(data.IDX), [idx1,idx2,idx3]=ind2sub(size(data.X),data.IDX);XYZ=data.matdim.mat(1:3,:)*[idx1,idx2,idx3,ones(size(idx1))]'; else XYZ=zeros(3,0);end
                cdata=spm_ss_mip(data.X(data.IDX)',XYZ,data.matdim.mat);
                axes(data.handles.ax1);cla;
                data.handles.im1=imagesc(cdata);axis equal;axis tight;
                data.handles.ax1=gca;set(data.handles.ax1,'xtick',[],'ytick',[],'xcolor','w','ycolor','w');
                data.handles.cb1=colorbar; set(data.handles.cb1,'units','norm','position',[.8,.4,.02,.4],'ydir','reverse','ytick',[0,64*8/9],'yticklabel',[]);
                set(data.handles.cb1,'ydir','reverse','ytick',[0,64*8/9],'yticklabel',{num2str(max(data.X(data.IDX)),'%0.3f'),num2str(min(data.X(data.IDX)),'%0.3f')});
                set([data.handles.cb1;data.handles.ax1;get(data.handles.ax1,'children')],'visible','on');
                if isfield(data.handles,'ax2'), set([data.handles.ax2;get(data.handles.ax2,'children')],'visible','off'); end
                drawnow;
            case 2,
                axes(data.handles.ax2);cla;
                x=repmat(1:numel(data.idx),[size(data.X,1)+2,1])+repmat(linspace(-1/2,1/2,size(data.X,1)+2)',[1,numel(data.idx)]);x=x(2:end-1,:);
                if ~isempty(data.idx),
                    h=[];for n1=1:size(x,1),h(n1)=bar(x(n1,:),data.X(n1,data.idx),1/(size(x,1)+2));hold on; end; hold off; %h=bar(x',data.X(:,data.idx)');
                    set(h,'facecolor',.75*[1,1,1],'edgecolor','none');
                end
                set(gca,'xtick',1:numel(data.idx),'xticklabel',[],'xlim',[0,numel(data.idx)+1]);
                if numel(data.idx)<=20, hx=1:numel(data.idx); else hx=round(linspace(1,numel(data.idx),20)); end
                if data.opts.display(1)==1,
                    stderr=data.stderr(:,data.idx);
                    hold on; errorbar(x,data.X(:,data.idx),2*data.stderr(:,data.idx),'k.'); hold off; %hold on; errorbar(1:numel(data.idx),data.X(:,data.idx),2*stderr,'k.'); hold off;
                end
                hold on; h=text(hx,get(gca,'ylim')*[1;0]+zeros(1,numel(hx)),cellstr([repmat('ROI # ',[numel(hx),1]),num2str(reshape(data.idx(hx),[numel(hx),1])),repmat('  -',[numel(hx),1])])'); hold off;
                set(h,'rotation',90,'horizontalalignment','right');
                data.handles.ax2=gca;set([data.handles.ax2;get(data.handles.ax2,'children')],'visible','on');
                if isfield(data.handles,'ax1'), set([data.handles.cb1;data.handles.ax1;get(data.handles.ax1,'children')],'visible','off'); end
                drawnow;
        end
        
        switch(data.ss.typen)
            case 1,
                if ~isempty(data.idx), [xt,yt,zt]=ind2sub(size(data.X),data.idx);xyz=[xt,yt,zt]';data.clusters.C=spm_clusters(xyz);nclusters=max(data.clusters.C); else xyz=zeros(3,0); data.clusters.C=[]; nclusters=0; end
                data.clusters.size=zeros(1,nclusters);data.clusters.peakxyz=cell(1,nclusters);data.clusters.peakidx=cell(1,nclusters);for nr=1:nclusters,idxnr=find(data.clusters.C==nr);data.clusters.size(nr)=numel(idxnr);data.clusters.peakidx{nr}=spm_get_lm(data.T,xyz(:,idxnr)); if isempty(data.clusters.peakidx{nr}), [nill,data.clusters.peakidx{nr}]=max(data.T(data.idx(idxnr))); else [nill,idxsort]=sort(-data.T(data.idx(idxnr(data.clusters.peakidx{nr})))); data.clusters.peakidx{nr}=data.clusters.peakidx{nr}(idxsort); end; data.clusters.peakxyz{nr}=xyz(:,idxnr(data.clusters.peakidx{nr})); data.clusters.peakidx{nr}=data.idx(idxnr(data.clusters.peakidx{nr})); end
                [nill,idxsort]=sort(-data.clusters.size);
                txt={};
                data.clusters.idxsort=[];data.clusters.idxsortpeak=[];
                for ncidx=1:nclusters,
                    nr=idxsort(ncidx);
                    for npeak=1:numel(data.clusters.peakidx{nr}),
                        data.clusters.idxsort(end+1)=nr;data.clusters.idxsortpeak(end+1)=npeak;
                        temp=['( ',sprintf('%+03.0f ',(data.matdim.mat(1:3,:)*[data.clusters.peakxyz{nr}(:,npeak);1])'),') '];
                        temp=[repmat(' ',[1,18-length(temp)]),temp];
                        k=data.P(data.clusters.peakidx{nr}(npeak));if k>0&&k<.001,lk=floor(log10(k));sk=[num2str(floor(k*(10^abs(lk)))),'e',num2str(lk)];else sk=sprintf('%.3f',k);sk=sk(2:end); end;temp1=sk;
                        k=Pfdr(data.clusters.peakidx{nr}(npeak));if k>0&&k<.001,lk=floor(log10(k));sk=[num2str(floor(k*(10^abs(lk)))),'e',num2str(lk)];else sk=sprintf('%.3f',k);sk=sk(2:end); end;temp2=sk;
                        if npeak==1,temp0=sprintf('%8d',data.clusters.size(nr));else temp0=repmat(' ',[1,8]); end
                        txt{end+1}=[...
                            temp0,...
                            sprintf('%7.3f',data.CON(data.clusters.peakidx{nr}(npeak))),...
                            sprintf('%7.2f',data.T(data.clusters.peakidx{nr}(npeak))),...
                            sprintf('%6.1f',data.DOF(data.clusters.peakidx{nr}(npeak))),...
                            sprintf('%7s',temp1),...
                            sprintf('%7s',temp2),...
                            temp];
                    end
                end
                set(data.handles.list,'string',txt,'fontname','monospaced','callback',{@spm_ss_results,'list'},'value',max(1,min(length(txt),get(data.handles.list,'value'))));
            case {2,3},
                txt={};
                for nr=1:numel(data.idx),
                    nroi=data.idx(nr);
                    temp=['( ',sprintf('%+03.0f ',(data.matdim.mat(1:3,:)*[data.xyz{nroi}';1])'),') '];
                    temp=[repmat(' ',[1,18-length(temp)]),temp];
                    k=data.p(nroi);if k>0&&k<.001,lk=floor(log10(k));sk=[num2str(floor(k*(10^abs(lk)))),'e',num2str(lk)];else sk=sprintf('%.3f',k);sk=sk(2:end); end;temp1=sk;
                    k=Pfdr(nroi);if k>0&&k<.001,lk=floor(log10(k));sk=[num2str(floor(k*(10^abs(lk)))),'e',num2str(lk)];else sk=sprintf('%.3f',k);sk=sk(2:end); end;temp2=sk;
                    datacon=data.con(:,nroi);if numel(datacon)>1,datacon=sqrt(sum(abs(datacon).^2));end
                    txt{end+1}=[...
                        sprintf('%-3d',nroi),...
                        sprintf('%8d',data.size(nroi)),...
                        sprintf('%7.3f',datacon),...
                        sprintf('%7.2f',data.t(nroi)),...
                        sprintf('%6.1f',data.dof(nroi)),...
                        sprintf('%7s',temp1),...
                        sprintf('%7s',temp2),...
                        temp];
                end
                set(data.handles.list,'string',txt,'fontname','monospaced','callback',{@spm_ss_results,'list'},'value',max(1,min(length(txt),get(data.handles.list,'value'))));
        end
        set(data.handles.fig,'userdata',data);
end

        
end

function [X,Y]=spm_ss_mip(Z,XYZ,M,units)
% Karl Friston et al.
% $Id: spm_mip.m 3348 2009-09-03 10:32:01Z guillaume $

%-Get units and grid scaling
%--------------------------------------------------------------------------
global defaults
try, Grid = defaults.grid; catch, Grid = 0.4;               end
try, units;                catch, units = {'mm' 'mm' 'mm'}; end
try, M;                    catch, M = 1; end

if size(M,1)   == 1, M   = speye(4,4)*M; end

%-Scale & offset point list values to fit in [0.25,1]
%==========================================================================
Z    = Z - min(Z);
mx   = max(Z);
Scal = 8;
if isempty(mx),
elseif isfinite(mx) && (numel(Z) ~= 1),
    Z = (1 + Scal*Z/max(eps,mx))/(Scal + 1);
else
    Z = ones(1,length(Z));
end

%-Display format
%==========================================================================
load('MIP.mat');

%-Single slice case
%--------------------------------------------------------------------------
if isempty(units{3}) && ~strcmp(units{2},'mm')
    %-2d case: Time-Frequency or Frequency-Frequency
    %----------------------------------------------------------------------
    mip = 4*grid_trans;
elseif isempty(units{3})
    %-2d case
    %----------------------------------------------------------------------
    mip = 4*grid_trans + mask_trans;
elseif strcmp(units{3},'ms') || strcmp(units{3},'Hz')
    %-3d case: Space-time
    %----------------------------------------------------------------------
    mip = 4*grid_time + mask_trans;
else
    %-3d case: Space
    %----------------------------------------------------------------------
    mip = 4*grid_all + mask_all;
end

% Load mip and create maximum intensity projection
%--------------------------------------------------------------------------
mip  = mip/max(mip(:));
c    = [0 0 0 ;
        0 0 1 ;
        0 1 0 ;
        0 1 1 ;
        1 0 0 ;
        1 0 1 ; 
        1 1 0 ; 
        1 1 1 ] - 0.5;
c    = c*M(1:3,1:3);
dim  = [(max(c) - min(c)) size(mip)];
d    = spm_project(Z,round(XYZ),dim);
mip  = max(d,Grid*mip);
X=(rot90((1 - mip)*64));
Y=(rot90(d));
end















