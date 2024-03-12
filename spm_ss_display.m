function spm_ss_display(activationfiles,roifiles,thresholded,roiindexes,roinames)

% SPM_SS_DISPLAY renders surface plots
% 

data.init=1;
if ~nargin
    if strcmp(mfilename,'spm_ss_display')
        str='Select spm_ss*.mat analysis file';
        disp(str);
        Pdefault='';objname=findobj('tag','spm_ss');if numel(objname)==1,objdata=get(objname,'userdata');if isfield(objdata,'files_spm_ss'),Pdefault=objdata.files_spm_ss;end;end;
        P=spm_select(1,'^SPM_ss.*\.mat$',str,{Pdefault});
    else P=[];
    end
    if ~isempty(P),
        if numel(objname)==1&&~isempty(P),objdata.files_spm_ss=P;set(objname,'userdata',objdata);end;
        load(P);
        ss.swd=fileparts(P);
        activationfiles={};
        thresholded={};
        posstr=1;
        
        if ss.typen>1,
            if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
            str='Display ROIs?';
            disp(str);
            dispact=spm_input(str,posstr,'m','display ROIs|none',[1,2], 1);posstr='+1';
            if dispact==1
                try
                    if isempty(spm_fileparts(ss.PM{1})), roifiles={fullfile(ss.swd,ss.PM{1})}; else roifiles=cellstr(ss.PM); end
                catch
                    if isempty(spm_fileparts(ss.PM)), roifiles={fullfile(ss.swd,ss.PM)}; else roifiles={ss.PM}; end
                end
            else roifiles={};
            end
        else roifiles={};
        end        
        roinames={}; 
        roiindexes=cell(1,numel(roifiles));

        if (ss.typen==1&&isfield(ss,'estimate')&&isfield(ss.estimate,'OVERLAP'))||(ss.typen>1&&isfield(ss,'PM1')),
            if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
            str='Display overlap of subject-specific masks?';
            disp(str);
            dispact=spm_input(str,posstr,'m','display overlap map|none',[1,2], 1);posstr='+1';
            if dispact==1,
                if ss.typen==1,  activationfiles=cat(2,activationfiles,{fullfile(ss.swd,ss.estimate.OVERLAP.fname)});
                else             activationfiles=cat(2,activationfiles,{fullfile(ss.swd,ss.PM1)});
                end
                thresholded=cat(2,thresholded,{nan});
            end
        end

        if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
        str='Display individual subject-specific masks?';
        disp(str);
        if ss.files_selectmanually, str2='thresholded maps (one image per subject and per partition/session) | thresholded maps (one image per subject)|none'; dispact=spm_input(str,posstr,'m',str2,[1,2,5], 2);posstr='+1'; 
        else str2='thresholded maps (one image per subject and per partition/session) | thresholded maps (one image per subject) | unthresholded maps (localizer contrast volumes, one or several per subject)| unthresholded maps and allow me to select a different threshold for these maps | none';disp(str2);dispact=spm_input(str,posstr,'m',str2,[1,2,3,4,5], 2);posstr='+1'; 
        end
        if dispact==1,
            activationfiles=cat(2,activationfiles,ss.PN); thresholded=cat(2,thresholded,repmat({.5},[1,numel(ss.PN)]));
        elseif dispact==2,
            if isfield(ss,'PL'), temp1=spm_vol(fullfile(ss.swd,ss.PL)); for nt=1:numel(temp1), activationfiles=cat(2,activationfiles,{fullfile(ss.swd,[ss.PL,',',num2str(nt)])}); thresholded=cat(2,thresholded,{.5}); end; clear temp1;
            else
                disp('Warning: This model had been run on an older version of spm_ss. Displaying thresholded masks (one image per subject and partition/session) instead.');
                activationfiles=cat(2,activationfiles,ss.PN); thresholded=cat(2,thresholded,repmat({.5},[1,numel(ss.PN)])); 
            end
        elseif (dispact==3||dispact==4)&&~ss.files_selectmanually,
            if dispact==4
                str='Enter threshold values (leave empty to display unthresholded data; you can also enter an expression (e.g. p<.05) if displaying an SPM''s T or F contrast)';
                disp(str);
                answer=inputdlg({ss.Localizer_contrasts{1,:}},str,size(ss.Localizer_contrasts,2),repmat({'pfdr<.05'},[1,size(ss.Localizer_contrasts,2)]));
            end
            disp('Importing contrast information from SPM.mat files. Please wait...');            
            for np=1:size(ss.EffectOfInterest_spm,2),%subjects
                spm_data=[];
                Ic2=[];Ec2=[];ok=1;
                for ncontrast=1:size(ss.Localizer_contrasts,2),%contrasts
                    current_spm=ss.Localizer_spm{min(size(ss.Localizer_spm,1),ncontrast),np};
                    [spm_data,SPM,Ec2(ncontrast)]=spm_ss_importspm(spm_data,current_spm);
                    Cnames={SPM.xCon(:).name};
                    for ncrossvalid=1:size(ss.Localizer_contrasts,1),%crossvalidation partitions
                        temp=strmatch(ss.Localizer_contrasts{ncrossvalid,ncontrast},Cnames,'exact');if numel(temp)~=1&&~isempty(strmatch('-not',ss.Localizer_contrasts{ncrossvalid,ncontrast})),temp=-strmatch(ss.Localizer_contrasts{ncrossvalid,ncontrast}(5:end),Cnames,'exact');end;if numel(temp)~=1,ok=0;else Ic2(ncrossvalid,ncontrast)=temp;end
                        if ~ok, error(['the target contrasts are not found inside ',current_spm]); 
                        else
                            activationfiles=cat(2,activationfiles,{fullfile(SPM.swd,SPM.xCon(abs(temp)).Vspm.fname)});
                            if dispact==4, thresholded=cat(2,thresholded,{answer{ncontrast}}); 
                            else thresholded=cat(2,thresholded,{nan}); end
                        end
                    end
                end
            end
        end
    else
        activationfiles={};
        roifiles={};
        roinames={}; 
        roiindexes=cell(1,numel(roifiles));
        thresholded=cell(1,numel(activationfiles));
    end
elseif isempty(activationfiles),
    activationfiles={};
    roifiles={};
    roinames={};
    roiindexes=cell(1,numel(roifiles));
    thresholded=cell(1,numel(activationfiles));
elseif isstruct(activationfiles),
    data=activationfiles;
    activationfiles=data.activationfiles;
    roifiles=data.roifiles;
    roinames=data.roinames;
    roiindexes=data.roiindexes;
    thresholded=data.thresholded;
elseif ischar(activationfiles),
    switch(lower(activationfiles)),
        case 'options',
            hfig=gcbf;if isempty(get(hfig,'userdata')), hfig=gcf; end
            data=get(hfig,'userdata');
            switch(lower(roifiles))
                case 'renderfile'
                    data.rendfile=thresholded; % file to render activations
                    my_spm_render;
                    data.plot_background=my_spm_render([],1,data.rendfile,'display');
                    data.plot_rois=cell(1,numel(data.plot_background));for nview=1:numel(data.plot_background),data.plot_rois{nview}=zeros([size(data.plot_background{nview},1),size(data.plot_background{nview},2),2]); end
                    data.init_activationfiles=1:numel(data.activationfiles);
                    data.init_roifiles=1;
                    data.init_display=1;
                otherwise
                    data.options.(roifiles)=thresholded;
                    set(hfig,'userdata',data);
            end
            %return
        case 'mosaic',
            hfig=gcbf;if isempty(get(hfig,'userdata')), hfig=gcf; end
            data=get(hfig,'userdata');
            n0=get(data.handles{1},'value');
            ins={[1,2,3,4,5,6],[1,2;1,3;2,3;2,2;2,1;2,4]};%btrlfb
            a=cell(2,4);
            for n1=1:numel(ins{1}),a{ins{2}(n1,1),ins{2}(n1,2)}=zeros(size(data.plot_background{n1}));end;
            a{1,1}=zeros(size(a{1,2},1),size(a{2,1},2),3);a{1,4}=zeros(size(a{1,3},1),size(a{2,4},2),3);
            cx=zeros(2,4);cy=zeros(2,4);for n1=1:numel(a),cx(n1)=size(a{n1},2);cy(n1)=size(a{n1},1);end;cx=[zeros(size(cx,1),1),cumsum(cx,2)];cy=[zeros(1,size(cy,2));cumsum(cy,1)];
            b=cell2mat(a);hfig2=figure('name',[mfilename,' mosaic view'],'numbertitle','off','units','norm','position',[0,.3,1,.6],'color','k','colormap',jet);him=imagesc(b);axis equal;axis off;hax=gca;set(hax,'ydir','normal'); figure(hfig);
            for n1=1:numel(ins{1}),set(data.handles{1},'value',ins{1}(n1));spm_ss_display('redraw');data=get(hfig,'userdata');
                a{ins{2}(n1,1),ins{2}(n1,2)}=get(data.handle_display,'cdata');b=cell2mat(a);set(him,'cdata',b); 
                if ~isempty(data.handles_contour)&&ishandle(data.handles_contour), hc(ins{2}(n1,1),ins{2}(n1,2))=data.handles_contour; set(data.handles_contour,'parent',hax,'xdata',get(data.handles_contour,'xdata')+cx(ins{2}(n1,1),ins{2}(n1,2)),'ydata',get(data.handles_contour,'ydata')+cy(ins{2}(n1,1),ins{2}(n1,2))); else hc(ins{2}(n1,1),ins{2}(n1,2))=nan; end
            end
            a{1,2}=a{1,2}(end:-1:1,end:-1:1,:); if ~isnan(hc(1,2)),set(hc(1,2),'xdata',cx(1,2)+size(a{1,2},2)-(get(hc(1,2),'xdata')-cx(1,2)),'ydata',cy(1,2)+size(a{1,2},1)-(get(hc(1,2),'ydata')-cy(1,2)));end
            b=cell2mat(a);set(him,'cdata',b);
            set(data.handles{1},'value',n0);spm_ss_display('redraw');
            return
        case 'slice'
            hfig=gcbf;if isempty(get(hfig,'userdata')), hfig=gcf; end
            data=get(hfig,'userdata');
            ma=spm_vol(fullfile(fileparts(which('spm')),'canonical','avg152T1.nii')); % image to superimpose the contours on
            mb=spm_read_vols(ma);
            direction=2; % 1:x 2:y 3:z
            plane=setdiff(1:3,direction);
            [x,y]=ndgrid(1:ma.dim(plane(1)),1:ma.dim(plane(2)));
                        
            nrois=0;for nfiles=1:numel(data.roifiles),nrois=nrois+numel(data.roiindexes{nfiles});end
            answer=inputdlg({'Number of slices?'},'Number of slices:',1,{'15'});
            nplots=str2num(answer{1});%15;
            nplots2=[floor(sqrt(nplots*3/5)),ceil(nplots/floor(sqrt(nplots*3/5)))];
            map=get(hfig,'colormap');%jet(nrois);
            if size(map,1)~=nrois, map=map(round(linspace(1,size(map,1),nrois+2)),:);map=map(2:end-1,:); end

            nfiles=data.nfiles; nfiles(isnan(nfiles))=[];
            mb3=repmat(mb,[1,1,1,3]);
            if ~isempty(data.plot_activation)&&~isempty(nfiles),
                [xx,yy,zz]=ndgrid(1:ma.dim(1),1:ma.dim(2),1:ma.dim(3));
                xyzxyz=[xx(:),yy(:),zz(:),ones(numel(zz),1)]';
                clear a1 b1;
                 for nfile=1:numel(nfiles),
                    a1(nfile)=spm_vol(data.activationfiles{nfiles(nfile)});
                    b1{nfile}=reshape(spm_get_data(a1(nfile),pinv(a1(nfile).mat)*ma.mat*xyzxyz),ma.dim(1:3));%spm_read_vols(a1(nfile));
                    if ischar(data.thresholded{nfiles(nfile)}), b1{nfile}=spm_ss_display_str2num(data.thresholded{nfiles(nfile)},b1{nfile},a1(nfile).descrip);
                    elseif ~isnan(data.thresholded{nfiles(nfile)}), b1{nfile}=double(b1{nfile}>data.thresholded{nfiles(nfile)}); end
                end                
                if numel(nfiles)==1, thisplot=zeros(size(mb3));thisplot(:,:,:,1)=max(0,b1{1}/max(b1{1}(:)));
                elseif isempty(nfiles), thisplot=mb3;
                elseif all(~isnan(nfiles))&&numel(unique(nfiles))==1,
                    nfile=unique(nfiles);
                    cmap=get(data.hfig,'colormap');
                    thisplot=reshape(max(0,min(1, ind2rgb(round(1+max(0,size(cmap,1)*b1{nfile}(:))/max(b1{nfile}(:))),[zeros(1,3);(cmap)]) )),size(mb3));
                else
%                     ntot=0;
%                     idx=1:3;if ~isnan(nfiles(1)), thisplot=max(0,b1{nfiles(1)})/max(b1{nfiles(1)}(:)); ntot=ntot+1; else thisplot=0; end
%                     for nfile=2:min(3,numel(nfiles)),
%                         idx=[idx(3),idx(1:2)];
%                         b1t=max(0,b1{nfile})/max(b1{nfile}(:));
%                         %if ~isnan(nfiles(nfile)), thisplot=max(thisplot,data.plot_activation{nfiles(nfile)}{nview}(:,:,idx)); ntot=ntot+1; end
%                         if ~isnan(nfiles(nfile)), thisplot=thisplot+b1t(:,:,idx); ntot=ntot+1; end  % mixing options
%                     end
%                     thisplot=thisplot/max(eps,ntot); %min(3,numel(nfiles));
                    thisplot=zeros(size(mb3)); 
                    for nfile=1:min(3,numel(nfiles)),
                        if ~isnan(nfiles(nfile)), thisplot(:,:,:,nfile)=max(0,b1{nfile})/max(b1{nfile}(:)); end
                    end
                end
            else thisplot=mb3; end                
            facti=1;idx=find(all(~thisplot,4));thisplot(idx)=mb(idx)/facti;thisplot(idx+numel(mb))=mb(idx)/facti;thisplot(idx+2*numel(mb))=mb(idx)/facti;
            
            hfig2=figure('name',[mfilename,' slice view'],'numbertitle','off','units','norm','position',[0,.3,1,.6],'color','k');
            h=[];hi=[];
            for n2=1:nplots,
                h(n2)=subplot(nplots2(1),nplots2(2),n2);
                slice=ceil(n2*ma.dim(direction)/(nplots+1));
                if direction==1,        mbt=permute(thisplot(slice,:,:,:),[2,3,4,1]); 
                                        hi(n2)=imagesc(max(0,min(1, permute(mbt(end:-1:1,end:-1:1,:),[2,1,3])) ));
                                        axis equal; axis off;
                                        hold on; ht=text(size(mbt,1),size(mbt,2)*.05,['x=',num2str(round(ma.mat(1,:)*[slice,0,0,1]'),'%0.0f')]); hold off; set(ht,'color','w','fontsize',9,'horizontalalignment','right');
                elseif direction==2,    mbt=permute(thisplot(:,slice,:,:),[1,3,4,2]); 
                                        hi(n2)=imagesc(max(0,min(1, permute(mbt(end:-1:1,end:-1:1,:),[2,1,3])) ));
                                        axis equal; axis off;
                                        hold on; ht=text(size(mbt,1),size(mbt,2)*.05,['y=',num2str(round(ma.mat(2,:)*[0,slice,0,1]'),'%0.0f')]); hold off; set(ht,'color','w','fontsize',9,'horizontalalignment','right');
                else                    mbt=permute(thisplot(:,:,slice,:),[1,2,4,3]); %reshape(spm_get_data(ma,pinv(ma.mat)*ma.mat*xyz'),ma.dim(1:2));
                                        hi(n2)=imagesc(max(0,min(1, permute(mbt(end:-1:1,end:-1:1,:),[2,1,3])) ));
                                        axis equal; axis off;
                                        hold on; ht=text(size(mbt,1),size(mbt,2)*.05,['z=',num2str(round(ma.mat(3,:)*[0,0,slice,1]'),'%0.0f')]); hold off; set(ht,'color','w','fontsize',9,'horizontalalignment','right');
                end
            end
            
            for nfiles=1:numel(data.roifiles),
                a1=spm_vol(data.roifiles{nfiles});
                for n2=1:nplots,
                    slice=ceil(n2*ma.dim(direction)/(nplots+1));
                    xyz=ones(numel(x),4);
                    xyz(:,plane)=[x(:),y(:)];
                    xyz(:,direction)=slice;
                    %xyz=[x(:),y(:),slice+zeros(numel(x),1),ones(numel(x),1)];
                    b=reshape(round(spm_get_data(a1,pinv(a1.mat)*ma.mat*xyz')),ma.dim(plane));
                    axes(h(n2));
                    call=0;
                    for n1=1:numel(data.roiindexes{nfiles}),
                        c=(b==data.roiindexes{nfiles}(n1));
                        call=call|c;
                        if any(c(:)>0), 
                            hold on; [nill,hf]=contourf(flipud(fliplr(c')),[.5,.5],'k','linewidth',1); hold off;
                            if data.coloredrois, set(get(hf,'children'),'facecolor',map(n1,:),'facealpha',.25); else set(get(hf,'children'),'facecolor','none'); end
                            if data.displaylabels, hold on; ht=text(mean(size(b,1)-x(c)+1),mean(size(b,2)-y(c)+1),10,data.roinames{n1});set(ht,'color','k','horizontalalignment','center','fontweight','bold','interpreter','none'); hold off; end
                        end
                    end
                    if data.intersectactivations,
                        if direction==1,    mbt=permute(thisplot(slice,:,:,:),[2,3,4,1]); 
                                            mb3t=permute(mb3(slice,:,:,:),[2,3,4,1]); 
                        elseif direction==2,mbt=permute(thisplot(:,slice,:,:),[1,3,4,2]); 
                                            mb3t=permute(mb3(:,slice,:,:),[1,3,4,2]); 
                        elseif direction==3,mbt=permute(thisplot(:,:,slice,:),[1,2,4,3]); 
                                            mb3t=permute(mb3(:,:,slice,:),[1,2,4,3]); 
                        end
                        call=repmat(call,[1,1,3]);mbt(~call)=mb3t(~call);
                        set(hi(n2),'cdata',max(0,min(1, permute(mbt(end:-1:1,end:-1:1,:),[2,1,3])) ));
                    end
                end
            end
            for n2=1:nplots, set(h(n2),'units','norm','position',[1/nplots2(2)*rem(n2-1,nplots2(2)),1-1/nplots2(1)*ceil(n2/nplots2(2)),1/nplots2(2),1/nplots2(1)]);end
            colormap gray;
            set(hfig2,'userdata',h,'resizefcn','h=get(gcbf,''userdata'');nplots=numel(h);x=get(gcbf,''position'');nplots2=[floor(sqrt(nplots*x(4)/x(3))),ceil(nplots/floor(sqrt(nplots*x(4)/x(3))))];for n2=1:nplots, set(h(n2),''units'',''norm'',''position'',[1/nplots2(2)*rem(n2-1,nplots2(2)),1-1/nplots2(1)*ceil(n2/nplots2(2)),1/nplots2(2),1/nplots2(1)]);axis(h(n2),''equal''); end');
            figure(hfig);
            return;
            
        case 'redraw',
            data=get(gcbf,'userdata'); if isempty(data), data=get(gcf,'userdata'); end
            data.nview=get(data.handles{1},'value');
            if data.nview==7,set(data.handles{1},'value',1);spm_ss_display('mosaic');return;
            elseif data.nview==8,spm_ss_display('slice');return;end
        case 'selectroifile',
            data=get(gcbf,'userdata');
            P=spm_select(inf,'image','Select ROI file(s)',data.roifiles);
            if ~isempty(P), data.roifiles=cellstr(P); else data.roifiles={}; end
            data.init_roifiles=1;
            data.roiindexes=cell(1,numel(data.roifiles));
            data.plot_rois_idx={};
            data.plot_rois=cell(1,numel(data.plot_background));for nview=1:numel(data.plot_background),data.plot_rois{nview}=zeros([size(data.plot_background{nview},1),size(data.plot_background{nview},2),2]); end
            data.roinames={};
        case 'selectactivationfile',
            data=get(gcbf,'userdata');
            str='Select Activation file(s) (note: multiple activation files can be browsed separately or assigned to individual R/G/B color channels)'; disp(str);
            P=spm_select(inf,'image','Select Activation file(s)',data.activationfiles);
            if ~isempty(P), data.activationfiles=cellstr(P); else data.activationfiles={}; end
            data.init_activationfiles=1:numel(data.activationfiles);
            data.thresholded=cell(1,numel(data.activationfiles));            
            data.nfiles=1:min(1,numel(data.activationfiles));
        case 'coloredrois',
            data=get(gcbf,'userdata');
            data.coloredrois=get(data.handles{6},'value');
        case 'intersectactivations',
            data=get(gcbf,'userdata');
            data.intersectactivations=get(data.handles{7},'value');
        case 'displaylabels',
            data=get(gcbf,'userdata');
            data.displaylabels=get(data.handles{14},'value');
        case 'prevactivationfile',
            data=get(gcbf,'userdata');
            data.nfiles=1+mod(data.nfiles-2,numel(data.activationfiles));
            data.init_display=0;
        case 'nextactivationfile',
            data=get(gcbf,'userdata');
            data.nfiles=1+mod(data.nfiles,numel(data.activationfiles));
            data.init_display=0;
        case 'mixactivationfile',
            data=get(gcbf,'userdata');
            idx1=listdlg('name',[mfilename, ' Activation mixer'],'promptstring','RED channel:','selectionmode','single','liststring',data.activationfiles,'initialvalue',1,'listsize',[400,200]);
            idx2=listdlg('name',[mfilename, ' Activation mixer'],'promptstring','GREEN channel:','selectionmode','single','liststring',data.activationfiles,'initialvalue',1,'listsize',[400,200]);
            idx3=listdlg('name',[mfilename, ' Activation mixer'],'promptstring','BLUE channel:','selectionmode','single','liststring',data.activationfiles,'initialvalue',1,'listsize',[400,200]);
            if isempty(idx1),idx1=nan; end;if isempty(idx2),idx2=nan; end;if isempty(idx3),idx3=nan; end
            data.nfiles=[idx1,idx2,idx3];
            data.init_display=0;
        otherwise,
            disp([activationfiles, 'not recognized']);
            return;
    end
else
    if nargin<5, roinames={}; end
    if nargin<4, roiindexes=cell(1,numel(roifiles)); end
    if nargin<3, thresholded=cell(1,numel(activationfiles)); end    
end

if data.init,
    if isfield(data,'hfig')&&ishandle(data.hfig), figure(data.hfig);
    else data.hfig=figure('units','norm','position',[.3,.2,.4,.7],'name',mfilename,'numbertitle','off','color','k');
    end
    if ~iscell(activationfiles), activationfiles=cellstr(activationfiles); end
    if ~iscell(roifiles), roifiles=cellstr(roifiles); end
    if ~iscell(roiindexes), roiindexes={roiindexes}; end
    data.rendfile=fullfile(fileparts(which('spm')),'rend','render_smooth_average.mat'); % file to render activations
    %data.rendfile=fullfile(fileparts(which('spm')),'rend','render_single_subj.mat'); % file to render activations
    data.plot_background=my_spm_render([],1,data.rendfile,'display');
    data.plot_activation={};
    data.plot_rois_idx={};
    data.plot_rois=cell(1,numel(data.plot_background));for nview=1:numel(data.plot_background),data.plot_rois{nview}=zeros([size(data.plot_background{nview},1),size(data.plot_background{nview},2),2]); end
    data.activationfiles=activationfiles;
    data.roifiles=roifiles;
    data.roinames=roinames;
    data.roiindexes=roiindexes;
    data.thresholded=thresholded;
    data.nfiles=1;
    data.nview=4;
    data.displaylabels=1;
    data.coloredrois=1;
    data.intersectactivations=0;
    data.fact=.02;%02;
    data.nroi_tot=0;
    data.init_activationfiles=1:numel(data.activationfiles);
    data.init_roifiles=1;
    data.init_display=1;
    data.init=0;
end

if ~isempty(data.init_activationfiles),
    hmsg=msgbox('Rendering activation file(s). Please wait...');delete(findobj(hmsg,'type','uicontrol'));
    anyask=[];anyasknames={};anyaskvalues={};anyaskdof={};
    for nfiles=data.init_activationfiles,
        [p,f,e]=spm_fileparts(data.activationfiles{nfiles});
        switch(lower(e))
            case {'.img','.nii'}
                if isempty(data.thresholded{nfiles}),
                    a1=spm_vol(data.activationfiles{nfiles});
                    b1=spm_read_vols(a1);
                    ub1=unique(b1(:));
                    if numel(ub1)==2, data.thresholded{nfiles}=mean(ub1);
                    else
                        anyask=[anyask,nfiles];
                        defvalue=sort(b1(~isnan(b1)));if ~isempty(defvalue), defvalue=defvalue(ceil(.99*numel(defvalue))); end;
                        anyaskvalues{end+1}=num2str(defvalue);
                        anyasknames{end+1}=data.activationfiles{nfiles};
                    end
                end
        end
    end
    if ~isempty(anyask),
        str='Enter threshold values (leave empty to display unthresholded data; you can also enter an expression (e.g. p<.05) if displaying an SPM''s T or F contrast)';
        disp(str);
        answer=inputdlg(anyasknames,str,1,anyaskvalues);
        for n1=1:numel(anyask),
            temp=str2num(answer{n1});
            if ~isempty(temp)||isempty(answer{n1}), data.thresholded{anyask(n1)}=temp;
            else data.thresholded{anyask(n1)}=answer{n1}; end
        end
    end
    for nfiles=data.init_activationfiles,
        set(hmsg,'name',[mfilename,' (',num2str(nfiles),'/',num2str(numel(data.init_activationfiles)),')']);
        [p,f,e]=spm_fileparts(data.activationfiles{nfiles});
        switch(lower(e))
            case {'.img','.nii'}
                a1=spm_vol(data.activationfiles{nfiles});
                b1=spm_read_vols(a1);
                if ischar(data.thresholded{nfiles}), [b1,data.thresholded{nfiles}]=spm_ss_display_str2num(data.thresholded{nfiles},b1,a1.descrip);
                elseif ~isnan(data.thresholded{nfiles}), b1=double(b1>data.thresholded{nfiles}); end
                idx=find(b1>0);
                if isempty(idx)
                    data.plot_activation{nfiles}=data.plot_background;
                    data.range_activation{nfiles}=nan;
                else
                    dat=struct('VOL',b1,'mat',a1.mat,'dim',a1.dim);
                    data.range_activation{nfiles}=max(b1(:));
                    if ~all(isnan(data.thresholded{nfiles})),
                        data.plot_activation{nfiles}=my_spm_render(dat,data.fact,data.rendfile,'locate');
                        for nview=1:numel(data.plot_activation{nfiles}),
                            thisplot=data.plot_activation{nfiles}{nview};
                            idx=find(thisplot>0);
                            %idx=find((1-thisplot(:,:,1))<.05 & thisplot(:,:,2)<.05 & thisplot(:,:,3)<.05);
                            k=zeros([size(thisplot,1),size(thisplot,2)]);
                            k(idx)=1;
                            %k=convn(convn(k,hanning(11)/sum(hanning(11)),'same'),hanning(11)'/sum(hanning(11)),'same');
                            cmap=[1,0,0]; for n1=1:3,thisplot(:,:,n1)=(1-k).*data.plot_background{nview}(:,:,n1)+k*cmap(n1); end
                            data.plot_activation{nfiles}{nview}=thisplot;
                        end
                    else
                        [data.plot_activation{nfiles},rangergb]=my_spm_render(dat,1,data.rendfile,'display');
                        data.range_activation{nfiles}=max(rangergb);
                    end
                end
            case '.mat'
                load(data.activationfiles{nfiles},'xSPM');
                dat=struct('XYZ',xSPM.XYZ,'t',xSPM.Z','mat',xSPM.M,'dim',xSPM.DIM);
                [data.plot_activation{nfiles},rangergb]=my_spm_render(dat,1,data.rendfile,'display');
                %data.range_activation{nfiles}=max(xSPM.Z(:));
                data.range_activation{nfiles}=max(rangergb);
        end
    end
    data.init_activationfiles=[];
    close(hmsg);
end

if data.init_roifiles,
    hmsg=msgbox('Rendering ROI file(s). Please wait...');delete(findobj(hmsg,'type','uicontrol'));
    data.nroi_tot=0;
    data.plot_rois=cell(1,numel(data.plot_background));for nview=1:numel(data.plot_background),data.plot_rois{nview}=zeros([size(data.plot_background{nview},1),size(data.plot_background{nview},2),2]); end
    for nfiles=1:numel(data.roifiles),
        set(hmsg,'name',[mfilename,' (',num2str(nfiles),'/',num2str(numel(data.roifiles)),')']);
        [p,f,e]=spm_fileparts(data.roifiles{nfiles});
        a1=spm_vol(data.roifiles{nfiles});
        b1=round(spm_read_vols(a1));
        if ~isempty(dir(fullfile(p,[f,'.txt']))),
            roinames=textread(fullfile(p,[f,'.txt']),'%s','delimiter','\n');
        else roinames={}; end
        if isempty(data.roiindexes{nfiles}),
            maxb1=max(b1(:));
            for n1=numel(roinames)+1:maxb1, roinames{n1}=num2str(n1); end
            if maxb1==1, data.roiindexes{nfiles}=1; 
            else
                str=['Select ROI(s) to display from file ',data.roifiles{nfiles}];
                disp(str);
                if isempty(roinames), roinames=cellstr([repmat('',[maxb1,1]),num2str((1:maxb1)')]); end
                idx=listdlg('promptstring',str,'selectionmode','multiple','liststring',roinames,'initialvalue',1);
                if isempty(idx)
                    b1=double(b1>0);
                    maxb1=1;
                    roinames={''};
                    idx=1;
                end
            	data.roiindexes{nfiles}=idx; 
            end
        elseif any(isinf(data.roiindexes{nfiles})),
            data.roiindexes{nfiles}=1:maxb1;
        end
        rois=data.roiindexes{nfiles};
        dat=struct('VOL',b1,'mat',a1.mat,'dim',a1.dim);
        plot_rois=my_spm_render(dat,data.fact,data.rendfile,'locate');
        for nroi=1:numel(rois),
            for nview=1:numel(plot_rois),
                data.plot_rois_idx{data.nroi_tot+nroi}{nview}=find(plot_rois{nview}==rois(nroi));
                data.plot_rois{nview}(data.plot_rois_idx{data.nroi_tot+nroi}{nview})=data.nroi_tot+nroi;
            end
            if numel(data.roinames)<data.nroi_tot+nroi, 
                if numel(data.roifiles)==1, data.roinames{data.nroi_tot+nroi}=roinames{rois(nroi)}; 
                elseif numel(rois)==1, data.roinames{data.nroi_tot+nroi}=f;%[num2str(nfiles)]; 
                else data.roinames{data.nroi_tot+nroi}=[f,'.',roinames{rois(nroi)}]; end; 
            end
        end
        for nview=1:numel(plot_rois),
            data.plot_rois{nview}(:,:,2)=(data.plot_rois{nview}(:,:,1)>0);
        end
        data.nroi_tot=data.nroi_tot+numel(rois);
    end
    
    data.init_roifiles=0;
    close(hmsg);
end


nfiles=data.nfiles;
nview=data.nview;
coloredrois=data.coloredrois;
displaylabels=data.displaylabels;

set(data.hfig,'pointer','watch');
if ~isempty(data.plot_activation)&&~isempty(nfiles), 
    if numel(nfiles)==1, thisplot=data.plot_activation{nfiles}{nview};
    elseif isempty(nfiles), thisplot=data.plot_background{nview}; 
    elseif all(~isnan(nfiles))&&numel(unique(nfiles))==1,
        nfile=unique(nfiles);
        emph=1;cmap=get(data.hfig,'colormap');
        thisplot=max(0,min(1, ind2rgb(round(1+((data.plot_activation{nfile}{nview}(:,:,1)>.05).*(1+emph*(size(cmap,1)-1)*data.plot_activation{nfile}{nview}(:,:,3)))),[zeros(1,3);flipud(cmap)]) ));
    else
        ntot=0;
        idx=1:3;if ~isnan(nfiles(1)), thisplot=data.plot_activation{nfiles(1)}{nview}; ntot=ntot+1; else thisplot=0; end
        for nfile=2:min(3,numel(nfiles)),
            idx=[idx(3),idx(1:2)];
            %if ~isnan(nfiles(nfile)), thisplot=max(thisplot,data.plot_activation{nfiles(nfile)}{nview}(:,:,idx)); ntot=ntot+1; end
            if ~isnan(nfiles(nfile)), thisplot=thisplot+data.plot_activation{nfiles(nfile)}{nview}(:,:,idx); ntot=ntot+1; end  % mixing options 
        end
        thisplot=thisplot/max(eps,ntot); %min(3,numel(nfiles));
    end
else thisplot=data.plot_background{nview}; end
if coloredrois,
%     if numel(data.plot_rois_idx)==1, cmap=gray(numel(data.plot_rois_idx)+4);cmap=cmap(4:end-1,:); 
%     else cmap=winter(numel(data.plot_rois_idx)+2);cmap=cmap(2:end-1,:);end
    cmap=get(data.hfig,'colormap');
    if size(cmap,1)~=numel(data.plot_rois_idx), cmap=cmap(round(linspace(1,size(cmap,1),numel(data.plot_rois_idx)+2)),:);cmap=cmap(2:end-1,:); end
    for nroi=1:numel(data.plot_rois_idx),
        idx=data.plot_rois_idx{nroi}{nview}; %find(data.plot_rois{nview}(:,:,1)==nroi);
        k=data.plot_rois{nview}(:,:,2); 
        for n1=1:3,thisplot(idx+(n1-1)*size(thisplot,1)*size(thisplot,2))=(1-k(idx)).*thisplot(idx+(n1-1)*size(thisplot,1)*size(thisplot,2)) + k(idx).*thisplot(idx+(n1-1)*size(thisplot,1)*size(thisplot,2)).*cmap(nroi,n1);end;
    end
end

figure(data.hfig);
temp=data.plot_rois{nview}(:,:,1).*(data.plot_rois{nview}(:,:,2)>.5);
%if any(temp(:)), temp=temp.*(sum(thisplot,3)>.1); end
if data.intersectactivations,
    idxnull=find(~temp);
    for n=0:2,thisplot(idxnull+n*size(thisplot,1)*size(thisplot,2))=data.plot_background{nview}(idxnull+n*size(thisplot,1)*size(thisplot,2));end
end
if data.init_display,
    clf;
    hi=imagesc(thisplot);set(gca,'ydir','normal');axis equal;axis tight;axis off;
    if ~isfield(data,'options')||~isfield(data.options,'linewidth'), data.options.linewidth=1; end
    if ~isfield(data,'options')||~isfield(data.options,'linecolor'), data.options.linecolor=.5*[1,1,1]; end
    if ~isfield(data,'options')||~isfield(data.options,'fontsize'), data.options.fontsize=9; end
    if ~isfield(data,'options')||~isfield(data.options,'fontcolor'), data.options.fontcolor='w'; end
    data.handle_display=hi;
    if any(temp(:)>0),
        hold on;
        [nill,hc]=contour(temp,.5:numel(data.plot_rois_idx)+.5,'k-');set(hc,'linecolor',data.options.linecolor,'linewidth',data.options.linewidth,'linestyle','-');
        data.handles_contour=hc;
        hold off;
    else data.handles_contour=[]; end
    if displaylabels,
        for nroi=1:numel(data.plot_rois_idx),
            idx=data.plot_rois_idx{nroi}{nview}; %find(temp==nroi);
            if ~isempty(idx),
                [idxx,idxy]=ind2sub(size(temp),idx);
                ht=text(mean(idxy),mean(idxx),data.roinames{nroi});
                set(ht,'horizontalalignment','center','fontweight','bold','color',data.options.fontcolor,'fontsize',data.options.fontsize);
            end
        end
    end
    set(gca,'units','norm','position',[0,.3,1,.65]);
    data.handles_axis=gca;
    cb1=1;cb2=.96;
    uicontrol('style','frame','units','norm','position',[0,0,1,.25],'backgroundcolor',cb1*[1,1,1]);
    data.handles{1}=uicontrol('style','popupmenu','units','norm','position',[0,.20,1,.05],'string',strvcat('View: Bottom','View: Top','View: Right','View: Left','View: Front','View: Back','View: Mosaic','View: Slices'),'value',data.nview,'callback',[mfilename,'(''redraw'');']);
    data.handles{2}=uicontrol('style','pushbutton','units','norm','position',[.01,.05,.48,.05],'string','Select ROI file(s)','callback',[mfilename,'(''selectroifile'');']);
    data.handles{3}=uicontrol('style','pushbutton','units','norm','position',[.51,.05,.48,.05],'string','Select Activation file(s)','callback',[mfilename,'(''selectactivationfile'');']);
    data.handles{4}=uicontrol('style','text','units','norm','position',[.01,.02,.48,.03],'string',[num2str(numel(data.roifiles)),' file(s) selected'],'foregroundcolor','r','backgroundcolor',cb1*[1,1,1]);if numel(data.roifiles)>0, set(data.handles{4},'foregroundcolor','g'); end
    data.handles{5}=uicontrol('style','text','units','norm','position',[.51,.02,.48,.03],'string',[num2str(numel(data.activationfiles)),' file(s) selected'],'foregroundcolor','r','backgroundcolor',cb1*[1,1,1]);if numel(data.activationfiles)>0, set(data.handles{5},'foregroundcolor','g'); end
    uicontrol('style','frame','units','norm','position',[.01,.105,.98,.0995],'backgroundcolor',cb2*[1,1,1]);
    data.handles{6}=uicontrol('style','checkbox','units','norm','position',[.25,.17,.5,.03],'string','filled ROI contours','value',data.coloredrois,'foregroundcolor','k','backgroundcolor',cb2*[1,1,1],'callback',[mfilename,'(''coloredrois'');']);
    data.handles{7}=uicontrol('style','checkbox','units','norm','position',[.25,.14,.5,.03],'string','show activation within ROIs only','value',data.intersectactivations,'foregroundcolor','k','backgroundcolor',cb2*[1,1,1],'callback',[mfilename,'(''intersectactivations'');']);
    data.handles{8}=uicontrol('style','pushbutton','units','norm','position',[.39,.26,.10,.04],'string','Prev','callback',[mfilename,'(''prevactivationfile'');']);
    data.handles{9}=uicontrol('style','pushbutton','units','norm','position',[.51,.26,.10,.04],'string','Next','callback',[mfilename,'(''nextactivationfile'');']);
    data.handles{10}=uicontrol('style','pushbutton','units','norm','position',[.89,.26,.10,.04],'string','Mixer','callback',[mfilename,'(''mixactivationfile'');']);
    data.handles{11}=uicontrol('style','text','units','norm','position',[.01,.96,.98,.03],'string','','backgroundcolor','k','foregroundcolor','r','fontsize',8,'horizontalalignment','left');
    data.handles{12}=uicontrol('style','text','units','norm','position',[.01,.93,.98,.03],'string','','backgroundcolor','k','foregroundcolor','g','fontsize',8,'horizontalalignment','left');
    data.handles{13}=uicontrol('style','text','units','norm','position',[.01,.90,.98,.03],'string','','backgroundcolor','k','foregroundcolor','b','fontsize',8,'horizontalalignment','left');
    data.handles{14}=uicontrol('style','checkbox','units','norm','position',[.25,.11,.5,.03],'string','display labels','value',data.displaylabels,'foregroundcolor','k','backgroundcolor',cb2*[1,1,1],'callback',[mfilename,'(''displaylabels'');']);
    
    if isempty(data.nfiles)||any(data.nfiles>numel(data.activationfiles))||all(isnan(data.nfiles)), set([data.handles{11},data.handles{12},data.handles{13}],'string',''); 
    else for n1=1:3, if numel(data.nfiles)<n1||isnan(data.nfiles(n1)), set(data.handles{10+n1},'string',''); else set(data.handles{10+n1},'string',data.activationfiles{data.nfiles(n1)}); end; end; 
    end
    mfiles=nfiles;idxmfiles=find(~isnan(mfiles));mfiles=unique(mfiles(idxmfiles));
    %colorbar
    if numel(mfiles)==1&&numel(data.thresholded)>=mfiles&&isempty(data.thresholded{mfiles})&&~isnan(data.range_activation{mfiles}),
        if all(~isnan(data.nfiles))&&numel(data.nfiles)>1&&numel(unique(data.nfiles))==1, cmap=permute((get(data.hfig,'colormap')),[1,3,2]); 
        else cmap0=[linspace(.85,1,64)',linspace(1,0,64)'];cmap=zeros([size(cmap0,1),1,3]);cmap(:,:,idxmfiles)=cmap0(:,1+zeros(1,numel(idxmfiles)));cmap(:,:,setdiff(1:3,idxmfiles))=cmap0(:,2+zeros(1,3-numel(idxmfiles)));end
        axes('units','norm','position',[.925,.5,.025,.3]);image(cmap);set(gca,'ydir','normal'); axis off; hold on; ht=text([1,1],[1-.1*size(cmap,1),size(cmap,1)*1.1],{['0'],[num2str(data.range_activation{mfiles},'%0.2f')]});set(ht,'color',.85*[1,1,1],'horizontalalignment','center');
        data.handles{15}=gca;
    elseif numel(mfiles)==1&&numel(data.thresholded)>=mfiles&&~isempty(data.thresholded{mfiles})&&~all(isnan(data.thresholded{mfiles}))&&~isnan(data.range_activation{mfiles}),
        cmap0=[[.85*ones(32,1);ones(32,1)],[.85*ones(32,1);zeros(32,1)]];cmap=zeros([size(cmap0,1),1,3]);cmap(:,:,idxmfiles)=cmap0(:,1+zeros(1,numel(idxmfiles)));cmap(:,:,setdiff(1:3,idxmfiles))=cmap0(:,2+zeros(1,3-numel(idxmfiles)));axes('units','norm','position',[.925,.5,.025,.3]);image(cmap);set(gca,'ydir','normal'); axis off; hold on; 
        if ischar(data.thresholded{mfiles}), ht=text(1,[size(cmap,1)*1.1],{data.thresholded{mfiles}}); 
        else ht=text([1,1],[1-.1*size(cmap,1),size(cmap,1)*1.1],{['<',num2str(data.thresholded{mfiles},'%0.2f')],['>',num2str(data.thresholded{mfiles},'%0.2f')]}); end
        set(ht,'color',.85*[1,1,1],'horizontalalignment','center');
        data.handles{15}=gca;
    else data.handles{15}=[]; end
    set(data.hfig,'userdata',data);
else
    set(data.handle_display,'cdata',thisplot); 
    data.init_display=1;
    if isempty(data.nfiles)||any(data.nfiles>numel(data.activationfiles))||all(isnan(data.nfiles)), set([data.handles{11},data.handles{12},data.handles{13}],'string',''); 
    else for n1=1:3, if numel(data.nfiles)<n1||isnan(data.nfiles(n1)), set(data.handles{10+n1},'string',''); else set(data.handles{10+n1},'string',data.activationfiles{data.nfiles(n1)}); end; end; 
    end
    if ~isempty(data.handles{15})&&all(ishandle(data.handles{15})),delete(data.handles{15}); end
    mfiles=nfiles;idxmfiles=find(~isnan(mfiles));mfiles=unique(mfiles(idxmfiles));
    if numel(mfiles)==1&&numel(data.thresholded)>=mfiles&&isempty(data.thresholded{mfiles})&&~isnan(data.range_activation{mfiles}),
        if all(~isnan(data.nfiles))&&numel(data.nfiles)>1&&numel(unique(data.nfiles))==1,cmap=permute(get(data.hfig,'colormap'),[1,3,2]); 
        else cmap0=[linspace(.85,1,64)',linspace(1,0,64)']; cmap=zeros([size(cmap0,1),1,3]);cmap(:,:,idxmfiles)=cmap0(:,1+zeros(1,numel(idxmfiles)));cmap(:,:,setdiff(1:3,idxmfiles))=cmap0(:,2+zeros(1,3-numel(idxmfiles))); end
        axes('units','norm','position',[.925,.5,.025,.3]);image(cmap);set(gca,'ydir','normal'); axis off; hold on; ht=text([1,1],[1-.1*size(cmap,1),size(cmap,1)*1.1],{['0'],[num2str(data.range_activation{mfiles},'%0.2f')]});set(ht,'color',.85*[1,1,1],'horizontalalignment','center');
        data.handles{15}=gca;
    elseif numel(mfiles)==1&&numel(data.thresholded)>=mfiles&&~isempty(data.thresholded{mfiles})&&~all(isnan(data.thresholded{mfiles}))&&~isnan(data.range_activation{mfiles}),
        cmap0=[[.85*ones(32,1);ones(32,1)],[.85*ones(32,1);zeros(32,1)]];cmap=zeros([size(cmap0,1),1,3]);cmap(:,:,idxmfiles)=cmap0(:,1+zeros(1,numel(idxmfiles)));cmap(:,:,setdiff(1:3,idxmfiles))=cmap0(:,2+zeros(1,3-numel(idxmfiles)));axes('units','norm','position',[.925,.5,.025,.3]);image(cmap);set(gca,'ydir','normal'); axis off; hold on; 
        if ischar(data.thresholded{mfiles}), ht=text(1,[size(cmap,1)*1.1],{data.thresholded{mfiles}}); 
        else ht=text([1,1],[1-.1*size(cmap,1),size(cmap,1)*1.1],{['<',num2str(data.thresholded{mfiles},'%0.2f')],['>',num2str(data.thresholded{mfiles},'%0.2f')]}); end
        set(ht,'color',.85*[1,1,1],'horizontalalignment','center');
        data.handles{15}=gca;
    else data.handles{15}=[]; end
    set(data.hfig,'userdata',data);
    drawnow; 
end
set(data.hfig,'pointer','arrow');










function [b1,strthr]=spm_ss_display_str2num(strthr,b1,strdescrip)
p=[];T=[];F=[];X=b1;x=X;
try
    if strcmp(strdescrip(1:4),'SPM{')&&strcmp(strdescrip(6:7),'_['),
        temp=regexp(strdescrip,'\[([0-9\.,]*)\]','tokens');
        dof=str2num(temp{1}{1});
        if strdescrip(5)=='T', T=b1; p=1-spm_Tcdf(T,dof);
        elseif strdescrip(5)=='F', F=b1; p=1-spm_Fcdf(F,dof(1),dof(2)); end
        pfdr=p;pfdr(:)=spm_ss_fdr(p(:));
    end
    temp=eval(strthr);
    if numel(temp)==1, b1=double(b1>temp);
    elseif numel(temp)==numel(b1), b1=double(temp);
    else
        disp(['unable to implement threshold expression ',strthr,' from file ',strdescrip,'.  Displaying unthresholded data instead']);
        strthr=nan;
    end
catch
    try
        temp=eval(strthr);
        if numel(temp)==1, b1=double(b1>temp);
        elseif numel(temp)==numel(b1), b1=double(temp);
        else
            disp(['unable to implement threshold expression ',strthr,' from file ',strdescrip,'.  Displaying unthresholded data instead']);
            strthr=nan;
        end
    catch
        disp(['unable to implement threshold expression ',strthr,' from file ',strdescrip,'. Displaying unthresholded data instead']);
        strthr=nan;
    end
end



function [RGB,RANGERGB]=my_spm_render(dat,brt,rendfile,rendertype)
RGB={};RANGERGB={};
% Render blobs on surface of a 'standard' brain
% FORMAT spm_render(dat,brt,rendfile)
%
% dat      - a struct array of length 1 to 3
%            each element is a structure containing:
%            - XYZ - the x, y & z coordinates of the transformed SPM{.}
%                    values in units of voxels.
%            - t   - the SPM{.} values.
%            - mat - affine matrix mapping from XYZ voxels to MNI.
%            - dim - dimensions of volume from which XYZ is drawn.
% brt      - brightness control:
%            If NaN, then displays using the old style with hot
%            metal for the blobs, and grey for the brain.
%            Otherwise, it is used as a ``gamma correction'' to
%            optionally brighten the blobs up a little.
% rendfile - the file containing the images to render on to (see also
%            spm_surf.m) or a surface mesh file.
%
% Without arguments, spm_render acts as its own UI.
%__________________________________________________________________________
% 
% spm_render prompts for details of up to three SPM{.}s that are then
% displayed superimposed on the surface of a 'standard' brain.
%
% The first is shown in red, then green then blue.
%
% The blobs which are displayed are the integral of all transformed t
% values, exponentially decayed according to their depth. Voxels that
% are 10mm behind the surface have half the intensity of ones at the
% surface.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_render.m 3289 2009-07-27 15:28:24Z guillaume $

persistent rend;
if nargin<1, rend=[]; return; end
DEPTH=.2;
SVNrev = '$Rev: 3289 $';

if nargin<4, rendertype='display'; end
prevrend = struct('rendfile',rendfile,...
    'brt',brt,...
    'col',[]);

%-Parse arguments, get data if not passed as parameters
%==========================================================================
if nargin < 1
    spm('FnBanner',mfilename,SVNrev);
    spm('FigName','Results: render');

    num   = spm_input('Number of sets',1,'1 set|2 sets|3 sets',[1 2 3],1);

    for i = 1:num
        [SPM,xSPM] = spm_getSPM;
        dat(i)    = struct( 'XYZ',  xSPM.XYZ,...
                    't',    xSPM.Z',...
                    'mat',  xSPM.M,...
                    'dim',  xSPM.DIM);
    end
    showbar = 1;
else
    num     = length(dat);
    showbar = 0;
end

%-Get surface
%--------------------------------------------------------------------------
if nargin < 3
    [rendfile, sts] = spm_select(1,'mesh','Render file'); % .mat or .gii file
    if ~sts, return; end
end
prevrend.rendfile = rendfile;

if isempty(rend),
[p,f,e] = fileparts(rendfile);
loadgifti = false;
if strcmpi(e,'.mat')
    load(rendfile);
    if ~exist('rend','var') && ~exist('Matrixes','var')
        loadgifti = true;
    end
end
if ~strcmpi(e,'.mat') || loadgifti
    try
        rend = export(gifti(rendfile),'patch');
    catch
        error('\nCannot read  render file "%s".\n', rendfile);
    end
    if num == 1
        col = hot(256);
    else
        col = eye(3);
        if spm_input('Which colours?','!+1','b',{'RGB','Custom'},[0 1],1)
            for k = 1:num
                col(k,:) = uisetcolor(col(k,:),sprintf('Colour of blob set %d',k));
            end
        end
    end
    surf_rend(dat,rend,col);
    return
end
end

%-Get brightness & colours
%--------------------------------------------------------------------------
if nargin < 2
    brt = 1;
    if num==1
        brt = spm_input('Style',1,'new|old',[1 NaN], 1);
    end

    if isfinite(brt)
        brt = spm_input('Brighten blobs',1,'none|slightly|more|lots',[1 0.75 0.5 0.25], 1);
        col = eye(3);
        % ask for custom colours & get rgb values
        %------------------------------------------------------------------
        if spm_input('Which colours?','!+1','b',{'RGB','Custom'},[0 1],1)
            for k = 1:num
                col(k,:) = uisetcolor(col(k,:),sprintf('Colour of blob set %d',k));
            end
        end
    else
        col = [];
    end
elseif isfinite(brt) && isempty(prevrend.col)
    col = eye(3);
elseif isfinite(brt)  % don't need to check prevrend.col again
    col = prevrend.col;
else
    col = [];
end
prevrend.brt = brt;
prevrend.col = col;

%-Perform the rendering
%==========================================================================
spm('Pointer','Watch');

if ~exist('rend','var') % Assume old format...
    rend = cell(size(Matrixes,1),1);
    for i=1:size(Matrixes,1),
        rend{i}=struct('M',eval(Matrixes(i,:)),...
            'ren',eval(Rens(i,:)),...
            'dep',eval(Depths(i,:)));
        rend{i}.ren = rend{i}.ren/max(max(rend{i}.ren));
    end
end

if showbar, spm_progress_bar('Init', size(dat,1)*length(rend),...
            'Formatting Renderings', 'Number completed'); end
for i=1:length(rend),
    rend{i}.max=0;
    rend{i}.data = cell(size(dat,1),1);
    if issparse(rend{i}.ren),
        % Assume that images have been DCT compressed
        % - the SPM99 distribution was originally too big.
        d = size(rend{i}.ren);
        B1 = spm_dctmtx(d(1),d(1));
        B2 = spm_dctmtx(d(2),d(2));
        rend{i}.ren = B1*rend{i}.ren*B2';
        % the depths did not compress so well with
        % a straight DCT - therefore it was modified slightly
        rend{i}.dep = exp(B1*rend{i}.dep*B2')-1;
        idx=find(rend{i}.dep>200);temp=rend{i}.dep(idx);rend{i}.dep(idx)=nan;rend{i}.dep=convn(rend{i}.dep,ones(5)/25,'same');rend{i}.dep(idx)=temp;
        %rend{i}.ren(idx)=.75;
    end
    rend{i}.ren(rend{i}.ren>=1) = 1;
    rend{i}.ren(rend{i}.ren<=0) = 0;
    if showbar, spm_progress_bar('Set', i); end
end
if showbar, spm_progress_bar('Clear'); end

if showbar, spm_progress_bar('Init', length(dat)*length(rend),...
            'Making pictures', 'Number completed'); end

mx = zeros(length(rend),1)+eps;
mn = zeros(length(rend),1);

for j=1:length(dat),
    if ~isfield(dat(j),'VOL')
        XYZ = dat(j).XYZ;
        t   = dat(j).t;
        VOL=zeros(dat(j).dim);
        VOL(sub2ind(size(VOL),XYZ(1,:),XYZ(2,:),XYZ(3,:)))=t(:)';
    else
        VOL=dat(j).VOL;
    end
    dim = dat(j).dim;
    mat = dat(j).mat;

    for i=1:length(rend),

          idx=find(rend{i}.dep<200);
          [x,y]=ind2sub(size(rend{i}.dep),idx);
          z=rend{i}.dep(idx); 
          xyz=inv(rend{i}.M)*[x(:),y(:),z(:),ones(numel(x),1)]'; 
          X=zeros(size(rend{i}.dep));Z=X;
          if strcmp(rendertype,'display'),
              %for n1=1:10, k=.8+.2*(n1-1)/9;
              for n1=1:10, k=1-DEPTH*(n1-1)/9;
                  xyzt=inv(dat.mat)*diag([k,k,k,1])*xyz;
                  Xt=spm_sample_vol(VOL,xyzt(1,:),xyzt(2,:),xyzt(3,:),0)';
                  X(idx)=max(X(idx),Xt);
                  %X(idx)=X(idx)+Xt;
              end
          else
              for n1=1:10, k=1-DEPTH*(n1-1)/9;
                  xyzt=inv(dat.mat)*diag([k,k,k,1])*xyz;
                  Xt=spm_sample_vol(VOL,xyzt(1,:),xyzt(2,:),xyzt(3,:),0)';
                  Z(idx(Xt>0))=k;
                  X(idx)=max(X(idx),Xt);
                  idx=idx(Xt==0);
                  xyz=xyz(:,Xt==0);
              end
              if ~isempty(which('poly2mask')),
                  uX=unique(X(X>0)); X2=zeros(size(X));
                  for n1=1:numel(uX),
                      c=contourc(double(X==uX(n1)),[.5,.5]);
                      n2=1;while n2<size(c,2),
                          xy=c(:,n2+(1:c(2,n2)));
                          nxy=size(xy,2);
                          if nxy>20,
                              xy=convn([xy,xy],ones(1,20)/20,'same');
                              xy=xy(:,ceil(nxy/2):(nxy+ceil(nxy/2)-1));
                              x2=poly2mask(xy(1,:),xy(2,:),size(X,1),size(X,2));
                              X2(x2>0)=uX(n1);
                          end
                          n2=n2+c(2,n2)+1;
                          %plot(xy(1,:),xy(2,:),'.-');pause;
                      end
                  end
                  X=X2;
              end
%               ux=unique(X(X>0));
%               for n1=1:numel(ux),
%                   idx=find(X==ux(n1));
%                   Xt=zeros(size(X));Xt(idx)=1;Xt=convn(Xt,ones(11)/121,'same');
%                   X(idx)=0;X(Xt>.5)=ux(n1);
%               end
          end
%         % transform from Talairach space to space of the rendered image
%         %------------------------------------------------------------------
%         M1  = rend{i}.M*mat;
%         zm  = sum(M1(1:2,1:3).^2,2).^(-1/2);
%         M2  = diag([zm' 1 1]);
%         M  = M2*M1;
%         cor = [1 1 1 ; dim(1) 1 1 ; 1 dim(2) 1; dim(1) dim(2) 1 ;
%                1 1 dim(3) ; dim(1) 1 dim(3) ; 1 dim(2) dim(3); dim(1) dim(2) dim(3)]';
%         tcor= M(1:3,1:3)*cor + M(1:3,4)*ones(1,8);
%         off = min(tcor(1:2,:)');
%         M2  = spm_matrix(-off+1)*M2;
%         M  = M2*M1;
%         xyz = (M(1:3,1:3)*XYZ + M(1:3,4)*ones(1,size(XYZ,2)));
%         d2  = ceil(max(xyz(1:2,:)'));
% 
%         % Calculate 'depth' of values
%         %------------------------------------------------------------------
%         dep = spm_slice_vol(rend{i}.dep,spm_matrix([0 0 1])*inv(M2),d2,1);
%         z1  = dep(round(xyz(1,:))+round(xyz(2,:)-1)*size(dep,1));
% 
%         if ~isfinite(brt), msk = find(xyz(3,:) < (z1+20) & xyz(3,:) > (z1-5));
%         else,      msk = find(xyz(3,:) < (z1+60) & xyz(3,:) > (z1-5)); end
% 
%         if ~isempty(msk),
% 
%             % Generate an image of the integral of the blob values.
%             %--------------------------------------------------------------
%             xyz = xyz(:,msk);
%             if ~isfinite(brt), t0  = t(msk);
%             else,   dst = xyz(3,:) - z1(msk);
%                 dst = max(dst,0);
%                 t0  = t(msk).*exp((log(0.5)/10)*dst)';
%             end
%             X0  = full(sparse(round(xyz(1,:)), round(xyz(2,:)), t0, d2(1), d2(2)));
%             hld = 1; if ~isfinite(brt), hld = 0; end
%             X   = spm_slice_vol(X0,spm_matrix([0 0 1])*M2,size(rend{i}.dep),hld);
%             msk = find(X<0);
%             X(msk) = 0;
%         else
%             X = zeros(size(rend{i}.dep));
%         end
% 
        % Brighten the blobs
        %------------------------------------------------------------------
        if strcmp(rendertype,'display'),
        if isfinite(brt), X = X.^brt; end
        end
        
        mx(j) = max([mx(j) max(max(X))]);
        mn(j) = min([mn(j) min(min(X))]);

        rend{i}.data{j} = X;

        if showbar, spm_progress_bar('Set', i+(j-1)*length(rend)); end
    end
end

mxmx = max(mx);
mnmn = min(mn);

if showbar, spm_progress_bar('Clear'); end
if ~nargout,
Fgraph = spm_figure('GetWin','Graphics');
spm_results_ui('Clear',Fgraph);

nrow = ceil(length(rend)/2);
if showbar, hght = 0.95; else, hght = 0.5; end
% subplot('Position',[0, 0, 1, hght]);
ax=axes('Parent',Fgraph,'units','normalized','Position',[0, 0, 1, hght],'Visible','off');
image(0,'Parent',ax);
set(ax,'YTick',[],'XTick',[]);
end

if ~isfinite(brt),
    % Old style split colourmap display.
    %----------------------------------------------------------------------
    load Split;
    colormap(split);
    for i=1:length(rend),
        ren = rend{i}.ren;
        X   = (rend{i}.data{1}-mnmn)/(mxmx-mnmn);
        msk = find(X);
        ren(msk) = X(msk)+(1+1.51/64);
        ax=axes('Parent',Fgraph,'units','normalized',...
            'Position',[rem(i-1,2)*0.5, floor((i-1)/2)*hght/nrow, 0.5, hght/nrow],...
            'Visible','off');
        image(ren*64,'Parent',ax);
        set(ax,'DataAspectRatio',[1 1 1], ...
            'PlotBoxAspectRatioMode','auto',...
            'YTick',[],'XTick',[],'XDir','normal','YDir','normal');
    end
else
    % Combine the brain surface renderings with the blobs, and display using
    % 24 bit colour.
    %----------------------------------------------------------------------
    if strcmp(rendertype,'display'),
        for i=1:length(rend),
            ren = rend{i}.ren;
            X = cell(3,1);
            for j=1:length(rend{i}.data),
                X{j} = (rend{i}.data{j}-mnmn)/(mxmx-mnmn);
            end
            for j=(length(rend{i}.data)+1):3
                X{j}=zeros(size(X{1}));
            end
            
            rgb = zeros([size(ren) 3]);
            for k=1:3,if isempty(X{k}),X{k}=0;end;end
            tmp = ren.*max(1-X{1}-X{2}-X{3},0);
            for k = 1:3
                rgb(:,:,k) = tmp + X{1}*col(1,k) + X{2}*col(2,k) +X{3}*col(3,k);
            end
            rgb(rgb>1) = 1;
            if nargout,
                RGB{end+1}=rgb;
                RANGERGB=[mnmn,mxmx];
            else,
                ax=axes('Parent',Fgraph,'units','normalized',...
                    'Position',[rem(i-1,2)*0.5, floor((i-1)/2)*hght/nrow, 0.5, hght/nrow],...
                    'nextplot','add', ...
                    'Visible','off');
                image(rgb,'Parent',ax);
                set(ax,'DataAspectRatio',[1 1 1], ...
                    'PlotBoxAspectRatioMode','auto',...
                    'YTick',[],'XTick',[],...
                    'XDir','normal','YDir','normal');
            end
        end
    else
        if nargout,
            for i=1:length(rend),
                RGB{end+1}=rend{i}.data{1};
            end
        end
    end
end

spm('Pointer','Arrow');

%==========================================================================
% function surf_rend(dat,rend,col)
%==========================================================================
function surf_rend(dat,rend,col)

%-Setup figure and axis
%--------------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');
spm_results_ui('Clear',Fgraph);
rdr = get(Fgraph,'Renderer');
set(Fgraph,'Renderer','OpenGL');

ax0 = axes(...
    'Parent',Fgraph,...
    'Units','normalized',...
    'Color', [1 1 1],...
    'XTick',[],...
    'YTick',[],...
    'Position',[-0.05, -0.05, 1.05, 0.555]);

ax = axes(...
    'Parent',Fgraph,...
    'Units','normalized',...
    'Position',[0.05, 0.05, 0.9, 0.4],...
    'Visible','off');

%-Project data onto surface mesh
%--------------------------------------------------------------------------
v = spm_mesh_project(rend, dat);

%-Compute mesh curvature texture
%--------------------------------------------------------------------------
curv = spm_mesh_curvature(rend) > 0;
curv = 0.5 * repmat(curv,1,3) + 0.3 * repmat(~curv,1,3);
%curv(:)=0.5;
%-Combine projected data and mesh curvature
%--------------------------------------------------------------------------
cdat = zeros(size(v,2),3);
if any(v(:))
    if size(col,1)>3
        cdat = squeeze(ind2rgb(floor(v(:)/max(v(:))*size(col,1)),col));
    else
        m = max(v(:));
        for i=1:size(v,1)
            cdat = cdat + v(i,:)'/m * col(i,:);
        end
    end
end
cdat = repmat(~any(v,1),3,1)' .* curv + repmat(any(v,1),3,1)' .* cdat;

%-Display the surface mesh with texture
%--------------------------------------------------------------------------
hp = patch(rend, 'Parent',ax,...
    'FaceVertexCData',cdat, ...
    'FaceColor', 'interp', ...
    'EdgeColor', 'none',...
    'FaceLighting', 'phong',...
    'SpecularStrength' ,0.7, 'AmbientStrength', 0.1,...
    'DiffuseStrength', 0.7, 'SpecularExponent', 10,...
    'DeleteFcn', {@mydeletefcn,Fgraph,rdr});

set(Fgraph,'CurrentAxes',ax);
view(ax,[-90 0]);
axis(ax,'image');

l = camlight; set(l,'Parent',ax);
material(Fgraph,'dull');
setappdata(ax,'camlight',l);

%-Setup context menu
%--------------------------------------------------------------------------
r = rotate3d(ax);
set(r,'enable','off');
cmenu = uicontextmenu;
c1 = uimenu(cmenu, 'Label', 'Inflate', 'Interruptible','off', 'Callback', @myinflate);
setappdata(c1,'patch',hp);
setappdata(c1,'axis',ax);
c2 = uimenu(cmenu, 'Label', 'Connected Components', 'Interruptible','off');
C = spm_mesh_label(hp);
setappdata(c2,'patch',hp);
setappdata(c2,'cclabel',C);
for i=1:length(unique(C))
    uimenu(c2,'Label',sprintf('Component %d',i), 'Checked','on', 'Callback', @mycclabel);
end
c3 = uimenu(cmenu, 'Label', 'Rotate', 'Checked','on','Separator','on','Callback', @myswitchrotate);
setappdata(c3,'rotate3d',r);
c4 = uimenu(cmenu, 'Label', 'View');
setappdata(c4,'axis',ax);
uimenu(c4,'Label','Go to Y-Z view (right)',  'Callback', {@myview, [90 0]});
uimenu(c4,'Label','Go to Y-Z view (left)',   'Callback', {@myview, [-90 0]});
uimenu(c4,'Label','Go to X-Y view (top)',    'Callback', {@myview, [0 90]});
uimenu(c4,'Label','Go to X-Y view (bottom)', 'Callback', {@myview, [-180 -90]});
uimenu(c4,'Label','Go to X-Z view (front)',  'Callback', {@myview, [-180 0]});
uimenu(c4,'Label','Go to X-Z view (back)',   'Callback', {@myview, [0 0]});
c5 = uimenu(cmenu, 'Label', 'Transparency');
setappdata(c5,'patch',hp);
uimenu(c5,'Label','0%',  'Checked','on',  'Callback', @mytransparency);
uimenu(c5,'Label','20%', 'Checked','off', 'Callback', @mytransparency);
uimenu(c5,'Label','40%', 'Checked','off', 'Callback', @mytransparency);
uimenu(c5,'Label','60%', 'Checked','off', 'Callback', @mytransparency);
uimenu(c5,'Label','80%', 'Checked','off', 'Callback', @mytransparency);
c6 = uimenu(cmenu, 'Label', 'Background Color');
setappdata(c6,'fig',ax0);
uimenu(c6,'Label','White', 'Callback', {@mybackgroundcolor, [1 1 1]});
uimenu(c6,'Label','Black', 'Callback', {@mybackgroundcolor, [0 0 0]});
uimenu(c6,'Label','Custom...', 'Callback', {@mybackgroundcolor, []});
c7 = uimenu(cmenu, 'Label', 'Save As...','Separator','on','Callback', @mysave);
setappdata(c7,'patch',hp);
setappdata(c7,'fig',Fgraph);
setappdata(c7,'axis',ax);
setappdata(c7,'ax0',ax0);
try, set(r,'uicontextmenu',cmenu); end
try, set(hp,'uicontextmenu',cmenu); end
set(r,'enable','on');
set(r,'ActionPostCallback',@mypostcallback);
try
    setAllowAxesRotate(r, setxor(findobj(Fgraph,'Type','axes'),ax), false);
end
    
%-Register with MIP
%--------------------------------------------------------------------------
try % meaningless when called outside spm_results_ui
    hReg = spm_XYZreg('FindReg',spm_figure('GetWin','Interactive'));
    xyz  = spm_XYZreg('GetCoords',hReg);
    hs   = mydispcursor('Create',ax,dat.mat,xyz);
    spm_XYZreg('Add2Reg',hReg,hs,@mydispcursor);
end
    
%==========================================================================
function myinflate(obj,evd)
spm_mesh_inflate(getappdata(obj,'patch'),Inf,1);
axis(getappdata(obj,'axis'),'image');

%==========================================================================
function mycclabel(obj,evd)
C = getappdata(get(obj,'parent'),'cclabel');
F = get(getappdata(get(obj,'parent'),'patch'),'Faces');
ind = sscanf(get(obj,'Label'),'Component %d');
V = get(getappdata(get(obj,'parent'),'patch'),'FaceVertexAlphaData');
Fa = get(getappdata(get(obj,'parent'),'patch'),'FaceAlpha');
if ~isnumeric(Fa)
    if ~isempty(V), Fa = max(V); else Fa = 1; end
    if Fa == 0, Fa = 1; end
end
if isempty(V) || numel(V) == 1
    Ve = get(getappdata(get(obj,'parent'),'patch'),'Vertices');
    if isempty(V) || V == 1
        V = Fa * ones(size(Ve,1),1);
    else
        V = zeros(size(Ve,1),1);
    end
end
if strcmpi(get(obj,'Checked'),'on')
    V(reshape(F(C==ind,:),[],1)) = 0;
    set(obj,'Checked','off');
else
    V(reshape(F(C==ind,:),[],1)) = Fa;
    set(obj,'Checked','on');
end
set(getappdata(get(obj,'parent'),'patch'), 'FaceVertexAlphaData', V);
if all(V)
    set(getappdata(get(obj,'parent'),'patch'), 'FaceAlpha', Fa);
else
    set(getappdata(get(obj,'parent'),'patch'), 'FaceAlpha', 'interp');
end
    
%==========================================================================
function myswitchrotate(obj,evd)
if strcmpi(get(getappdata(obj,'rotate3d'),'enable'),'on')
    set(getappdata(obj,'rotate3d'),'enable','off');
    set(obj,'Checked','off');
else
    set(getappdata(obj,'rotate3d'),'enable','on');
    set(obj,'Checked','on');
end

%==========================================================================
function myview(obj,evd,varargin)
view(getappdata(get(obj,'parent'),'axis'),varargin{1});
axis(getappdata(get(obj,'parent'),'axis'),'image');
camlight(getappdata(getappdata(get(obj,'parent'),'axis'),'camlight'));

%==========================================================================
function mytransparency(obj,evd)
t = 1 - sscanf(get(obj,'Label'),'%d%%') / 100;
set(getappdata(get(obj,'parent'),'patch'),'FaceAlpha',t);
set(get(get(obj,'parent'),'children'),'Checked','off');
set(obj,'Checked','on');

%==========================================================================
function mybackgroundcolor(obj,evd,varargin)
if isempty(varargin{1})
    c = uisetcolor(getappdata(get(obj,'parent'),'fig'), ...
        'Pick a background color...');
    if numel(c) == 1, return; end
else
    c = varargin{1};
end
set(getappdata(get(obj,'parent'),'fig'),'Color',c);

%==========================================================================
function mysave(obj,evd)
[filename, pathname, filterindex] = uiputfile({...
    '*.gii' 'GIfTI files (*.gii)'; ...
    '*.png' 'PNG files (*.png)';...
    '*.dae' 'Collada files (*.dae)'}, 'Save as');
if ~isequal(filename,0) && ~isequal(pathname,0)
    [pth,nam,ext] = fileparts(filename);
    switch ext
        case '.gii'
            filterindex = 1;
        case '.png'
            filterindex = 2;
        case '.dae'
            filterindex = 3;
        otherwise
            switch filterindex
                case 1
                    filename = [filename '.gii'];
                case 2
                    filename = [filename '.png'];
                case 3
                    filename = [filename '.dae'];
            end
    end
    switch filterindex
        case 1
            g = gifti;
            g.vertices = get(getappdata(obj,'patch'),'Vertices');
            g.faces = get(getappdata(obj,'patch'),'Faces');
            g.cdata = get(getappdata(obj,'patch'),'FaceVertexCData');
            save(g,fullfile(pathname, filename));
        case 2
            ax = getappdata(obj,'axis');
            u = get(ax,'units');
            set(ax,'units','pixels');
            p = get(ax,'Position');
            r = get(getappdata(obj,'fig'),'Renderer');
            h = figure('Position',p+[0 0 10 10], ...
                'InvertHardcopy','off', ...
                'Color',get(getappdata(obj,'ax0'),'Color'), ...
                'Renderer',r);
            copyobj(getappdata(obj,'axis'),h);
            set(ax,'units',u);
            set(get(h,'children'),'visible','off');
            %a = get(h,'children');
            %set(a,'Position',get(a,'Position').*[0 0 1 1]+[10 10 0 0]);       
            print(h, '-dpng', '-opengl', fullfile(pathname, filename));
            close(h);
            set(getappdata(obj,'fig'),'renderer',r);
        case 3
            g = gifti;
            g.vertices = get(getappdata(obj,'patch'),'Vertices');
            g.faces = get(getappdata(obj,'patch'),'Faces');
            g.cdata = get(getappdata(obj,'patch'),'FaceVertexCData');
            save(g,fullfile(pathname, filename),'collada');
    end
end

%==========================================================================
function mypostcallback(obj,evd)
try, camlight(getappdata(evd.Axes,'camlight')); end

%==========================================================================
function mydeletefcn(obj,evd,varargin)
try, rotate3d(get(obj,'parent'),'off'); end
set(varargin{1},'Renderer',varargin{2});

%==========================================================================
function varargout = mydispcursor(varargin)

switch lower(varargin{1})
    %======================================================================
    case 'create'
    %======================================================================
    % hMe = mydispcursor('Create',ax,M,xyz)
    ax  = varargin{2};
    M   = varargin{3};
    xyz = varargin{4};
    
    [X,Y,Z] = sphere;
    vx = sqrt(sum(M(1:3,1:3).^2));
    X = X*vx(1) + xyz(1);
    Y = Y*vx(2) + xyz(2);
    Z = Z*vx(3) + xyz(3);
    hold(ax,'on');
    hs = surf(X,Y,Z,'parent',ax,...
        'EdgeColor','none','FaceColor',[1 0 0],'FaceLighting', 'phong');
    set(hs,'UserData',xyz);
    
    varargout = {hs};
    
    %=======================================================================
    case 'setcoords'    % Set co-ordinates
    %=======================================================================
    % [xyz,d] = mydispcursor('SetCoords',xyz,hMe,hC)
    hMe  = varargin{3};
    pxyz = get(hMe,'UserData');
    xyz  = varargin{2};
    
    set(hMe,'XData',get(hMe,'XData') - pxyz(1) + xyz(1));
    set(hMe,'YData',get(hMe,'YData') - pxyz(2) + xyz(2));
    set(hMe,'ZData',get(hMe,'ZData') - pxyz(3) + xyz(3));
    set(hMe,'UserData',xyz);
    
    varargout = {xyz,[]};
    
    %=======================================================================
    otherwise
    %=======================================================================
    error('Unknown action string')

end
