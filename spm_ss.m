function varargout=spm_ss(varargin)
% SPM_SS Subject-specific analysis GUI
%
% see SPM_SS_DESIGN, SPM_SS_ESTIMATE, SPM_SS_RESULTS, SPM_SS_DISPLAY
%

% other:
% spm_ss qa <foldername>
%   creates and displays QA plots from QA_*.nii files generated during first-level estimation step
%

spm_ss_ver='18.a';

if nargin>0, 
    %disp(varargin); 
    if ischar(varargin{1}),
        switch(lower(varargin{1}))
            case 'init'
                if ~isempty(which('evlab17')), evlab17(varargin{:}); 
                else fprintf('warning: no evlab17 package found. SPM/CONN/SPM_SS version-checks not performed\n');
                end
            case {'ver','version'}
                varargout={spm_ss_ver};
            case {'qa','qacreate'}
                if isempty(which('conn'))&&~isempty(which('evlab17')), evlab17 init silent; end
                assert(~isempty(which('conn')),'no CONN toolbox found. Please install/addpath CONN to create QA plots');
                if usejava('awt')
                    if numel(varargin)<2||isempty(varargin{2}), cwd=pwd;
                    else cwd=conn_prepend('',varargin{2},'');
                    end
                    f=conn_dir(fullfile(cwd,'QA_effects.*.nii'),'-ls');
                    next={'_localizer.','_parcels.'};
                    for n=1:numel(f)
                        for inext=1:numel(next)
                            fh=conn_slice_display(regexprep(f{n},'_effects\.',next{inext}),f{n},'',[0 -inf]);
                            fh('contour_transparency',1);
                            fh('title',regexprep(f{n},'^.*QA_effects\.(.*?)\.nii$','$1'));
                            fh('multisliceset',1,100,6);
                            fh('togglegui',1);
                            filename=conn_prepend('',regexprep(f{n},'_effects\.',next{inext}),'_plot.jpg');
                            fh('print',filename,'-r150','-nopersistent','-nogui');
                            state=fh('getstate');
                            conn_args={'slice_display',state};
                            save(conn_prepend('',filename,'.mat'),'conn_args');
                            fh('close');
                        end
                    end
                    if strcmpi(varargin{1},'qa'), conn_qaplotsexplore(cwd); end
                else fprintf('warning: unable to generate QA plots (possibly missing graphic display capabilities). Please use "spm_ss qa ''%s''" syntax to create these plots at a later time\n',cwd);
                end
            otherwise
                if ischar(varargin{1})&&~isempty(which(sprintf('spm_ss_%s',varargin{1}))),
                    fh=eval(sprintf('@spm_ss_%s',varargin{1}));
                    if ~nargout, feval(fh,varargin{2:end});
                    else [varargout{1:nargout}]=feval(fh,varargin{2:end});
                    end
                else
                    disp(sprintf('unrecognized option %s or spm_ss_%s function',varargin{1},varargin{1}));
                end
                
        end
    elseif ishandle(varargin{1})&&numel(varargin)>=3
        switch(varargin{3}),
            case 'cv',
                spm_ss_crossvalidate_sessions;
            case 'lm',
                spm_ss_createlocalizermask;
            case 'surface',
                spm_ss_display;
            case 'preparePSTH',
                spm_ss_preparePSTH;
            case 'design',
                spm_ss_design;
            case 'estimate',
                spm_ss_estimate;
            case 'results',
                [ss,Ic]=spm_ss_selectcontrast;
                ss=spm_ss_contrast(ss,Ic,0);
                spm_ss_results(ss,Ic);
        end
    end
else
    h=findobj('tag',mfilename);
    if isempty(h),
        data.handles.fig=figure('units','norm','position',[.2,.5,.22,.2],'color','w','tag',mfilename,'menubar','none','numbertitle','off','name',[mfilename,' ',spm_ss_ver]);
        uicontrol('units','norm','position',[0,.8,1,.2],'style','frame','backgroundcolor','k');
        data.handles.uimenu=uimenu(data.handles.fig,'label','tools');
        data.handles.uimenu1=uimenu(data.handles.uimenu,'label','Creates Cross-validated contrasts','callback',{@spm_ss,'cv'});
        data.handles.uimenu2=uimenu(data.handles.uimenu,'label','Creates Localizer masks','callback',{@spm_ss,'lm'});
        data.handles.uimenu3=uimenu(data.handles.uimenu,'label','Renders subject-specific activations on brain surface','callback',{@spm_ss,'surface'});
        data.handles.uimenu4=uimenu(data.handles.uimenu,'label','Prepares PSTH analyses','callback',{@spm_ss,'preparePSTH'});
        data.handles.txt=uicontrol('units','norm','position',[.1,.85,.8,.1],'style','text','string','Subject-specific analyses','backgroundcolor','k','foregroundcolor','w','fontweight','bold','horizontalalignment','center');
        data.handles.button_design=uicontrol('units','norm','position',[.1,.525,.8,.2],'style','pushbutton','string','Specify 2nd-level','callback',{@spm_ss,'design'},'tooltipstring','GLM setup for subject-specific analyses');
        data.handles.button_estimate=uicontrol('units','norm','position',[.1,.325,.8,.2],'style','pushbutton','string','Estimate','callback',{@spm_ss,'estimate'},'tooltipstring','estimates parameters of a specified model');
        data.handles.button_results=uicontrol('units','norm','position',[.1,.075,.8,.2],'style','pushbutton','string','Results','callback',{@spm_ss,'results'},'tooltipstring','inference and regional responses etc.');
        set(data.handles.fig,'userdata',data);
    else
        data=get(h,'userdata');
    end
end


