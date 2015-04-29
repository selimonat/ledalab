function leda_batchanalysis(varargin)

global leda2

%parse batch-mode arguments and check their validity
[valid_options, pathname, open_datatype, filter_settings, downsample_factor, smooth_settings, analysis_method, do_optimize, ...
    export_era_settings, export_scrlist_settings,export_eta_settings, do_save_overview] = parse_arguments(varargin{:});

if ~valid_options || ~(downsample_factor > 1 || analysis_method || do_optimize || any(export_era_settings) || any(export_scrlist_settings) || do_save_overview) %invalid option or no option
    disp('No valid operations for Batch-mode defined.')
    return;
end

%this is a cell array that contains the list of pathnames to the files to be analyzed.
dirL  = varargin{1};
nFile = length(dirL);

add2log(1,['Starting Ledalab batch for ',pathname,' (',num2str(nFile),' file/s)'],1,0,0,1)
for iFile = 1:nFile

leda2.current.batchmode.file = [];
leda2.current.batchmode.command.pathname = pathname;
leda2.current.batchmode.command.datatype = open_datatype;
leda2.current.batchmode.command.filter_settings = filter_settings;
leda2.current.batchmode.command.downsample = downsample_factor;
leda2.current.batchmode.command.smooth = smooth_settings;
leda2.current.batchmode.command.method = analysis_method;
leda2.current.batchmode.command.optimize = do_optimize;
leda2.current.batchmode.command.overview = do_save_overview;
leda2.current.batchmode.command.export_era = export_era_settings;
leda2.current.batchmode.command.export_scrlist = export_scrlist_settings;
leda2.current.batchmode.start = datestr(now, 21);
leda2.current.batchmode.version = leda2.intern.version;
leda2.current.batchmode.settings = leda2.set;
tic

%ledalab is designed to run across different files located in
%a single folder during batch mode, this is largely incompatible with our
%data tree.
    filename = dirL{iFile};
    [pathname filename ext]= fileparts(filename);
    filename = [filename ext];
    %
    leda2.current.batchmode.file(iFile).name = filename;
    disp(' '); add2log(1,['Batch-Analyzing ',filename],1,0,0,1)
    
%     try
        %Open
        if strcmp(open_datatype,'leda')
            open_ledafile(0, pathname, filename);
        else
            import_data(open_datatype, pathname, filename);
        end
        if ~leda2.current.fileopen_ok
            disp('Unable to open file!');
            continue;
        end
        
        %Filter, MB: 14.05.2014
        if filter_settings(1) > 0
            leda_filter(filter_settings);
        end
        
        %Downsample
        if downsample_factor > 1
            leda_downsample(downsample_factor, 'mean');  %MB 11.06.2013
        end
        
        %Smooth
        if iscell(smooth_settings)
            if strcmpi(smooth_settings{1},'adapt')
                adaptive_smoothing;
            else
                smooth_data(smooth_settings{2}, smooth_settings{1})
            end
        end
        
        %Analysis
        if analysis_method > 0
            delete_fit;
            if analysis_method == 1
                sdeco(do_optimize);
            elseif analysis_method == 2
                nndeco(do_optimize);
            end
            leda2.current.batchmode.file(iFile).tau = leda2.analysis.tau;
            leda2.current.batchmode.file(iFile).error = leda2.analysis.error;
        end
        
        %Export ERA
        if any(export_era_settings)
            leda2.set.export.SCRstart = export_era_settings(1);
            leda2.set.export.SCRend = export_era_settings(2);
            leda2.set.export.SCRmin = export_era_settings(3);
            if length(export_era_settings) > 3
                leda2.set.export.savetype = export_era_settings(4);
            else
                leda2.set.export.savetype = 1;
            end
            export_era('savePeaks')
        end
        
        %Export ETA
        if any(export_eta_settings)
            %right now the temporal window parameters are taken from the era
            %analysis, in case era analysis is not required, this will give
            %an error.
            leda2.set.split.start  = export_era_settings(1);
            leda2.set.split.end    = export_era_settings(2);            
            %
            leda2.set.split.plot     = 0;
            for variables2split = {'driver' 'tonicDriver' 'phasicData' 'tonicData' 'phasicDriverRaw'}
                leda2.set.split.variable = variables2split{1};
                leda_split('split');
            end            
        end
        
        %Export Scrlist
        if any(export_scrlist_settings)
            leda2.set.export.SCRmin = export_scrlist_settings(1);
            if length(export_scrlist_settings) > 1
                leda2.set.export.savetype = export_scrlist_settings(2);
            else
                leda2.set.export.savetype = 1;
            end
            export_scrlist('saveList')
        end
        
        %Save
        if do_save_overview
            analysis_overview;
        end
        
        if downsample_factor > 0 || analysis_method  || iscell(smooth_settings)
            save_ledafile([leda2.file.filename(1:end-4) '_results']);
        end
        
%     catch
%         add2log(1,'ERROR !!!',1,0,0,1)
%     end
    
end

leda2.current.batchmode.processing_time = toc;
protocol = leda2.current.batchmode;
save([pathname,filesep,'batchmode_protocol'],'protocol');



function [valid_options, wdir, open_datatype, filter_settings, downsample_factor, smooth_settings, analysis_method, do_optimize, ...
    export_era_settings, export_scrlist_settings, export_eta_settings, do_save_overview] = parse_arguments(varargin)

wdir = [];%we don't have a single working directory so it is empty.

valid_options = 1;
%default options
open_datatype = 'leda'; %open
filter_settings = [0 0];
downsample_factor = 0;
smooth_settings = 0;
analysis_method = 0;
do_optimize = 2;
%do_export_scr = 0;
export_era_settings = [0 0 0 0];
export_scrlist_settings = [0 0];
do_save_overview = 0;

%valid_datatypeL = {'leda','mat','text','cassylab','biopac','biopacmat','biotrace','visionanalyzer','userdef'};
%datatype_extL = {'*.mat','*.mat','*.txt','*.txt','*.txt','',''};

if nargin > 1
    vars = varargin(2:end);
    
    while length(vars) >= 2
        thisvar = vars(1:2);
        vars = vars(3:end);
        
        option_name = thisvar{1};
        option_arg = thisvar{2};
        
        switch option_name
            case 'open',
                %if ischar(option_arg) && any(strcmp(option_arg, valid_datatypeL))
                open_datatype = option_arg;
                %wdir = wdir(1:end-5);  %remove default value *.mat
                %wdir = [wdir(1:end-5), datatype_extL{strcmp(option_arg, valid_datatypeL)}];
                %else
                %    disp(['Unknown datatype: ',option_arg])
                %    return;
                %end
                
            case 'filter'
                if isnumeric(option_arg)
                    filter_settings = option_arg;
                else
                    valid_options = 0;
                    disp('Filter settings require 2 numeric arguments (filter order, and lower cutoff, e.g. [1 5])')
                    return;
                end
                
            case 'downsample'
                if isnumeric(option_arg)
                    downsample_factor = option_arg;
                else
                    valid_options = 0;
                    disp('Downsample option requires numeric argument (downsample factor)')
                    return;
                end
                
            case 'smooth'       %by Christoph Berger
                if iscell(option_arg) && any(strcmpi(option_arg{1},{'hann','mean','gauss','adapt'}))
                    smooth_settings = option_arg;
                else
                    valid_options = 0;
                    disp('Smooth option requires cell; first argument is ''hann'', ''mean'', ''gauss'', or ''adapt'', second argument is width')
                    return;
                end
                
            case 'analyze'
                if any(strcmpi(option_arg,{'CDA','DDA'}))
                    analysis_method = find(strcmpi(option_arg, {'CDA','DDA'}));  % 1 = CDA, 2 = DDA
                else
                    valid_options = 0;
                    disp('Method option should either bei ''CDA'' or ''DDA''.')
                    return;
                end
                
            case 'optimize'
                if isnumeric(option_arg)
                    do_optimize = option_arg;
                else
                    valid_options = 0;
                    disp('Optimize option requires numeric argument (# of initial values for optimization)')
                    return;
                end
                
            case 'export_era'
                if isnumeric(option_arg) && any(option_arg) && length(option_arg) >= 3
                    export_era_settings = option_arg;
                else
                    valid_options = 0;
                    disp('Export requires numeric argument (respwin_start respwin_end amp_threshold [filetype]) ')
                    return;
                end
                
            case 'export_scrlist'
                if isnumeric(option_arg) && any(option_arg) && length(option_arg) >= 1
                    export_scrlist_settings = option_arg;
                else
                    valid_options = 0;
                    disp('Export requires numeric argument (amp_threshold [filetype]) ')
                    return;
                end
            
            case 'export_eta'
                if ~isempty(option_arg)
                    export_eta_settings = option_arg;
                else
                    valid_options = 0;
                    disp('Export requires numeric argument (amp_threshold [filetype]) ')
                    return;
                end    
                
            case 'overview'
                if isnumeric(option_arg)
                    do_save_overview = option_arg;
                else
                    valid_options = 0;
                    disp('Overview option requires numeric argument (1 = yes, 0 = no)')
                    return;
                end
                
            otherwise
                valid_options = 0;
                disp(['Could not parse batch-mode option: ',option_name])
                
        end %switch
        
    end %while
    
end %if nargin > 1



function analysis_overview
global leda2

t = leda2.data.time.data;
analysis = leda2.analysis;
events = leda2.data.events;
%correct for extended data range of older versions
if leda2.file.version < 3.12
    n_offset = length(analysis.time_ext);
    remainder = analysis.remainder(n_offset+1:end);
    driver = leda2.analysis.driver(n_offset+1:end);
else
    remainder = analysis.remainder;
    driver = leda2.analysis.driver;
end


figure('Units','normalized','Position',[0 0.05 1 .9],'MenuBar','none','NumberTitle','off');

%Decomposition
subplot(2,1,1);
cla; hold on;
title('SC Data')
if leda2.file.version < 3.12 || strcmp(leda2.analysis.method,'nndeco')
    if length(analysis.phasicRemainder) * length(analysis.tonicData) < 4*10^6
        for i = 2:length(analysis.phasicRemainder)
            plot(t, analysis.tonicData + analysis.phasicRemainder{i})
        end
    end
end

plot(t, leda2.data.conductance.data, 'k','Linewidth',2);
plot(t, analysis.tonicData + analysis.phasicData,'k:','Linewidth',2);
plot(t, analysis.tonicData,'Color',[.6 .6 .6],'Linewidth',2);
%plot(analysis.groundtime, analysis.groundlevel,'o','LineWidth',2,'MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[.9 .9 .9],'MarkerSize',3)

%ensure minimum scaling of 2 muS
yl = get(gca,'YLim');
if abs(diff(yl)) < 2
    yl(2) = yl(1) + 2;
end
set(gca, 'XLim', [t(1), t(end)],'Ylim',yl)
%Events
yl = ylim;
for i = 1:events.N
    plot([events.event(i).time, events.event(i).time], yl, 'r')
end
set(gca,'YLim',yl)
if strcmp(analysis.method,'nndeco')
    l = legend('SC Data','Decomposition Fit','Tonic Data', sprintf('tau = %4.2f, %4.2f,  dist0 = %4.4f',analysis.tau, analysis.dist0), sprintf('RMSE = %4.2f', analysis.error.RMSE));
else
    l = legend('SC Data','Decomposition Fit','Tonic Data', sprintf('tau = %4.2f, %4.2f',analysis.tau), sprintf('RMSE = %4.2f', analysis.error.RMSE));
end
set(l, 'FontSize',8,'Location','NorthEast');
xlabel('Time [s]'); ylabel('[\muS]')
legend boxoff;

%Driver
subplot(2,1,2);
cla; hold on;
title('Phasic Driver')
plot(t, driver,'k','LineWidth',1);
plot(t, -2*remainder,'b','LineWidth',1);
set(gca, 'XLim', [t(1), t(end)], 'YLim', [min(min(driver), min(-2*remainder))*1.2, max(1, max(driver)*1.2)])
%Events
yl = ylim;
for i = 1:events.N
    plot([events.event(i).time, events.event(i).time], yl, 'r')
end
set(gca,'YLim',yl)
l = legend('Driver', 'Remainder', sprintf('Error-compound = %5.2f',analysis.error.compound), sprintf('Error-discr = %4.2f,  %4.2f', analysis.error.discreteness), sprintf('Error-neg = %4.2f', analysis.error.negativity)); %SouthOutside
set(l, 'FontSize',8,'Location','NorthEast');
xlabel('Time [s]'); ylabel('[\muS]')

legend boxoff;
subject_file = leda2.file.filename;
supertitle(subject_file, 1, 'interpreter','none');

SaveFigure([leda2.file.filename(1:end-4) '.png']);

drawnow;
pause(0.5)
close(gcf);
drawnow;


