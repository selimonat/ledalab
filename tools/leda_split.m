function leda_split(action)
% LEDA_SPLIT(ACTION)
%
% LEDA_SPLIT splits continuous data into epochs according to the time of
%   events. Each category of events is used to define a condition. The
%   epoched and averaged data can be plotted with or without the standard
%   error of the mean.
%
% NOTE: To get the same number of samples for each epoche, the time stamps
%   of the events are moved to the time of the closest sample. If sampling
%   rates are extremely low and the times of the events are NOT in
%   accordance with the times the samples were recorded, the resulting mean
%   (and standard error) might be distorted.
%
% analysis.driver and analysis.phasicData is available, 2012-04-16.

% Please adress questions, comments, bug reports related to this function to christoph.huber-huber@univie.ac.at
% Christoph Huber-Huber, 2012.



if nargin < 1,
    action = 'start';
end

switch action,
    case 'start', start;
    case 'take_settings', take_settings;
    case 'split', split;
    case 'select_all', select_all;
    case 'deselect_select_all', deselect_select_all;
end

end


function start
global leda2

if ~leda2.file.open
    if leda2.intern.prompt
        msgbox('No open File!','Split','error')
    end
    return
end
if leda2.data.events.N < 1
    if leda2.intern.prompt
        msgbox('File has no Events!','Split','error')
    end
    return
end
if isempty(leda2.analysis) || ~strcmp(leda2.analysis.method,'sdeco')
    if leda2.intern.prompt
        msgbox('This function requires that data have been analyzed with Continuous Decomposition Analysis!','Split','error')
    end
    return;
end

conditionnames = unique({leda2.data.events.event.name});
nrconds = numel(conditionnames);
nrboxlines = ceil(numel(conditionnames) / 3);   % reserve space for three columns of codition names
nrlines = nrboxlines + 4;
dy = 1 / (nrlines + 6);     % height of one line (consider some more space between lines => +6)

leda2.gui.split = [];

leda2.gui.split.fig = figure('Units','normalized','Position', ...
    [.25 .125 .5 .5], ...
    'Name','Plot event-related data (Split and Average)','MenuBar','none','NumberTitle','off');

leda2.gui.split.text_WindowLimits = uicontrol('Style','text','Units','normalized', ...
    'Position',[.05 1-2*dy .30 dy],'BackgroundColor',get(gcf,'Color'), ...
    'String','Window around events (from BEFORE to AFTER event) [sec]:', ...
    'HorizontalAlignment','left');

% START value
leda2.gui.split.edit_WindowStart = uicontrol('Style','edit','Units','normalized', ...
    'Position',[.4 1-2*dy .1 dy],'BackgroundColor',[1 1 1], ...
    'String',num2str(leda2.set.split.start,'%1.2f'));
% END value
leda2.gui.split.edit_WindowEnd   = uicontrol('Style','edit','Units','normalized', ...
    'Position',[.55 1-2*dy .1 dy],'BackgroundColor',[1 1 1], ...
    'String',num2str(leda2.set.split.end,'%1.2f'));

% VARIABLE
leda2.gui.split.text_SplitVariable = uicontrol('Style','text','Units','normalized', ...
    'Position',[.05 1-4*dy .30 dy],'BackgroundColor',get(gcf,'Color'), ...
    'String','Variable:','HorizontalAlignment','right');
% make CURRENT DATA available!!!!
leda2.gui.split.edit_SplitVariable = uicontrol('Style','popupmenu','Units','normalized', ...
    'Position',[.4 1-4*dy .25 dy], ...
    'String',leda2.set.split.variables, ...
    'Value',leda2.set.split.var);

% standard error
leda2.gui.split.stderr = uicontrol('Style','checkbox', 'Units','normalized', ...
    'Position', [.70 1-4*dy .25 dy], 'Min', 0, 'Max', 1, ...
    'String', 'Include standard error', 'HorizontalAlignment','Left','BackgroundColor',get(gcf,'Color'));

leda2.gui.split.text_Conditions = uicontrol('Style','text','Units','normalized', ...
    'Position',[.05 1-6.5*dy .15 dy],'BackgroundColor',get(gcf,'Color'), ...
    'String','Conditions:','HorizontalAlignment','left');

leda2.gui.split.all_Conditions = uicontrol('Style','checkbox', 'Units','normalized', ...
    'Position', [.25 1-6.5*dy .25 dy], 'Min', 0, 'Max', 1, ...
    'String', 'select all', 'HorizontalAlignment','Left','BackgroundColor',get(gcf,'Color'), 'Callback', 'leda_split(''select_all'')');

% create 'checkboxes' to enable selecting a subset of all conditions
for i = 1:(nrboxlines-1)
    for j = 1:3
        leda2.gui.split.edit_Conditions((i-1)*3+j) = ...
            uicontrol('Style','checkbox', 'Units','normalized', 'Position', ...
            [.03+.333*(j-1) (1-((7+i)*dy)) .333 dy], ... [6+condnamelength*(j-1) 10+(i*2) condnamelength 36], ...
            'Min', 0, 'Max', 1, ...
            'String',conditionnames{(i-1)*3+j}, 'HorizontalAlignment','Left','BackgroundColor',get(gcf,'Color'), ...
            'Callback', 'leda_split(''deselect_select_all'')');
    end
end
for j = 1:(nrconds-(nrboxlines-1)*3)
    leda2.gui.split.edit_Conditions((nrboxlines-1)*3+j) = ...
        uicontrol('Style','checkbox', 'Units','normalized', 'Position', ...
        [.03+.333*(j-1) (1-((7+nrboxlines)*dy)) .333 dy], ... [6+condnamelength*(j-1) 10+(i*2) condnamelength 36], ...
        'Min', 0, 'Max', 1, ...
        'String',conditionnames{(nrboxlines-1)*3+j}, 'HorizontalAlignment','Left','BackgroundColor',get(gcf,'Color'), ...
        'Callback', 'leda_split(''deselect_select_all'')');
end

% button
leda2.gui.split.butt_split = uicontrol('Units','normalized','Position', [.333 dy .333 dy], ...
    'String','Plot event-related data','Callback','leda_split(''take_settings'')');

end


function select_all
global leda2
if get(leda2.gui.split.all_Conditions, 'Value') == 0
    for i = 1:numel(leda2.gui.split.edit_Conditions)
        set(leda2.gui.split.edit_Conditions(i), 'Value', 0);
    end
else
    for i = 1:numel(leda2.gui.split.edit_Conditions)
        set(leda2.gui.split.edit_Conditions(i), 'Value', 1);
    end
end
end

function deselect_select_all
global leda2
if get(leda2.gui.split.all_Conditions, 'Value') == 1
    set(leda2.gui.split.all_Conditions, 'Value', 0);
else
    if all(cell2mat(get(leda2.gui.split.edit_Conditions, 'Value')))
        set(leda2.gui.split.all_Conditions, 'Value', 1);
    end
end
end

function take_settings
global leda2

leda2.set.split.start = str2double(get(leda2.gui.split.edit_WindowStart, 'String'));
leda2.set.split.end = str2double(get(leda2.gui.split.edit_WindowEnd, 'String'));

% Get variable
leda2.set.split.variable = leda2.set.split.variables{...
    get(leda2.gui.split.edit_SplitVariable, 'Value')};

% Get selected conditions
if get(leda2.gui.split.all_Conditions, 'Value') == 0
    conditionnames = unique({leda2.data.events.event.name});
    leda2.set.split.selectedconditions = ...
        conditionnames(logical(cell2mat(get(leda2.gui.split.edit_Conditions, 'Value'))));
    if isempty(leda2.set.split.selectedconditions)
        mbox = msgbox('No condition selected!', 'Condition error');
        return
    end
else
    leda2.set.split.selectedconditions = unique({leda2.data.events.event.name});
end

% get whether standard error should be plotted
leda2.set.split.stderr = get(leda2.gui.split.stderr, 'Value');

% check settings
if leda2.set.split.start > leda2.set.split.end
    mbox = msgbox('BEFORE event value is greater than AFTER event value!','Event window error');
    return
end
if leda2.set.split.start == leda2.set.split.end
    mbox = msgbox('Interval of 0 sec not allowed!','Event window error');
    return
end

split;

close(leda2.gui.split.fig)

end


function split
global leda2
%%
% name of variable to be processed
variable        = leda2.set.split.variable;
target_field    = ['split_' variable ];
% take all the conditions in a cell...
[leda2.set.split.selectedconditions i]   = unique({leda2.data.events.event.name},'first');
% get names of all events
alleventnames = {leda2.data.events.event.name};
condnames     = leda2.set.split.selectedconditions;
% number of conditions
nrcond                          = numel(condnames);
leda2.set.split.condition       = [leda2.data.events.event(i).nid];

% number of instances per condition
npcond = zeros(1,nrcond);
for i = 1:nrcond
    npcond(i) = sum(strcmp(leda2.set.split.selectedconditions{i}, alleventnames));
end
% convert the event's times (in sec) to sample nr. ...
etsamplenr = round([leda2.data.events.event.time] *leda2.data.samplingrate + 1);

% get FROM and TO and convert to sample number.
% NOTE: START and END are presumably seconds.
%       And START is usually negative.
from = leda2.set.split.start * leda2.data.samplingrate;
to   = leda2.set.split.end * leda2.data.samplingrate;
% lock them to sampling rate, but move both in same direction!
if from < to
    if round(from) < from
        to = floor(to);
    elseif round(from) > from
        to = ceil(to);
    else    % round(floor) == floor
        to = round(to);
    end
    from = round(from);
end
%% split data of the requested variable
i       = 0;
count   = zeros(1,nrcond);
for co = 1:nrcond%across conditions    
    for j = find(strcmp(condnames(co), alleventnames))%across trials
        i = i + 1;
        %trials
        leda2.analysis.(target_field).y(:,i)     = leda2.analysis.(variable)((etsamplenr(j) + from):(etsamplenr(j) + to));   % FROM is negative! 
        %condition indices
        leda2.analysis.(target_field).c(:,i)     = co;
        %
        count(co)                                = count(co)  + 1;
        leda2.analysis.(target_field).rank(i)    = count(co);
    end
    %compute the statistics and store them in an aligned way
    tc                                       = leda2.analysis.(target_field).y(:,leda2.analysis.(target_field).c == co);
    leda2.analysis.(target_field).mean(:,co) = mean(  tc,2);
    leda2.analysis.(target_field).med(:,co)  = median(tc,2);
    leda2.analysis.(target_field).std(:,co)  = std(   tc,1,2);
    leda2.analysis.(target_field).sem(:,co)  = leda2.analysis.(target_field).std(:,co)./sqrt(count(co));    
    leda2.analysis.(target_field).n(:,co)    = count(co);
end
leda2.analysis.(target_field).cond          = leda2.set.split.condition;
%time axis
x                                           = (from:to)/leda2.data.samplingrate;
leda2.analysis.(target_field).x             = repmat(x(:),1,size(leda2.analysis.(target_field).c,2));
%
%%   
% plot mean and indicate the event for all conditions in one plot
if ~leda2.intern.batchmode || (leda2.intern.batchmode && leda2.set.split.plot)
    figure;
    hold on
    leg = cell(1,nrcond);   % for the legend
    yvals = nan(nrcond, numel(leda2.analysis.(['split_' variable]).mean(:,1)));
%     if leda2.set.split.stderr
%         stderrvals = nan(nrcond, numel(leda2.analysis.(['split_' variable]).sem(:,1))   );
%     end
    for c = 1:nrcond        % gather data to plot in one variable
        leg{c} = sprintf('%s (n = %d)', condnames{c}, npcond(aligned_ca(c)));
        yvals(c,:) = leda2.analysis.(['split_' variable]).mean(:,c);
%         stderrvals(c,:) = [leda2.data.split(c).mean + leda2.data.split(c).stderr, ...
%             leda2.data.split(c).mean(end:-1:1) - leda2.data.split(c).stderr(end:-1:1)];
    end
    time = x;
    plt = plot(time,yvals);    % plot it
%     if leda2.set.split.stderr
%         col = get(plt, 'color');
%         for c = 1:nrcond
%             fill([time, time(end:-1:1)], stderrvals(c,:), ...
%                 col{c}, 'facealpha', .25, 'edgealpha', .25);
%         end
%     end
    plot([0; 0], get(gca, 'YLim')', 'r--');     % time of the event
    xlabel('time [sec]');
    ylabel(variable);
    legend(leg, 'location', 'best');
    % %
    legend boxoff
    [dummy subject_file] = fileparts(leda2.file.filename);
    supertitle(subject_file,1,'interpreter','none');    
    % %
    xlim([leda2.set.split.start, leda2.set.split.end])
    hold off
end

end