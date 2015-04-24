function import_data(datatype, pathname, filename)

global leda2
leda2.current.fileopen_ok = 0;

switch datatype
    case 'biotrace', ext = {'*.txt'};
    case 'biopac', ext = {'*.acq'};
    case 'biopacmat', ext = {'*.mat'};
    case 'cassylab', ext = {'*.lab'};
    case 'varioport', ext = {'*.vpd'};
    case 'visionanalyzer', ext = {'*.mat'};
    case 'vitaport', ext = {'*.asc'};
    case 'portilab', ext = {'*.txt'};
    case 'psychlab', ext = {'*.txt'};
    case 'userdef', ext = {'*.txt'};
    case 'mat', ext = {'*.mat'};
    case 'text', ext = {'*.txt'};
    case 'text2', ext = {'*.txt'};
    case 'text3', ext = {'*.txt'};

    otherwise
        if leda2.intern.prompt
            msgbox('Unknown filetype.','Info')
        end
        return
end

if nargin < 3
    [filename, pathname] = uigetfile(ext, ['Choose a ',datatype,' data-file']);
    if all(filename == 0) || all(pathname == 0) %Cancel
        return
    end
end
file = fullfile(pathname, filename);

%Try to import the selected data-file
try
    switch datatype
        case 'mat',
            load(file);
            conductance = data.conductance;
            time = data.time;
            event = data.event;

        case 'text'
            [time, conductance, event] = gettextdata(file);
            
        case 'text2'
            [time, conductance, event] = gettext2data(file);

        case 'text3'
            [time, conductance, event] = gettext3data(file);
            
        case 'biotrace'
            [time, conductance, event] = getBiotraceData(file);
            
        case 'biopac'
            [time, conductance, event] = getBiopacData(file);
            
        case 'biopacmat'
            [time, conductance, event] = getBiopacMatData(file);

        case 'cassylab',
            [time, conductance, event] = getcassydata(file);

        case 'varioport'       
            [time, conductance, event] = getVarioportData(file);

        case 'visionanalyzer'
            [time, conductance, event] = getVisionanalyzerData(file);
            
        case 'vitaport'
            [time, conductance, event] = getVitaportData(file);

        case 'portilab'
            [time, conductance, event] = getPortilabData(file);
            
        case 'psychlab'
            [time, conductance, event] = getPsychlabData(file);
                        
        case 'userdef'
            [time, conductance, event] = getuserdefdata(file);

    end
    
catch
   add2log(0,['Unable to import ',file,'.'],1,1,0,1,0,1)
   return
end

if isempty(conductance)
    return;
end


time = time(:)'; %force data in row
conductance = conductance(:)';

timeoffset = time(1);
time = time - timeoffset;


close_ledafile;
if leda2.file.open, return; end %closing failed

leda2.data = data;
%Load data
leda2.data.conductance.data = conductance;
leda2.data.time.data = time;
leda2.data.time.timeoff = timeoffset;
conductanceerror = sqrt(mean(diff(conductance).^2)/2);

%Get events
if ~isempty(event)
    leda2.data.events.event = event;  %event must contain at least event.time
    leda2.data.events.N = length(event);
    event_fields = fieldnames(event);

    %set dummy values for missing event-info
    for ev = 1:leda2.data.events.N

        leda2.data.events.event(ev).time = leda2.data.events.event(ev).time - timeoffset;
        if ~any(strcmp(event_fields, 'name'))
            leda2.data.events.event(ev).name = '';
        end
        if ~any(strcmp(event_fields, 'nid'))
            leda2.data.events.event(ev).nid = [];
        end
        if ~any(strcmp(event_fields, 'userdata'))
            leda2.data.events.event(ev).userdata = [];
        end
    end

end

leda2.file.filename = filename;
leda2.file.pathname = pathname;
leda2.file.date = clock;
leda2.file.log = {};
leda2.file.version = leda2.intern.version;
leda2.intern.current_dir = leda2.file.pathname;
cd(pathname);
leda2.file.open = 1;
file_changed(1);
add2log(1,[' Imported ',datatype,'-file ',file,' successfully.'],1,1,1);

refresh_data(0);    %Data statistics
update = 0;
%Positive SC values?
if leda2.data.conductance.min < 0 && ~leda2.intern.batchmode
    cmd = questdlg('Data shows negative values. This will complicate the analysis. Do you wish to correct this issue by adding a constant value?','Warning','Yes','No','Yes');
    if strcmp(cmd, 'Yes')
        update = 1;
        leda2.data.conductance.data = leda2.data.conductance.data - leda2.data.conductance.min + 1;
    end
end

%Downsample?
if (leda2.data.samplingrate > 32 || leda2.data.N > 36000) && ~leda2.intern.batchmode
    cmd = questdlg('Data is quite large. Do you wish to downsample your data in order to speed up the analysis?','Warning','Yes','No','Yes');
    if strcmp(cmd, 'Yes')
        update = 1;
        downsample;
    end
end
if update
    refresh_data(0);    %Data statistics
end

leda2.current.fileopen_ok = 1;
if leda2.intern.batchmode
    return;
end

plot_data; 