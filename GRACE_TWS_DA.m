function varargout = GRACE_TWS_DA(varargin)
% GRACE_TWS_DA MATLAB code for GRACE_TWS_DA.fig
%      GRACE_TWS_DA, by itself, creates a new GRACE_TWS_DA or raises the existing
%      singleton*.
%
%      H = GRACE_TWS_DA returns the handle to a new GRACE_TWS_DA or the handle to
%      the existing singleton*.
%
%      GRACE_TWS_DA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GRACE_TWS_DA.M with the given input arguments.
%
%      GRACE_TWS_DA('Property','Value',...) creates a new GRACE_TWS_DA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GRACE_TWS_DA_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GRACE_TWS_DA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GRACE_TWS_DA

% Last Modified by GUIDE v2.5 28-Aug-2020 09:16:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GRACE_TWS_DA_OpeningFcn, ...
                   'gui_OutputFcn',  @GRACE_TWS_DA_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GRACE_TWS_DA is made visible.
function GRACE_TWS_DA_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GRACE_TWS_DA (see VARARGIN)

% Choose default command line output for GRACE_TWS_DA
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Add the path to PACE: Principal Analysis by Conditional Expectation
addpath(genpath('C:\TWSModel\PACEV217\release2.17'));

% UIWAIT makes GRACE_TWS_DA wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GRACE_TWS_DA_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pbBrowseOutputFile.
function pbBrowseOutputFile_Callback(~, ~, handles)
% hObject    handle to pbBrowseOutputFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {  '*.hdr','Specify output data ENVI header file (*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify output data ENVI header file');
 if  filename > 0  
    set(handles.txtOutputFile, 'String', [pathname filename]);
 end 


function txtOutputFile_Callback(~, eventdata, handles)
% hObject    handle to txtOutputFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtOutputFile as text
%        str2double(get(hObject,'String')) returns contents of txtOutputFile as a double


% --- Executes during object creation, after setting all properties.
function txtOutputFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtOutputFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbClose.
function pbClose_Callback(hObject, ~, handles)
% hObject    handle to pbClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close;

% --- Executes on button press in pbComputeGRACETWSAandMasconIds.
function pbComputeGRACETWSAandMasconIds_Callback(hObject, eventdata, handles)
% hObject    handle to pbComputeGRACETWSAandMasconIds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% record the start time
tic;

%% get GRACE TWS input data file name
gtws_file = get(handles.txtGraceTwsFile, 'String');
if isempty(gtws_file)
    msgbox('Specify the GRACE TWS Anomaly NetCDF file.');
    return;
end

%% get the output data ENVI header file
hdr_file = get(handles.txtOutputFile, 'String');
out_file = get(handles.txtOutputFile, 'String');
if isempty(out_file)
    msgbox('Specify the output data ENVI header file.');
    return;
end
%% retrive the output data filder
[output_path, outfname, ext] = fileparts(out_file);

%% get EALCO TWS anomaly data file list
daily_files = get(handles.popFileType, 'Value');
if daily_files == 1
    fname = 'daily';
else
    fname = 'monthly';
end
file_list = get(handles.lstTerresTwsTimeSeries,'String');
if iscell(file_list)
    nf = length(file_list);
else
    nf = 1;
end
if nf<1
    msgbox('Add at least one LSM TWS data file to the LSM TWS Time Series list!');
    return;
end
%% create a time and time string list for the EALCO TWS input file list
time_list = zeros(nf,1);
time_str_list = cell(nf,1);
for i = 1:nf
    if ischar(file_list)
        ttwsi_file = file_list;
    else
        ttwsi_file = char(file_list{i});
    end
    if daily_files == 1
        [time_list(i), time_str_list{i}] = daysFromFileName(ttwsi_file,1);
    else
        [time_list(i), time_str_list{i}] = daysFromFileName(ttwsi_file,3);
    end
end
%% note: the time stamps is in days from Jan. 1, 2002

%% get the meta data parameters and processing options
switch get(handles.popDataType, 'Value')
    case 1      % daily data in int16
        data_type =  'int16';
    case 2      % monthly data in float32
        data_type =  'float32';
end
switch get(handles.popNoDataValue, 'Value')
    case 1      
        no_data =  -32760;
    case 2      
        no_data =  9999;
end

row = str2double(get(handles.txtRow, 'String'));    % row number of the GRACE TWS data
col = str2double(get(handles.txtCol, 'String'));    % col number of the GRACE TWS data
band = str2double(get(handles.txtBand, 'String'));  % the number of months 

%% add geo reference coordinates to EALCO TWS
out_lcc_latlon_file = [output_path '\Canada_Lat_Lon_lcc_grid.mat'];
if ~exist(out_lcc_latlon_file,'file')
    lamb = Lambert(...          % Define Canada LCC projection parameters
        [49.0, 77.0], ...       % Latitude of the first and second standard parallel
        0.0, -95, ...           % Central latitude and longitude
        0.0, 0.0, ...           % False easting and northing
        6378137.0, 298.257222100882711243 ...     % Ellipsoid: Semi-major axis and inverse flattening
    );

    lcc_ulx_e = -2600000.0;     
    lcc_uly_n = 10500000.0;
    % lcc_lrx_e = 3100000.0;      
    % lcc_lry_n = 5700000.0;
    lcc_dx = 5000.0;            
    lcc_dy = 5000.0;
    lcc_x = zeros(row, col);   
    lcc_y = zeros(row, col);
    lcc_lat = zeros(row, col);   
    lcc_lon = zeros(row, col);

    for i = 1:row
        for j =1:col
            lcc_x(i,j) = lcc_ulx_e + j*lcc_dx - lcc_dx/2;
            lcc_y(i,j) = lcc_uly_n - i*lcc_dy + lcc_dy/2;
            [lcc_lat(i,j), lcc_lon(i,j)] = lamb.cartesian2geographic(lcc_x(i,j), lcc_y(i,j));
        end
    end

    save(out_lcc_latlon_file,'lcc_x', 'lcc_y', 'lcc_lat', 'lcc_lon');
else
    load(out_lcc_latlon_file);
end

%% determine the geolocation area of CANADA
idx = find(lcc_lon < 0);
lcc_lon(idx) = lcc_lon(idx)+360;
minLat = min(lcc_lat(:));   maxLat = max(lcc_lat(:));
minLon = min(lcc_lon(:));   maxLon = max(lcc_lon(:));

%% read the GRACE TWS anomaly data from NETCDF file
% ncdisp(gtws_file)
gtws = ncread(gtws_file,'lwe_thickness');
gtws_std = ncread(gtws_file,'uncertainty');
% lat = ncread(gtws_file,'lat');
% lon = ncread(gtws_file,'lon');
day = ncread(gtws_file,'time');
lat_bounds = ncread(gtws_file,'lat_bounds');
lon_bounds = ncread(gtws_file,'lon_bounds');
day_bounds = ncread(gtws_file,'time_bounds');

%% determine the GRACE TWS anomaly data index range for CANADA area
[idx_lon,idx_lat,b] = size(gtws);
for i = 1:idx_lat
    if minLat >= lat_bounds(1,i) && minLat < lat_bounds(2,i)
        idx_lat_start = i;
        break;
    end
end
for i = 1:idx_lat
    if maxLat >= lat_bounds(1,i) && maxLat < lat_bounds(2,i)
        idx_lat_end = i;
        break;
    end
end
for i = 1:idx_lon
    if minLon >= lon_bounds(1,i) && minLon < lon_bounds(2,i)
        idx_lon_start = i;
        break;
    end
end
for i = 1:idx_lon
    if maxLon >= lon_bounds(1,i) && maxLon < lon_bounds(2,i)
        idx_lon_end = i;
        break;
    end
end
%% retreive GRACE TWS Anomaly and their uncerntainty data for Canada area
gtws = gtws(idx_lon_start:idx_lon_end,idx_lat_start:idx_lat_end,:);
gtws_std = gtws_std(idx_lon_start:idx_lon_end,idx_lat_start:idx_lat_end,:);
% lat = lat(idx_lat_start:idx_lat_end);
% lon = lon(idx_lon_start:idx_lon_end);
lat_bounds = lat_bounds(:,idx_lat_start:idx_lat_end);
lon_bounds = lon_bounds(:,idx_lon_start:idx_lon_end);


%% open log file to record computation progresses
logfile = [output_path '\GRACE_TWSA_Process.log'];
if exist(logfile,'file')
    delete logfile;
end
flog = fopen(logfile,'wt');

logmsg(flog,handles,'Processing started.');

%% initialize the output data container for all monthly datasets
%% the output data is presented in LCC coordinate system
if band > nf
    band = nf;
end

masconid = zeros(row, col);

%% determine the size of GRACE TWS anomaly input data
%% note: the size of gtws is for Canada coverage only
[idx_lon,idx_lat,b] = size(gtws);  % b may be >= band

%% compute or read the baseline values 
out_grace_baseline_file = [output_path '\Canada_GRACE_TWS_anomaly_' num2str(band) '_month_baseline.mat'];
if ~exist(out_grace_baseline_file,'file')
    gbase = zeros(idx_lon,idx_lat);
    for i = 1:idx_lon
        for j = 1:idx_lat
            gbase(i,j) = mean(gtws(i,j,1:band));
        end
    end
    save(out_grace_baseline_file, 'gbase');
else
    load(out_grace_baseline_file);
end

%% process data band by band for each monscon
%% note: each band represnts a monthly dataset
for i = 1:band         
    %% determine the day ranges for the monthly dataset i
    day_bounds_i = day_bounds(:,i);
    day_start = day_bounds_i(1);
    day_end = day_bounds_i(2);
    dayi = day(i);

    logmsg(flog,handles,' ');
    logmsg(flog,handles,['Processing GRACE TWS anomaly for the day: ' num2str(dayi)]);
    logmsg(flog,handles,' ');    
    
    %% read or compute monthly averaged EALCO TWS anomaly data
    if daily_files == 1
        %% retreive the daily EALCO TWS files for the GRACE TWS monthly dataset   
        idx_days = find(time_list>=day_start & time_list<=day_end);       
        if isempty(idx_days)        
            logmsg(flog,handles,['No matched daily LSM TWS files found for the days from ' num2str(day_start) ' to ' num2str(day_end)]);
            continue;
        else
            idx_middle = fix(length(idx_days)/2+0.5);
            %% determine a time string for output file name
            time_str = char(mtime_str_list(idx_middle));
        end
    else
        % idx_month = find(time_list>day_start+1 & time_list<day_end-1);
        idx_month = find(time_list>dayi-5 & time_list<dayi+5);
        if isempty(idx_month)        
            logmsg(flog,handles,['No matched monthly LSM TWS files found for the days from ' num2str(day_start) ' to ' num2str(day_end)]);
            continue;
        else
            %% determine a time string for output file name
            time_str = char(time_str_list(idx_month));
        end
    end
   
    %% retreive GRACE TWS data for the same month
    gtwsi = gtws(:,:,i)';
    gtwsi_std = gtws_std(:,:,i)';
    %% deduct the base line value
    gtwsi = (gtwsi - gbase');  
    %% convert the unit from cm to mm
    gtwsi = gtwsi*10;
    gtwsi_std = gtwsi_std*10;
    
    %% initialize the output data 
    gtwsio = ones(row, col)*no_data;
    gtwsio_std = ones(row, col)*no_data;
            
    %% start a loop to downscale and assimilate GRACE TWS to EALCO TWS
    continue_process_mascon = 1; 
    idx_lat1 = 0; idx_lon1 = 0;
    next_lat = 1;
    while continue_process_mascon > 0 && idx_lat1 <= idx_lat-1
        %% look for a GRACE TWS mascon data and its Lat Lon bounds
        if next_lat == 1
            idx_lat0 = idx_lat1+1; idx_lon0 = idx_lon1 + 1;
            mgtws = gtwsi(idx_lat0, idx_lon0);
            mgstd = gtwsi_std(idx_lat0,idx_lon0);
            for j=idx_lat0:idx_lat
                if gtwsi(j,idx_lon0) ~= mgtws
                    idx_lat1 = j-1;
                    break;
                elseif j == idx_lat
                    idx_lat1 = j;
                end
            end
        else
            idx_lon0 = idx_lon1 + 1;
            mgtws = gtwsi(idx_lat0, idx_lon0);
            mgstd = gtwsi_std(idx_lat0,idx_lon0);
        end
        for j=idx_lon0:idx_lon
            if gtwsi(idx_lat0,j) ~= mgtws
                idx_lon1 = j-1;
                break;
            elseif j == idx_lon
                idx_lon1 = j;
            end
        end
        lat0 = lat_bounds(1,idx_lat0);
        lat1 = lat_bounds(2,idx_lat1);
        lon0 = lon_bounds(1,idx_lon0);
        lon1 = lon_bounds(2,idx_lon1);
        %% check all valuse within the mascon are same (for debug test only)
        mg = gtwsi(idx_lat0:idx_lat1,idx_lon0:idx_lon1);
        G = unique(mg);
        ng = length(G);
        logmsg(flog,handles,[num2str(ng) ' GRACE TWS Values found for the mascon id ' num2str(continue_process_mascon)]);
        countg = zeros(ng,1);
        for k=1:ng
            countg(k) = length(find(mgtws == G(k)));
            logmsg(flog,handles,['GRACE TWS Value: ' num2str(G(k)) ' n = ' num2str(countg(k))]);
        end 
        
        %% determine the mascon id        
        m_id = continue_process_mascon; 
        
        %% retreive EALCO TWS within the mascon
        midx = find((lcc_lat>=lat0 & lcc_lat<lat1) & (lcc_lon>=lon0 & lcc_lon<lon1));
        
        %% set the GRACE TWS value to match EALCO TWS
        gtwsio(midx) = mgtws;
        gtwsio_std(midx) = mgstd; 

        %% assign a valid mascon id 
        masconid(midx) = m_id;
        
        %% contune to next mascon
        continue_process_mascon = continue_process_mascon + 1;
        if idx_lon1 == idx_lon && idx_lat1 == idx_lat
            continue_process_mascon = 0;
        elseif idx_lon1 == idx_lon && idx_lat1 < idx_lat
            next_lat = 1;
            idx_lon1 = 0;
        else
            idx_lon0 = idx_lon1+1;
            next_lat = 0;
        end                                          
    end
    %% completed processing a time stamp of one month

        
    %% save the results to the output file   
    out_gtwsi_file = [output_path '\GRC_MTWSA_' time_str '.dat'];
    fid = fopen(out_gtwsi_file,'w');
        fwrite(fid,gtwsio','float32');
        fwrite(fid,gtwsio_std','float32');
    fclose(fid);
    copyfile(hdr_file,[output_path '\GRC_MTWSA_' time_str '.hdr']);
       
    %% save the mask id into file
    out_mask_file = [output_path '\Canada_JPL_mascon_id_lcc_grid.dat'];
    if ~exist(out_mask_file,'file')
        fid = fopen(out_mask_file,'w');
        fwrite(fid,masconid','float32');
        fclose(fid);
        copyfile(hdr_file,[output_path '\Canada_JPL_mascon_id_lcc_grid.hdr']);
    end
end

ttime = toc;

disp(['>>> Elapsed time: ', num2str(ttime),' seconds. <<<'])

msgbox('Completed GRACE TWSA process successfully!');

logmsg(flog,handles,['>>> Elapsed time: ', num2str(ttime),' seconds. <<<']); 
logmsg(flog,handles,'Completed GRACE TWSA process successfully!'); 

fclose(flog);




function txtGraceTwsStd_Callback(hObject, eventdata, handles)
% hObject    handle to txtGraceTwsStd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtGraceTwsStd as text
%        str2double(get(hObject,'String')) returns contents of txtGraceTwsStd as a double


% --- Executes during object creation, after setting all properties.
function txtGraceTwsStd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtGraceTwsStd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtTerresTwsStd_Callback(~, eventdata, handles)
% hObject    handle to txtTerresTwsStd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtTerresTwsStd as text
%        str2double(get(hObject,'String')) returns contents of txtTerresTwsStd as a double


% --- Executes during object creation, after setting all properties.
function txtTerresTwsStd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtTerresTwsStd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popUnitTypeGraceTwsStd.
function popUnitTypeGraceTwsStd_Callback(~, eventdata, handles)
% hObject    handle to popUnitTypeGraceTwsStd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popUnitTypeGraceTwsStd contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popUnitTypeGraceTwsStd


% --- Executes during object creation, after setting all properties.
function popUnitTypeGraceTwsStd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popUnitTypeGraceTwsStd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popUnitTypeTerresTwsStd.
function popUnitTypeTerresTwsStd_Callback(hObject, eventdata, handles)
% hObject    handle to popUnitTypeTerresTwsStd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popUnitTypeTerresTwsStd contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popUnitTypeTerresTwsStd


% --- Executes during object creation, after setting all properties.
function popUnitTypeTerresTwsStd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popUnitTypeTerresTwsStd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popSolutionMethod.
function popSolutionMethod_Callback(hObject, eventdata, handles)
% hObject    handle to popSolutionMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popSolutionMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popSolutionMethod


% --- Executes during object creation, after setting all properties.
function popSolutionMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popSolutionMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtGraceTwsFile_Callback(hObject, eventdata, handles)
% hObject    handle to txtGraceTwsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtGraceTwsFile as text
%        str2double(get(hObject,'String')) returns contents of txtGraceTwsFile as a double


% --- Executes during object creation, after setting all properties.
function txtGraceTwsFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtGraceTwsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowseGraceTwsFile.
function pbBrowseGraceTwsFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbBrowseGraceTwsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {  '*.nc','Specify GRACE TWS Anomaly NetCDF file (*.nc)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify GRACE TWS Anomaly NetCDF file');
 if  filename > 0  
    set(handles.txtGraceTwsFile, 'String', [pathname filename]);
 end 


function txtGraceTwsIdMaskFile_Callback(~, eventdata, handles)
% hObject    handle to txtGraceTwsIdMaskFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtGraceTwsIdMaskFile as text
%        str2double(get(hObject,'String')) returns contents of txtGraceTwsIdMaskFile as a double


% --- Executes during object creation, after setting all properties.
function txtGraceTwsIdMaskFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtGraceTwsIdMaskFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowseGraceTwsIDMaskFile.
function pbBrowseGraceTwsIDMaskFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbBrowseGraceTwsIDMaskFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {  '*JPL_mascon_id*.*','Specify JPL Mascon Id for Canada LCC file (*JPL_mascon_id*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify JPL Mascon Id for Canada LCC file');
 if  filename > 0  
    set(handles.txtGraceTwsIdMaskFile, 'String', [pathname filename]);
 end 

% --- Executes on button press in pbAddTerresTwsFiles.
function pbAddTerresTwsFiles_Callback(hObject, eventdata, handles)
% hObject    handle to pbAddTerresTwsFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% get the files already in the list
listed_files = get(handles.lstTerresTwsTimeSeries,'String');

%% get the new files
[filenames, folder] = uigetfile( ...
    {  '*tws*.dat','Input image files (*tws*.dat)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Multiselect','on',...
    'Select trrestrial model TWS data files');
%% add the new files to the list
if ~isempty(filenames)
     filenames = fullfile(folder, filenames);
     if  isempty(listed_files) 
         listed_files = filenames;
     else
         listed_files = {listed_files; str2cell(filenames)};
     end
     set(handles.lstTerresTwsTimeSeries,'String',listed_files,'Value',1); 
end 

% --- Executes on button press in pbRemoveTwsFiles.
function pbRemoveTwsFiles_Callback(hObject, eventdata, handles)
% hObject    handle to pbRemoveTwsFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index_selected = get(handles.lstTerresTwsTimeSeries,'Value');
file_list = get(handles.lstTerresTwsTimeSeries,'String');
if iscell(file_list)
    file_list(index_selected)=[];
else
    file_list = [];
end
set(handles.lstTerresTwsTimeSeries,'String',file_list,'Value',1);

% --- Executes on selection change in lstTerresTwsTimeSeries.
function lstTerresTwsTimeSeries_Callback(hObject, ~, handles)
% hObject    handle to lstTerresTwsTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lstTerresTwsTimeSeries contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lstTerresTwsTimeSeries


% --- Executes during object creation, after setting all properties.
function lstTerresTwsTimeSeries_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lstTerresTwsTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtRow_Callback(hObject, eventdata, handles)
% hObject    handle to txtRow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtRow as text
%        str2double(get(hObject,'String')) returns contents of txtRow as a double


% --- Executes during object creation, after setting all properties.
function txtRow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtRow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtCol_Callback(hObject, eventdata, handles)
% hObject    handle to txtCol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtCol as text
%        str2double(get(hObject,'String')) returns contents of txtCol as a double


% --- Executes during object creation, after setting all properties.
function txtCol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtCol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, ~)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtBand_Callback(hObject, eventdata, handles)
% hObject    handle to txtBand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtBand as text
%        str2double(get(hObject,'String')) returns contents of txtBand as a double


% --- Executes during object creation, after setting all properties.
function txtBand_CreateFcn(hObject, eventdata, ~)
% hObject    handle to txtBand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbAddGraceTwsFiles.
function pbAddGraceTwsFiles_Callback(hObject, eventdata, handles)
% hObject    handle to pbAddGraceTwsFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% get the files already in the list
listed_files = get(handles.lstGraceTwsTimeSeries,'String');

%% get the new files
[filenames, folder] = uigetfile( ...
    {  '*tws*.dat','Input GRACE TWS files (*tws*.dat)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Multiselect','on',...
    'Select GRACE TWS data files');
%% add the new files to the list
if ~isempty(filenames)
     filenames = fullfile(folder, filenames);
     if  isempty(listed_files) 
         listed_files = filenames;
     else
         listed_files = {listed_files; str2cell(filenames)};
     end
     set(handles.lstGraceTwsTimeSeries,'String',listed_files,'Value',1); 
end 

% --- Executes on button press in pbRemoveGraceTwsFiles.
function pbRemoveGraceTwsFiles_Callback(hObject, eventdata, handles)
% hObject    handle to pbRemoveGraceTwsFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index_selected = get(handles.lstGraceTwsTimeSeries,'Value');
file_list = get(handles.lstGraceTwsTimeSeries,'String');
if iscell(file_list)
    file_list(index_selected)=[];
else
    file_list = [];
end
set(handles.lstGraceTwsTimeSeries,'String',file_list,'Value',1);

% --- Executes on selection change in lstGraceTwsTimeSeries.
function lstGraceTwsTimeSeries_Callback(hObject, eventdata, handles)
% hObject    handle to lstGraceTwsTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lstGraceTwsTimeSeries contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lstGraceTwsTimeSeries


% --- Executes during object creation, after setting all properties.
function lstGraceTwsTimeSeries_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lstGraceTwsTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popFileType.
function popFileType_Callback(hObject, eventdata, handles)
% hObject    handle to popFileType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popFileType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popFileType
% if get(handles.popFileType, 'Value')==1
%     set(handles.popDataType, 'Value',1);
%     set(handles.popNoDataValue, 'Value',1);
% else
%     set(handles.popDataType, 'Value',2);
%     set(handles.popNoDataValue, 'Value',2);
% end

% --- Executes during object creation, after setting all properties.
function popFileType_CreateFcn(hObject, eventdata, ~)
% hObject    handle to popFileType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popDataType.
function popDataType_Callback(hObject, eventdata, handles)
% hObject    handle to popDataType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popDataType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popDataType


% --- Executes during object creation, after setting all properties.
function popDataType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popDataType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popNoDataValue.
function popNoDataValue_Callback(hObject, eventdata, handles)
% hObject    handle to popNoDataValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popNoDataValue contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popNoDataValue


% --- Executes during object creation, after setting all properties.
function popNoDataValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popNoDataValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popWeightMethod.
function popWeightMethod_Callback(hObject, eventdata, handles)
% hObject    handle to popWeightMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popWeightMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popWeightMethod


% --- Executes during object creation, after setting all properties.
function popWeightMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popWeightMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbSpatialDataAssimilationBySCVCM.
function pbSpatialDataAssimilationBySCVCM_Callback(hObject, eventdata, handles)
% hObject    handle to pbSpatialDataAssimilationBySCVCM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% record the start time
tic;

%% get GRACE TWS input data file name
gtws_file = get(handles.txtGraceTwsFile, 'String');
if isempty(gtws_file)
    msgbox('Specify the GRACE TWSA data file.');
    return;
end

%% get JPL mascon grid id data file name
mask_file = get(handles.txtGraceTwsIdMaskFile, 'String');
if isempty(mask_file)
    msgbox('Specify the mascon id mask data file.');
    return;
end

%% get the EALCO baseline data file
baseline_file = get(handles.txtEALCOBaselineDataFile, 'String');
if isempty(baseline_file)
    msgbox('Specify the EALCO Baseline data file.');
    return;
end

%% get the output header file in ENVI data format
hdr_file = get(handles.txtOutputFile, 'String'); 
out_file = get(handles.txtOutputFile, 'String');
if isempty(out_file)
    msgbox('Specify the output header file.');
    return;
end
[output_path, outfname, ext] = fileparts(out_file);
%% create two new output folders for the downscaled TWSA and GSWSA
sd_twsa_output_path = [output_path '\SD_TWSA'];
sd_gswsa_output_path = [output_path '\SD_GSWSA'];
if exist(sd_twsa_output_path, 'dir') ~= 7
    mkdir(sd_twsa_output_path);
end
if exist(sd_gswsa_output_path, 'dir') ~= 7
    mkdir(sd_gswsa_output_path);
end

%% get EALCO model TWS data file list
daily_files = get(handles.popFileType, 'Value');
file_list = get(handles.lstTerresTwsTimeSeries,'String');
if iscell(file_list)
    nf = length(file_list);
else
    nf = 1;
end
if nf<1
    msgbox('Add at least one EALCO TWS data file to the EALCO TWS Time Series list!');
    return;
end

%% create a time and time string list for the EALCO TWS input file list
time_list = zeros(nf,1);
time_str_list = cell(nf,1);
for i = 1:nf
    if ischar(file_list)
        ttwsi_file = file_list;
    else
        ttwsi_file = char(file_list{i});
    end
    if daily_files == 1
        [time_list(i), time_str_list{i}] = daysFromFileName(ttwsi_file,1);
    else
        [time_list(i), time_str_list{i}] = daysFromFileName(ttwsi_file,3);
    end
end
%% note: the time stamps is in days from Jan. 1, 2002

%% get GRACE TWS anomaly data file list
gfile_list = get(handles.lstGraceTwsTimeSeries,'String');
if iscell(gfile_list)
    ngf = length(gfile_list);
else
    ngf = 1;
end
if ngf<1
    msgbox('Add at least one GRACE TWS data file to the list!');
    return;
end
%% create time and time string list for GRACE TWS anomaly files
gtime_list = zeros(ngf,1);
gtime_str_list = cell(ngf,1);
for i = 1:ngf
    if ischar(gfile_list)
        gtwsi_file = gfile_list;
    else
        gtwsi_file = char(gfile_list{i});
    end
    [gtime_list(i), gtime_str_list{i}] = daysFromFileName(gtwsi_file,3);
end

%% get the meta data parameters and processing options
switch get(handles.popDataType, 'Value')
    case 1      % daily data in int16
        data_type =  'int16';
    case 2      % monthly data in float32
        data_type =  'float32';
end
switch get(handles.popNoDataValue, 'Value')
    case 1      
        no_data =  -32760;
    case 2      
        no_data =  9999;
end

weight = get(handles.popWeightMethod, 'Value');

row = str2double(get(handles.txtRow, 'String'));    % row number of the GRACE TWS data
col = str2double(get(handles.txtCol, 'String'));    % col number of the GRACE TWS data
band = str2double(get(handles.txtBand, 'String'));  % the number of months 

%% get the precalculated scale factors file name for Canada LCC grid
apply_gain_factors = get(handles.chkApplyGainFactors, 'Value');
scale_file = get(handles.txtGainFactorFile, 'String');

%% read in the pixel id mask input data
fid = fopen(mask_file);
    mask = (fread(fid,[col, row],'float32')');
fclose(fid);

%% read in EALCO baseline data
fid = fopen(baseline_file);
    baseline = (fread(fid,[col, row],'float32')');
fclose(fid);

%% read the GRACE TWS data from NETCDF file
%ncdisp(gtws_file)
% gtws = ncread(gtws_file,'lwe_thickness');
% gtws_std = ncread(gtws_file,'uncertainty');
% lat = ncread(gtws_file,'lat');
% lon = ncread(gtws_file,'lon');
% lat_bounds = ncread(gtws_file,'lat_bounds');
% lon_bounds = ncread(gtws_file,'lon_bounds');
day = ncread(gtws_file,'time');
day_bounds = ncread(gtws_file,'time_bounds');

%% open log file to record computation progresses
logfile = [output_path '\SD_Processing.log'];
if exist(logfile,'file')
    delete logfile;
end
flog = fopen(logfile,'wt');

%% open a statistic log file to record computation results
st_logfile = [output_path '\SD_Statistics.log'];
if exist(st_logfile,'file')
    delete st_logfile;
end
fstlog = fopen(st_logfile,'wt');

logmsg(flog,handles,'Processing started.');
logmsg(fstlog,handles,'Processing Statistics:');

%% initialize the output data container for all monthly datasets
%% the output data is presented in LCC coordinate system
all_atws = ones(row, col, band)*no_data;
all_ttws = ones(row, col, band)*no_data;
all_gtws = ones(row, col, band)*no_data;
all_stws = ones(row, col, band)*no_data;
all_gws = ones(row, col, band)*no_data;

if daily_files == 1 && apply_gain_factors == 1 && isempty(scale_file)
    all_scal = ones(row, col, band)*no_data;
end

%% read the gain factors if selected to consider the leakage errors
if apply_gain_factors == 1 && ~isempty(scale_file)
    fid = fopen(scale_file);
        scale = (fread(fid,[col, row], data_type)');
        scale_std = (fread(fid,[col, row], data_type)');
    fclose(fid);   
end

%% determe how many mascons are in the selected test area
maskvalues = unique(mask);
% maskvalues(maskvalues==0)=[];   % remove the zero 
% maskvalues = [234,235,236,271,272,305,306,307,339,340]';
% maskvalues = [273,274,308]';
% maskvalues = [162,163,203,204,205,206]';

maskvalues(maskvalues==0)=[];
nmask = length(maskvalues);

%% process data band by band and monscon by monscon
if band > nf 
    band = nf;
end
for i=1:band         %% note: each band represnts a monthly dataset
    %% determine the day ranges for the monthly dataset i
    day_bounds_i = day_bounds(:,i);
    day_start = day_bounds_i(1);
    day_end = day_bounds_i(2);
    dayi = day(i);
    
    idx_gmonth = find(gtime_list>dayi-5 & gtime_list<dayi+5);
    if isempty(idx_gmonth)
        continue;
    else %% read the monthly GRACE TWS from file
        %% read in monthly average GRACE TWS data directly from file
        if ischar(gfile_list)
            gtwsi_file = gfile_list;
        else
            gtwsi_file = char(gfile_list{idx_gmonth});
        end
        %% read in the GRACE TWS from the monthly file
        fid = fopen(gtwsi_file);
            gtwsi = (fread(fid,[col, row], data_type)');
            gstdi = (fread(fid,[col, row], data_type)');
        fclose(fid);
    end

    logmsg(flog,handles,' ');
    logmsg(flog,handles,['Processing GRACE TWS for the day: ' num2str(dayi)]);
    logmsg(flog,handles,' ');    
    
    %% read or compute monthly averaged EALCO TWS data
    if daily_files == 1
        %% retreive the daily EALCO TWS files for the GRACE TWS monthly dataset   
        idx_days = find(time_list>=day_start & time_list<=day_end);       
        if isempty(idx_days)        
            logmsg(flog,handles,['No matched daily EALCO TWS files found for the days from ' num2str(day_start) ' to ' num2str(day_end)]);
            continue;
        else   
            %% continue processing the matched GRACE and EALCO data    
            logmsg(flog,handles,[num2str(length(idx_days)) ' daily EALCO TWS files found for the days from ' num2str(day_start) ' to ' num2str(day_end)]);
            mfile_list = file_list(idx_days);
            %% compute monthly average of EALCO TWS from the input file list
            if apply_gain_factors == 1 && isempty(scale_file)
                [ttwsi,tstdi,dtwsi] = averageFromInputFiles(mfile_list, row, col, data_type, no_data);
            else
                [ttwsi,tstdi] = averageFromInputFiles(mfile_list, row, col, data_type, no_data);
            end
            idx_middle = fix(length(idx_days)/2+0.5);
            %% determine output file name according to the time
            time_str = char(time_str_list(idx_middle));
        end
    else
        % idx_month = find(time_list>day_start+1 & time_list<day_end-1);
        idx_month = find(time_list>dayi-5 & time_list<dayi+5);        
        if isempty(idx_month)        
            logmsg(flog,handles,['No matched monthly EALCO TWS files found for the days from ' num2str(day_start) ' to ' num2str(day_end)]);
            continue;
        else   
            %% continue processing the matched GRACE and EALCO data    
            logmsg(flog,handles,[num2str(length(idx_month)) ' monthly EALCO TWS file found for the days from ' num2str(day_start) ' to ' num2str(day_end)]);
            %% read in monthly average EALCO TWS data directly from file
            if ischar(file_list)
                ttwsi_file = file_list;
            else
                ttwsi_file = char(file_list{idx_month});
            end
            %% read in EALCO TWS data from the data file
            fid = fopen(ttwsi_file);
                ttwsi = (fread(fid,[col, row], data_type)');
                tstdi = (fread(fid,[col, row], data_type)');
            fclose(fid);
            %% determine output file name according to the time
            time_str = char(time_str_list(idx_month));
        end
    end
   
    fprintf(fstlog,'%s\n', ' ');
    fprintf(fstlog,'%s\n',['Processing EALCO TWS for the days from ' num2str(day_start) ' to ' num2str(day_end)]);
    fprintf(fstlog,'%s\n','   ID  #OfTTWS    MinTWSA     MaxTWSA    MeanTWSA      GTWSA     MeanSCAL   MeanSTWSA   MeanATWSA   MeanGSWSA      TTStd0      TTStd1      GTStd0      GTStd1       STStd      ATStd       GWStd        MinD        MaxD       MeanD    Iteration');

    %% initialize the output data after assimilation vy copy the ttws
    atwsio = ones(row, col)*no_data;
    gtwsio = ones(row, col)*no_data;
    ttwsio = ones(row, col)*no_data;
    stwsio = ones(row, col)*no_data;
    gwsio = ones(row, col)*no_data;

    atwsio_std = ones(row, col)*no_data;
    gwsio_std = ones(row, col)*no_data;
    
    if daily_files == 1 && apply_gain_factors == 1 && isempty(scale_file)
        scalio = ones(row, col)*no_data;
        scalio_std = ones(row, col)*no_data;
    end
        
    %% start a loop to downscale and assimilate GRACE TWS to EALCO TWS according to the mascon id
    for j = 1:nmask
        logmsg(flog,handles,['Processing mascon id: ' num2str(maskvalues(j))]);
        %% find grid cells with the mascon id maskvalues(j)
        midx = find(mask == maskvalues(j));
        %% retrive the corresponding GRACE and EALCO TWS values
        mgtws = gtwsi(midx);
        mgstd = gstdi(midx);
        mttws = ttwsi(midx);
        mtstd = tstdi(midx);
        
        basev = baseline(midx); % baseline values
        
%         if maskvalues(j) == 142
%             disp('check data!');
%         end

        %% check and remove possible invalid data points in EALCO TWS data
        %inval = find(mttws == 9999 | mttws == -32760);
        inval = find(abs(mttws) >= 9999); % cover invalid values <-32760 and 9999
        if ~isempty(inval)
            midx(inval) = [];
            mttws(inval) = [];
            mtstd(inval) = [];
            mgtws(inval) = [];
            mgstd(inval) = [];
            basev(inval) = [];
        end
        if isempty(midx)
            logmsg(flog,handles,['No valid EALCO TWS Values found for the mascon id ' num2str(maskvalues(j))]);
            continue;
        else
            logmsg(flog,handles,[num2str(length(midx)) ' valid EALCO TWS Values found for the mascon id ' num2str(maskvalues(j))]);
        end
        
        %% check unique GRACE TWS value for the mascon id
        G = unique(mgtws); G_std = unique(mgstd);
        ng = length(G); ng_std = length(G_std);
        logmsg(flog,handles,[num2str(ng) ' GRACE TWS Values found for the mascon id ' num2str(maskvalues(j))]);
        logmsg(flog,handles,[num2str(ng_std) ' GRACE TWS STD found for the mascon id ' num2str(maskvalues(j))]);
        if ng > 1 
            countg = zeros(ng,1);
            for k=1:ng
                countg(k) = length(find(G == G(k)));
                logmsg(flog,handles,['GRACE TWS Value: ' num2str(G(k)) ' n = ' num2str(countg(k))]);
            end
            [~, idx2] = max(countg);
            G = G(idx2);
            G_std= G_std(idx2);
            logmsg(flog,handles,['Apply GRACE TWS Value ' num2str(G) ' for the mascon id ' num2str(maskvalues(j))]);
            logmsg(flog,handles,['Apply GRACE TWS STD ' num2str(G_std) ' for the mascon id ' num2str(maskvalues(j))]);
        end
        
        %% calculate or retreive scale factors
        if daily_files == 1 && apply_gain_factors == 1 && isempty(scale_file)
            [mscal,msstd] = scaleFactorFromEALCO(dtwsi,baseline,midx);   % important: to include baseline values for scale factor estimation 
            scalio(midx) = mscal;
            scalio_std(midx) = msstd;
        elseif apply_gain_factors == 1 && ~isempty(scale_file)
            mscal = scale(midx);
            msstd = scale_std(midx);      
        end
        
        %% keep a copy of the original GRACE TWSA
        gtwsio(midx) = mgtws;
        
        %% calculate scaled GRACE TWS 
        %% to include the base value for the scaled GRACE TWS
        % mgtws = mscal*(G + mean(basev));        
        %% to exclude the base value for the assimilation
        % mgtws = mgtws - basev;
        
        mgtws = mscal*G;
        mgstd = mscal*G_std;
        
        %% call assimilation function
        [matws, mastd, mgsws, mgstd, mstws, tstd, gstd, iteration] = assimilateMasconBySCVCM(mgtws,mgstd,mttws,mtstd,weight);
        logmsg(flog,handles,['SCVCM iteration number: ' num2str(iteration) ]);
        
        %% put the assimilated TWS to the initialized output dataset
        atwsio(midx) = matws; atwsio_std(midx) = mastd;
        ttwsio(midx) = mttws;
        stwsio(midx) = mstws;
        gwsio(midx) = mgsws; gwsio_std(midx) = mgstd;

        %% check the differences between the TTWS and ATWS
        delta = matws - mttws - mgsws;       
        fprintf(fstlog,'%+5s', num2str(maskvalues(j),'% 5d'));
        fprintf(fstlog,'%+8s', num2str(length(mttws),'% 8d'));
        fprintf(fstlog,'%+12s', num2str(min(mttws),'% 12.4f'));
        fprintf(fstlog,'%+12s', num2str(max(mttws),'% 12.4f'));
        
        fprintf(fstlog,'%+12s', num2str(mean(mttws),'% 12.4f'));
        fprintf(fstlog,'%+12s', num2str(G,'% 12.4f'));
        fprintf(fstlog,'%+12s', num2str(mean(mscal),'% 12.4f'));        
        fprintf(fstlog,'%+12s', num2str(mean(mstws),'% 12.4f'));
        fprintf(fstlog,'%+12s', num2str(mean(matws),'% 12.4f'));
        fprintf(fstlog,'%+12s', num2str(mean(mgsws),'% 12.4f'));
        
        fprintf(fstlog,'%+12s', num2str(mean(mtstd),'% 12.4f'));
        fprintf(fstlog,'%+12s', num2str(tstd,'% 12.4f'));
        fprintf(fstlog,'%+12s', num2str(G_std,'% 12.4f'));
        fprintf(fstlog,'%+12s', num2str(gstd,'% 12.4f'));
        fprintf(fstlog,'%+12s', num2str(mean(msstd),'% 12.4f'));        
        fprintf(fstlog,'%+12s', num2str(mean(mastd),'% 12.4f'));
        fprintf(fstlog,'%+12s', num2str(mean(mgstd),'% 12.4f'));

        fprintf(fstlog,'%+12s', num2str(min(delta),'% 12.4f'));
        fprintf(fstlog,'%+12s', num2str(max(delta),'% 12.4f'));
        fprintf(fstlog,'%+12s', num2str(mean(delta),'% 12.4f'));

        fprintf(fstlog,'%+8s\n', num2str(iteration,'% 8d'));

    end
    %% completed processing one monthly dataset
        
    %% put the monthly data into a channel
    all_atws(:,:,i) = atwsio;
    all_gtws(:,:,i) = gtwsio;
    all_ttws(:,:,i) = ttwsio;
    all_gws(:,:,i) = gwsio;
    all_stws(:,:,i) = stwsio;
       
    if daily_files == 1 && apply_gain_factors == 1 && isempty(scale_file)
        all_scal(:,:,i) = scalio;
        % out_scali_file = [output_path '\E_SCALE_' time_str '.dat'];
        % fid = fopen(out_scali_file,'w');
        %     fwrite(fid,scalio','float32');
        % fclose(fid);
        % copyfile(hdr_file,[output_path '\E_SCALE_' time_str '.hdr']);
        % 
        % out_scali_std_file = [output_path '\E_SCALE_std_' time_str '.dat'];
        % fid = fopen(out_scali_std_file,'w');
        %     fwrite(fid,scalio_std','float32');
        % fclose(fid);
        % copyfile(hdr_file,[output_path '\E_SCALE_std_' time_str '.hdr']);
    end

    %% save the results to the output file
    out_atwsi_file = [sd_twsa_output_path '\SD_MTWSA_' time_str '.dat'];
    fid = fopen(out_atwsi_file,'w');
        fwrite(fid,atwsio','float32');
        fwrite(fid,atwsio_std','float32');
    fclose(fid);
    copyfile(hdr_file,[sd_twsa_output_path '\SD_MTWSA_' time_str '.hdr']);
    
    out_gwsi_file = [sd_gswsa_output_path '\SD_MGSWSA_' time_str '.dat'];
    fid = fopen(out_gwsi_file,'w');
        fwrite(fid,gwsio','float32');
        fwrite(fid,gwsio_std','float32');
    fclose(fid);
    copyfile(hdr_file,[sd_gswsa_output_path '\SD_MGSWSA_' time_str '.hdr']);   
end

%% Output all assimilated monthly results into one file
out_file = [output_path '\SD_All_TWSA.dat'];
fid = fopen(out_file,'w');
for k=1:band
    fwrite(fid,all_atws(:,:,k)','float32');
end
fclose(fid);

out_file = [output_path '\GR_All_TWSA.dat'];
fid = fopen(out_file,'w');
for k=1:band
    fwrite(fid,all_gtws(:,:,k)','float32');
end
fclose(fid);

out_file = [output_path '\EA_All_TWSA.dat'];
fid = fopen(out_file,'w');
for k=1:band
    fwrite(fid,all_ttws(:,:,k)','float32');
end
fclose(fid);


out_file = [output_path '\SD_All_GSWSA.dat'];
fid = fopen(out_file,'w');
for k=1:band
    fwrite(fid,all_gws(:,:,k)','float32');
end
fclose(fid);    
out_file = [output_path '\SD_All_STWSA.dat'];
fid = fopen(out_file,'w');
for k=1:band
    fwrite(fid,all_stws(:,:,k)','float32');
end
fclose(fid);

if daily_files == 1 && apply_gain_factors == 1 && isempty(scale_file)
    out_file = [output_path '\SD_ALL_SCALE.dat'];
    fid = fopen(out_file,'w');
    for k=1:band
        fwrite(fid,all_scal(:,:,k)','float32');
    end
    fclose(fid);
end

ttime = toc;

disp(['>>> Elapsed time: ', num2str(ttime),' seconds. <<<'])

msgbox('Completed data assimilation successfully!');

logmsg(flog,handles,['>>> Elapsed time: ', num2str(ttime),' seconds. <<<']); 
logmsg(flog,handles,'Completed data assimilation successfully!'); 

fclose(flog);
fclose(fstlog);


% --- Executes on button press in pbComputeEALCOGriddedGainScaleFactors.
function pbComputeEALCOGriddedGainScaleFactors_Callback(hObject, eventdata, handles)
% hObject    handle to pbComputeEALCOGriddedGainScaleFactors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% record the start time
tic;

%% get GRACE TWS input data file name
gtws_file = get(handles.txtGraceTwsFile, 'String');
if isempty(gtws_file)
    msgbox('Specify the GRACE TWS Anomaly NetCDF file.');
    return;
end

%% get the EALCO baseline data file
baseline_file = get(handles.txtEALCOBaselineDataFile, 'String');
if isempty(baseline_file)
    msgbox('Specify the EALCO Baseline data file.');
    return;
end

%% get the output data ENVI header file
hdr_file = get(handles.txtOutputFile, 'String');
out_file = get(handles.txtOutputFile, 'String');
if isempty(out_file)
    msgbox('Specify the output data ENVI header file.');
    return;
end
%% retrive the output data filder
[output_path, outfname, ext] = fileparts(out_file);

%% get EALCO TWS anomaly data file list
daily_files = get(handles.popFileType, 'Value');
file_list = get(handles.lstTerresTwsTimeSeries,'String');
if iscell(file_list)
    nf = length(file_list);
else
    nf = 1;
end
if nf<1
    msgbox('Add at least one LSM TWS data file to the LSM TWS Time Series list!');
    return;
end

%% get the meta data parameters and processing options
switch get(handles.popDataType, 'Value')
    case 1      % daily data in int16
        data_type =  'int16';
    case 2      % monthly data in float32
        data_type =  'float32';
end
switch get(handles.popNoDataValue, 'Value')
    case 1      
        no_data =  -32760;
    case 2      
        no_data =  9999;
end

row = str2double(get(handles.txtRow, 'String'));    % row number of the GRACE TWS data
col = str2double(get(handles.txtCol, 'String'));    % col number of the GRACE TWS data
band = str2double(get(handles.txtBand, 'String'));  % the number of months 

%% add geo reference coordinates to EALCO TWS
out_lcc_latlon_file = [output_path '\Canada_Lat_Lon_lcc_grid.mat'];
if ~exist(out_lcc_latlon_file,'file')
    lamb = Lambert(...          % Define Canada LCC projection parameters
        [49.0, 77.0], ...       % Latitude of the first and second standard parallel
        0.0, -95, ...           % Central latitude and longitude
        0.0, 0.0, ...           % False easting and northing
        6378137.0, 298.257222100882711243 ...     % Ellipsoid: Semi-major axis and inverse flattening
    );

    lcc_ulx_e = -2600000.0;     
    lcc_uly_n = 10500000.0;
    % lcc_lrx_e = 3100000.0;      
    % lcc_lry_n = 5700000.0;
    lcc_dx = 5000.0;            
    lcc_dy = 5000.0;
    lcc_x = zeros(row, col);   
    lcc_y = zeros(row, col);
    lcc_lat = zeros(row, col);   
    lcc_lon = zeros(row, col);

    for i = 1:row
        for j =1:col
            lcc_x(i,j) = lcc_ulx_e + j*lcc_dx - lcc_dx/2;
            lcc_y(i,j) = lcc_uly_n - i*lcc_dy + lcc_dy/2;
            [lcc_lat(i,j), lcc_lon(i,j)] = lamb.cartesian2geographic(lcc_x(i,j), lcc_y(i,j));
        end
    end

    save(out_lcc_latlon_file,'lcc_x', 'lcc_y', 'lcc_lat', 'lcc_lon');
else
    load(out_lcc_latlon_file);
end

%% determine the geolocation area of CANADA
idx = find(lcc_lon < 0);
lcc_lon(idx) = lcc_lon(idx) + 360;
minLat = min(lcc_lat(:));   maxLat = max(lcc_lat(:));
minLon = min(lcc_lon(:));   maxLon = max(lcc_lon(:));

%% read the GRACE TWS anomaly data from NETCDF file
% ncdisp(gtws_file)
gtws = ncread(gtws_file,'lwe_thickness');
% gtws_std = ncread(gtws_file,'uncertainty');
% lat = ncread(gtws_file,'lat');
% lon = ncread(gtws_file,'lon');
% day = ncread(gtws_file,'time');
lat_bounds = ncread(gtws_file,'lat_bounds');
lon_bounds = ncread(gtws_file,'lon_bounds');
% day_bounds = ncread(gtws_file,'time_bounds');

%% determine the GRACE TWS anomaly data index range for CANADA area
[idx_lon,idx_lat,b] = size(gtws);
for i = 1:idx_lat
    if minLat >= lat_bounds(1,i) && minLat < lat_bounds(2,i)
        idx_lat_start = i;
        break;
    end
end
for i = 1:idx_lat
    if maxLat >= lat_bounds(1,i) && maxLat < lat_bounds(2,i)
        idx_lat_end = i;
        break;
    end
end
for i = 1:idx_lon
    if minLon >= lon_bounds(1,i) && minLon < lon_bounds(2,i)
        idx_lon_start = i;
        break;
    end
end
for i = 1:idx_lon
    if maxLon >= lon_bounds(1,i) && maxLon < lon_bounds(2,i)
        idx_lon_end = i;
        break;
    end
end

%% retreive GRACE TWS Anomaly and lat lon boundaries for Canada area
gtws = gtws(idx_lon_start:idx_lon_end,idx_lat_start:idx_lat_end,:);
% gtws_std = gtws_std(idx_lon_start:idx_lon_end,idx_lat_start:idx_lat_end,:);
% lat = lat(idx_lat_start:idx_lat_end);
% lon = lon(idx_lon_start:idx_lon_end);
lat_bounds = lat_bounds(:,idx_lat_start:idx_lat_end);
lon_bounds = lon_bounds(:,idx_lon_start:idx_lon_end);

%% open log file to record computation progresses
logfile = [output_path '\LSM_SF_process.log'];
if exist(logfile,'file')
    delete logfile;
end
flog = fopen(logfile,'wt');

%% open a statistic log file to record some statistics of results
st_logfile = [output_path '\LSM_SF_Statistics.log'];
if exist(st_logfile,'file')
    delete st_logfile;
end
fstlog = fopen(st_logfile,'wt');

logmsg(flog,handles,'Processing started.');
logmsg(fstlog,handles,'Processing Statistics:');

%% initialize the output data containers 
%% the output data is presented in LCC coordinate system
dtws = ones(row, col, nf)*no_data;
scal = ones(row, col)*no_data;
sstd = ones(row, col)*no_data;

%% read the EALCO TWS time series data
for i = 1:nf
    dfile = char(file_list(i));
    fid = fopen(dfile);
        dtws(:,:,i) = (fread(fid,[col, row], data_type))';
    fclose(fid);
end

%% read in EALCO baseline data
fid = fopen(baseline_file);
    baseline = (fread(fid,[col, row],'float32')');
fclose(fid);

%% retreive one dataset for evaluation of data availability
gtwsi = gtws(:,:,1)'; 
ttwsi = dtws(:,:,1);    % used to evaluate availability of EALCO TWS

[idx_lon,idx_lat,b] = size(gtws);  % b may be >= band

fprintf(fstlog,'%s\n', ' ');
fprintf(fstlog,'%s\n','Statistics of the scale factors calculated from LSM TWS');
fprintf(fstlog,'%s\n','  ID                 Bounds           #OfCells #OfValid    MeanSF       MinSF       MaxSF     MeanSTD      MinSTD      MaxSTD');
       
%% start a loop to calculate the scale factors from EALCO TWS time series
continue_process_mascon = 1; 
idx_lat1 = 0; idx_lon1 = 0;
next_lat = 1;
while continue_process_mascon > 0 && idx_lat1 <= idx_lat-1
    %% look for Lat Lon boundary of a mascon
    if next_lat == 1
        idx_lat0 = idx_lat1+1; idx_lon0 = idx_lon1 + 1;
        mgtws = gtwsi(idx_lat0, idx_lon0);
        for j=idx_lat0:idx_lat
            if gtwsi(j,idx_lon0) ~= mgtws
                idx_lat1 = j-1;
                break;
            elseif j == idx_lat
                idx_lat1 = j;
            end
        end
    else
        idx_lon0 = idx_lon1 + 1;
        mgtws = gtwsi(idx_lat0, idx_lon0);
    end
    for j=idx_lon0:idx_lon
        if gtwsi(idx_lat0,j) ~= mgtws
            idx_lon1 = j-1;
            break;
        elseif j == idx_lon
            idx_lon1 = j;
        end
    end
    lat0 = lat_bounds(1,idx_lat0);
    lat1 = lat_bounds(2,idx_lat1);
    lon0 = lon_bounds(1,idx_lon0);
    lon1 = lon_bounds(2,idx_lon1);
    %% check all valuse within the mascon are same (for debug test only)
    % mg = gtwsi(idx_lat0:idx_lat1,idx_lon0:idx_lon1);
    % clear mg;  
    %% determine the mascon id        
    m_id = continue_process_mascon;
    %% determine the index of EALCO TWS data within the mascon
    midx = find((lcc_lat>=lat0 & lcc_lat<lat1) & (lcc_lon>=lon0 & lcc_lon<lon1));
    mttws = ttwsi(midx);
    total_cells = length(midx);
    %% remove invalid EALCO TWS
    % inval = find(abs(mttws) == -32760 | mttws == 9999 | isnan(mttws));
    inval = find(abs(mttws) >= 9999 | isnan(mttws));
    if ~isempty(inval)
        midx(inval) = [];
        mttws(inval) = [];
    end
    valid_cells = length(midx);
    if isempty(mttws)
        logmsg(flog,handles,' ');
        logmsg(flog,handles,['No valid LSM TWS found within the mascon ' num2str(continue_process_mascon) ': lat ' num2str(lat0) '-' num2str(lat1) ' lon ' num2str(lon0) '-' num2str(lon1)]);
        %% contune to next mascon
        continue_process_mascon = continue_process_mascon + 1;
    else            
        logmsg(flog,handles,' ');
        logmsg(flog,handles,['Found ' num2str(valid_cells) ' LSM TWS within the mascon ' num2str(continue_process_mascon) ': lat ' num2str(lat0) '-' num2str(lat1) ' lon ' num2str(lon0) '-' num2str(lon1)]);            
        %% calculate scale factors for a mascon
        [mscal,msstd] = scaleFactorFromEALCO(dtws,baseline,midx); % important: to include baseline value for scale factor estimation
        %% put the scale factors to the output dataset
        scal(midx) = mscal;
        sstd(midx) = msstd;

        fprintf(fstlog,'%+5s', num2str(m_id,'% 6d'));
        fprintf(fstlog,'%-32s', ['  lat ' num2str(lat0) '-' num2str(lat1) ' lon ' num2str(lon0) '-' num2str(lon1)]);
        fprintf(fstlog,'%+8s', num2str(total_cells,'% 8d'));
        fprintf(fstlog,'%+8s', num2str((valid_cells),'% 8d'));
        fprintf(fstlog,'%+12s', num2str(mean(mscal),'% 12.4f'));
        fprintf(fstlog,'%+12s', num2str(min(mscal),'% 12.4f'));
        fprintf(fstlog,'%+12s', num2str(max(mscal),'% 12.4f'));
        fprintf(fstlog,'%+12s', num2str(mean(msstd),'% 12.4f'));
        fprintf(fstlog,'%+12s', num2str(min(msstd),'% 12.4f'));
        fprintf(fstlog,'%+12s\n', num2str(max(msstd),'% 12.4f')); 

        %% contune to next mascon
        continue_process_mascon = continue_process_mascon + 1;
    end
    if idx_lon1 == idx_lon && idx_lat1 == idx_lat
        continue_process_mascon = 0;
    elseif idx_lon1 == idx_lon && idx_lat1 < idx_lat
        next_lat = 1;
        idx_lon1 = 0;
    else
        idx_lon0 = idx_lon1+1;
        next_lat = 0;
    end                                          
end
%% completed processing scale factors
        
%% save the results to the output file
if daily_files == 1
    fname = 'daily';
else
    fname = 'monthly';
end
out_scale_file = [output_path '\LSM_scale_gain_factors_lcc_grid_from_' fname '_EALCO.dat'];
fid = fopen(out_scale_file,'w');
    fwrite(fid,scal','float32');
    fwrite(fid,sstd','float32');
fclose(fid);
copyfile(hdr_file,[output_path '\LSM_scale_gain_factors_lcc_grid_from_' fname '_EALCO.hdr']);
    
ttime = toc;

disp(['>>> Elapsed time: ', num2str(ttime),' seconds. <<<'])

msgbox('Completed scale factor computation successfully!');

logmsg(flog,handles,['>>> Elapsed time: ', num2str(ttime),' seconds. <<<']); 
logmsg(flog,handles,'Completed scale gain factor computation successfully!'); 

fclose(flog);
fclose(fstlog);



function txtGainFactorFile_Callback(hObject, eventdata, handles)
% hObject    handle to txtGainFactorFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtGainFactorFile as text
%        str2double(get(hObject,'String')) returns contents of txtGainFactorFile as a double


% --- Executes during object creation, after setting all properties.
function txtGainFactorFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtGainFactorFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowseGainFactorFile.
function pbBrowseGainFactorFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbBrowseGainFactorFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile( ...
    {  '*Scale*.*','Specify a Canada scale factor LCC grid file (*Scale*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify a Canada scale factor LCC grid file');
 if  filename > 0  
    set(handles.txtGainFactorFile, 'String', [pathname filename]);
 end 

% --- Executes on button press in chkApplyGainFactors.
function chkApplyGainFactors_Callback(hObject, eventdata, handles)
% hObject    handle to chkApplyGainFactors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkApplyGainFactors


% --- Executes on button press in pbCalculateEALCOMonthlyTWSA.
function pbCalculateEALCOMonthlyTWSA_Callback(hObject, eventdata, handles)
% hObject    handle to pbCalculateEALCOMonthlyTWSA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% record the start time
tic;

%% get GRACE TWS input data file name
gtws_file = get(handles.txtGraceTwsFile, 'String');
if isempty(gtws_file)
    msgbox('Specify the GRACE TWS Anomaly NetCDF file.');
    return;
end

%% get the output folder and ENVI header file
hdr_file = get(handles.txtOutputFile, 'String');
out_file = get(handles.txtOutputFile, 'String');
if isempty(out_file)
    msgbox('Specify the output data ENVI header file.');
    return;
end
%% retrive the output data filder
[output_path, outfname, ext] = fileparts(out_file);

%% get EALCO TWS daily absolute anomaly data file list
file_list = get(handles.lstTerresTwsTimeSeries,'String');
if iscell(file_list)
    nf = length(file_list);
else
    nf = 1;
end
if nf < 1
    msgbox('Add at least one LSM daily TWS data file to the LSM TWS Time Series list!');
    return;
end

%% create a time and time string list for output file names
daily_files = get(handles.popFileType, 'Value'); % used to determine date formate in the file name
time_list = zeros(nf,1);
time_str_list = cell(nf,1);
for i = 1:nf
    if ischar(file_list)
        ttwsi_file = file_list;
    else
        ttwsi_file = char(file_list{i});
    end
    [time_list(i), time_str_list{i}] = daysFromFileName(ttwsi_file,daily_files);
end

%% get the meta data parameters and processing options
switch get(handles.popDataType, 'Value')
    case 1      % daily data in int16
        data_type =  'int16';
    case 2      % monthly data in float32
        data_type =  'float32';
end

row = str2double(get(handles.txtRow, 'String'));    % row number of the GRACE TWS data
col = str2double(get(handles.txtCol, 'String'));    % col number of the GRACE TWS data
band = str2double(get(handles.txtBand, 'String'));  % the number of months 

switch get(handles.popNoDataValue, 'Value')
    case 1      
        no_data =  -32760;
    case 2      
        no_data =  9999;
end

%% open log file to record computation progresses
logfile = [output_path '\LSM_monthly_TWSA_process.log'];
if exist(logfile,'file')
    delete logfile;
end
flog = fopen(logfile,'wt');

logmsg(flog,handles,'Processing started.');

%% initialize the output data containers 
dtws = ones(row, col, nf)*no_data;
base = ones(row, col)*no_data;
base_days = ones(row, col)*no_data;

%% read the LSM TWS time series data
for i = 1:nf
    dfile = char(file_list(i));
    fid = fopen(dfile);
        dtws(:,:,i) = (fread(fid,[col, row], data_type))';
    fclose(fid);
end
logmsg(flog,handles,'Read data from all daily data files sucessfully.');

%% calculate the baseline values from 2002 April 5 to Dec 31 2016
[r,c,b] = size(dtws);
for i = 1:r
    for j = 1:c
        tsij = dtws(i,j,:);
        inval = find(tsij == -32760 | tsij < 0);
        if length(inval) >= b/3
            base(i,j) = -32760;
            base_days(i,j) = 0;
        else
            tsij(inval) = [];
            base(i,j) = mean(tsij(:));
            base_days(i,j) = b - length(inval);
        end
    end
end

% %% save the baseline data
% out_base_file = [output_path '\LSM_baseline_2002-04-05_2016-12-31_lcc_grid.dat'];
% fid = fopen(out_base_file,'w');
%     fwrite(fid,base','float32');
%     fwrite(fid,base_days','float32');    
% fclose(fid);
% copyfile(hdr_file,[output_path '\LSM_baseline_2002-04-05_2016-12-31_lcc_grid.hdr']);
% logmsg(flog,handles,'Completed LSM baseline value grid computation sucessfully.');

%% calculate EALCO daily TWSA 
for i = 1:nf
    tsi = dtws(:,:,i);
    idx = find(tsi > 0);
    tsi(idx) = tsi(idx) - base(idx);
    out_daily_file = [output_path '\LSM_DTWSA_' time_str_list{i} '.dat'];
    fid = fopen(out_daily_file,'w');
        fwrite(fid,tsi','float32');   
    fclose(fid);
    % copyfile(hdr_file,[output_path '\LSM_DTWSA_' time_str_list{i} '.hdr']);
end
logmsg(flog,handles,'Completed daily anomaly computations sucessfully.');

%% calculate the monthy EALCO TWSA
%% read the GRACE TWS anomaly data from NETCDF file
day = ncread(gtws_file,'time');
day_bounds = ncread(gtws_file,'time_bounds');
for i = 1:length(day)
    day_bounds_i = day_bounds(:,i);
    day_start = day_bounds_i(1);
    day_end = day_bounds_i(2);
    
    twsi = ones(r, c)*no_data;
    stdi = ones(r, c)*no_data;

    idx_days = find(time_list>=day_start & time_list<=day_end);       
    if ~isempty(idx_days)
        % mfile_list = file_list(idx_days);
        mtime_str_list = time_str_list(idx_days);
        tsij = dtws(:,:,idx_days);
        %% calculate the mean and std of each grid cell
        for k = 1:r
            for j = 1:c
                %% remove possible invalid values
                dij = tsij(k,j,:);
                idx = find(dij == -32760 | dij == 9999);
                if~isempty(idx)
                    dij(idx) = [];
                end
                if ~isempty(dij)
                    twsi(k,j) = mean(dij);
                    stdi(k,j) = std(dij);
                end
            end
        end
        idx = find(twsi > 0);
        twsi(idx) = twsi(idx) - base(idx);

        %% determine output file name according to the time
        time_str = [mtime_str_list{1} '_' mtime_str_list{end}];
        out_monthly_file = [output_path '\LSM_MTWSA_' time_str '.dat'];
        fid = fopen(out_monthly_file,'w');
            fwrite(fid,twsi','float32');    % save the EALCO monthly TWSA
            fwrite(fid,stdi','float32');    % save the STD of the EALCO monthly TWS
        fclose(fid);
        copyfile(hdr_file,[output_path '\LSM_MTWSA_' time_str '.hdr']);
    end
end
logmsg(flog,handles,'Completed LSM monthly TWSA computations sucessfully.');

ttime = toc;

disp(['>>> Elapsed time: ', num2str(ttime),' seconds. <<<'])

msgbox('Completed LSM monthly TWSA computation successfully!');

logmsg(flog,handles,['>>> Elapsed time: ', num2str(ttime),' seconds. <<<']); 
logmsg(flog,handles,'Completed LSM monthly TWSA computation successfully!'); 

fclose(flog);


% --- Executes on button press in pbTemporalDataAssimilationBySCVCM.
function pbTemporalDataAssimilationBySCVCM_Callback(hObject, eventdata, handles)
% hObject    handle to pbTemporalDataAssimilationBySCVCM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% record the start time
tic;

%% get GRACE TWS input data file name
gtws_file = get(handles.txtGraceTwsFile, 'String');
if isempty(gtws_file)
    msgbox('Specify the GRACE TWS data file.');
    return;
end

%% get JPL mascon grid id data file name
mask_file = get(handles.txtGraceTwsIdMaskFile, 'String');
if isempty(mask_file)
    msgbox('Specify the mascon id mask data file.');
    return;
end

%% get the EALCO baseline data file
baseline_file = get(handles.txtEALCOBaselineDataFile, 'String');
if isempty(baseline_file)
    msgbox('Specify the EALCO Baseline data file.');
    return;
end

%% get the output header file in ENVI data format
hdr_file = get(handles.txtOutputFile, 'String'); 
out_file = get(handles.txtOutputFile, 'String');
if isempty(out_file)
    msgbox('Specify the output header file.');
    return;
end
[output_path, outfname, ext] = fileparts(out_file);
std_twsa_output_path = [output_path '\STD_TWSA'];
std_gswsa_output_path = [output_path '\STD_GSWSA'];
if exist(std_twsa_output_path, 'dir') ~= 7
    mkdir std_twsa_output_path;
end
if exist(std_gswsa_output_path, 'dir') ~= 7
    mkdir std_gswsa_output_path;
end

%% get EALCO model TWS data file list
daily_files = get(handles.popFileType, 'Value');
if daily_files ~= 1
    msgbox('Temporal data assimilation is for daily LSM TWSA! Specify daily LSM TWSA files.');
    return;
end
file_list = get(handles.lstTerresTwsTimeSeries,'String');
if iscell(file_list)
    nf = length(file_list);
else
    nf = 1;
end
if nf<1
    msgbox('Add at least one EALCO TWS data file to the EALCO TWS Time Series list!');
    return;
end

%% create a time and time string list for the EALCO TWS input file list
time_list = zeros(nf,1);
time_str_list = cell(nf,1);
for i = 1:nf
    if ischar(file_list)
        ttwsi_file = file_list;
    else
        ttwsi_file = char(file_list{i});
    end
    if daily_files == 1
        [time_list(i), time_str_list{i}] = daysFromFileName(ttwsi_file,1);
    else
        [time_list(i), time_str_list{i}] = daysFromFileName(ttwsi_file,3);
    end
end
%% note: the time stamps is in days from Jan. 1, 2002

%% get GRACE TWS anomaly data file list
gfile_list = get(handles.lstGraceTwsTimeSeries,'String');
if iscell(gfile_list)
    ngf = length(gfile_list);
else
    ngf = 1;
end
if ngf<1
    msgbox('Add at least one GRACE TWS data file to the list!');
    return;
end
%% create time and time string list for GRACE TWS anomaly files
gtime_list = zeros(ngf,1);
gtime_str_list = cell(ngf,1);
for i = 1:ngf
    if ischar(gfile_list)
        gtwsi_file = gfile_list;
    else
        gtwsi_file = char(gfile_list{i});
    end
    [gtime_list(i), gtime_str_list{i}] = daysFromFileName(gtwsi_file,3);
end

%% get the meta data parameters and processing options
switch get(handles.popDataType, 'Value')
    case 1      % daily data in int16
        data_type =  'int16';
    case 2      % monthly data in float32
        data_type =  'float32';
end
switch get(handles.popNoDataValue, 'Value')
    case 1      
        no_data =  -32760;
    case 2      
        no_data =  9999;
end

weight = get(handles.popWeightMethod, 'Value');

row = str2double(get(handles.txtRow, 'String'));    % row number of the GRACE TWS data
col = str2double(get(handles.txtCol, 'String'));    % col number of the GRACE TWS data
band = str2double(get(handles.txtBand, 'String'));  % the number of months 

%% read in the mascon id mask data
fid = fopen(mask_file);
    maskid = (fread(fid,[col, row],'float32')');
fclose(fid);

%% read in EALCO baseline data
fid = fopen(baseline_file);
    baseline = (fread(fid,[col, row],'float32')');
fclose(fid);

%% read the GRACE TWS data from NETCDF file
%ncdisp(gtws_file)
% gtws = ncread(gtws_file,'lwe_thickness');
% gtws_std = ncread(gtws_file,'uncertainty');
% lat = ncread(gtws_file,'lat');
% lon = ncread(gtws_file,'lon');
% lat_bounds = ncread(gtws_file,'lat_bounds');
% lon_bounds = ncread(gtws_file,'lon_bounds');
day = ncread(gtws_file,'time');
day_bounds = ncread(gtws_file,'time_bounds');

%% open log file to record computation progresses
logfile = [output_path '\TDA_Processing.log'];
if exist(logfile,'file')
    delete logfile;
end
flog = fopen(logfile,'wt');

%% open a statistic log file to record computation results
st_logfile = [output_path '\TDA_Statistics.log'];
if exist(st_logfile,'file')
    delete st_logfile;
end
fstlog = fopen(st_logfile,'wt');

logmsg(flog,handles,'Processing started.');
logmsg(fstlog,handles,'Processing Statistics:');

%% initialize the output data container for all monthly datasets
%% the output data is presented in LCC coordinate system
%all_atws = ones(row, col, nf)*no_data;
%all_gws = ones(row, col, nf)*no_data;
%all_atws_std = ones(row, col, nf)*no_data;
%all_gws_std = ones(row, col, nf)*no_data;

%% determe how many mascons are in the selected test area
maskvalues = unique(maskid);
maskvalues(maskvalues==0)=[];   % remove the zero 
% maskvalues = [234,235,236,271,272,305,306,307,339,340]';
% maskvalues = [234,235,236,271,272,305,306,307,339,340]';
% % maskvalues = [273,274,307,308,341]';
% % maskvalues = [273,274,308]';
% % maskvalues = [162,163,203,204,205,206]';
% nmask = length(maskvalues);

if band > ngf
    band = ngf;
end

%% process data band by band and monscon by monscon
for i=1:band         %% note: each band represnts a monthly dataset
    %% determine the day ranges for the monthly dataset i
    day_bounds_i = day_bounds(:,i);
    day_start = day_bounds_i(1);
    day_end = day_bounds_i(2);
    dayi = day(i);
    
    idx_gmonth = find(gtime_list>dayi-5 & gtime_list<dayi+5);
    if isempty(idx_gmonth)
        continue;
    else %% read the monthly GRACE TWS from file
        %% read in monthly average GRACE TWS data directly from file
        if ischar(gfile_list)
            gtwsi_file = gfile_list;
        else
            gtwsi_file = char(gfile_list{idx_gmonth});
        end
        %% read in the GRACE TWS from the monthly file
        fid = fopen(gtwsi_file);
            gtwsi = (fread(fid,[col, row], data_type)');
            gstdi = (fread(fid,[col, row], data_type)');
        fclose(fid);
%         %% read in the GSWSA from the monthly file
%         gswsi_file = replace(gtwsi_file,'MTWSA','MGWSA');
%         fid = fopen(gswsi_file);
%             gwsi = (fread(fid,[col, row], data_type)');
%             stdi = (fread(fid,[col, row], data_type)');
%         fclose(fid);
    end

    logmsg(flog,handles,' ');
    logmsg(flog,handles,['Processing GRACE TWS for the day: ' num2str(dayi)]);
    logmsg(flog,handles,' ');    
    
    %% retreive the daily EALCO TWS files for the GRACE TWS monthly dataset   
    idx_days = find(time_list>=day_start & time_list<=day_end);  % The start and end date should be included     
    if isempty(idx_days)        
        logmsg(flog,handles,['No matched daily EALCO TWS files found for the days from ' num2str(day_start) ' to ' num2str(day_end)]);
        continue;
    else   
        %% continue processing the matched GRACE and EALCO data    
        logmsg(flog,handles,[num2str(length(idx_days)) ' daily EALCO TWS files found for the days from ' num2str(day_start) ' to ' num2str(day_end)]);
        mfile_list = file_list(idx_days);
        mtimestr = time_str_list(idx_days);
        
        %% initialize the output data array
        mf = length(idx_days);
        m_atws = ones(row, col, mf)*no_data;
        m_atws_std = ones(row, col, mf)*no_data;
        m_gws = ones(row, col, mf)*no_data;
        m_gws_std = ones(row, col, mf)*no_data;   
        
        %% compute monthly average of EALCO TWS and the sacal factors
        [ttwsi,tstdi,dtwsi] = averageFromInputFiles(mfile_list, row, col, data_type, no_data);
        
        for pi = 1:row
            for pj =1:col
                if ismember(maskid(pi,pj),maskvalues) %% process the specified mascon id only
                    if ttwsi(pi,pj) == no_data || gtwsi(pi,pj) == no_data  % if no monthly EALCO TWSA data available
                        continue;
                    else
                        mttws = ttwsi(pi,pj);               % monthly EALCO TWSA
                        mgtws = gtwsi(pi,pj);               % spatially downscaled GRACE TWSA
                        mgstd = gstdi(pi,pj);               % std of the spatially downscaled GRACE TWSA
                        
                        basev = baseline(pi,pj);            % baseline value
                        
                        % mgws  = gwsi(pi,pj);              % monthly GSWSA. not used
                        % mstd  = stdi(pi,pj);              % std of the monthly GSWSA, not used

                        month_data = dtwsi(pi,pj,:);        % daily EALCO TWSA (one month)
                        month_data = month_data(:);
                        %vaild_idx = find(month_data ~= no_data);
                        vaild_idx = find(abs(month_data) <= 9999);
                        if length(vaild_idx) ~= length(month_data)
                            disp("Start debug!");
                        end
                        temp = month_data(vaild_idx);       % exclude possible no data points
%                         %scale = temp(:)./mttws;            % scale factors
%                         scale = abs(temp(:))./mean(abs(temp));   % scale factors
%                         mgtws = scale*mgtws;                % scaled (spatially downscaled) GRACE TWSA

                        scale = (temp+basev)./mean(temp+basev);
                        mgtws = scale*(mgtws+basev) - basev;  % scaled (spatially downscaled) GRACE TWSA
                        
                        mgstd = scale*mgstd;                % std of the scaled (spatially downscaled) GRACE TWSA
                        mttws = temp;                       % valid daily EALCO TWSA
                        mtstd = ones(length(mttws),1)*15;   % assign std to the daily EALCO TWSA

                        % mgws  = scale*mgws;               % scaled GSWSA
                        % mstd  = scale*mstd;               % std of the scaled GSWSA

                        mgws  = scale*0;                    % vertual zero observations for GSWSA 
                        mstd  = sqrt(mgstd.^2+mtstd.^2);    % std of the vertual zero observations for GSWSA 

                        %% downscale and assimilation
                        [amgtws, amgstd, amgws, amstd, iteration] = temporalAssimilateMasconBySCVCM2(mgtws, mgstd, mttws, mtstd, mgws, mstd, weight); 
                        logmsg(flog,handles,['Processed pixels (' num2str(pi) ' ' num2str(pj)  ') for the days from ' num2str(day_start) ' to ' num2str(day_end) ' SCVCM iteration number: ' num2str(iteration) ]);

                        if max(abs(amgws)>200)
                            disp("Start debug2!");
                        end
                        
                        %% put the assimilated results to output arrays
                        month_data(vaild_idx) = amgtws;
                        m_atws(pi,pj,:) = month_data;
                        %all_atws(pi,pj,idx_days) = month_data;
                        
                        month_data(vaild_idx) = amgstd;  
                        m_atws_std(pi,pj,:) = month_data;
                        %all_atws_std(pi,pj,idx_days) = month_data; 
                        
                        month_data(vaild_idx) = amgws;
                        m_gws(pi,pj,:) =  month_data;
                        %all_gws(pi,pj,idx_days) =  month_data;
                        
                        month_data(vaild_idx) = amstd;
                        m_gws_std(pi,pj,:) =  month_data;
                        %all_gws_std(pi,pj,idx_days) =  month_data;
                    end
                end
    
            end
        end
        %% save the output data to files, one file per day
        
        for k=1:mf
            time_str = char(mtimestr(k));
            out_file = [std_twsa_output_path '\STD_DTWSA_' time_str '.dat'];
            fid = fopen(out_file,'w');
                fwrite(fid,m_atws(:,:,k)','float32');
                fwrite(fid,m_atws_std(:,:,k)','float32');
            fclose(fid);
            copyfile(hdr_file,[std_twsa_output_path '\STD_DTWSA_'  time_str '.hdr']);
            
            out_file = [std_gswsa_output_path '\STD_DGWSA_' time_str '.dat'];
            fid = fopen(out_file,'w');
                fwrite(fid,m_gws(:,:,k)','float32');
                fwrite(fid,m_gws_std(:,:,k)','float32');
            fclose(fid);
            copyfile(hdr_file,[std_gswsa_output_path '\STD_DGSWSA_'  time_str '.hdr']);
        end
    end
end
   

%% Output all assimilated daily results into one file
% out_file = [output_path '\DA_All_DTWSA.dat'];
% fid = fopen(out_file,'w');
% for k=1:nf
%     fwrite(fid,all_atws(:,:,k)','float32');
% end
% fclose(fid);

% out_file = [output_path '\DA_All_DTWSA_STD.dat'];
% fid = fopen(out_file,'w');
% for k=1:nf
%     fwrite(fid,all_atws_std(:,:,k)','float32');
% end
% fclose(fid);


% out_file = [output_path '\DA_All_DGWSA.dat'];
% fid = fopen(out_file,'w');
% for k=1:nf
%     fwrite(fid,all_gws(:,:,k)','float32');
% end
% fclose(fid);

% out_file = [output_path '\DA_All_DGWSA_STD.dat'];
% fid = fopen(out_file,'w');
% for k=1:nf
%     fwrite(fid,all_gws_std(:,:,k)','float32');
% end
% fclose(fid); 

ttime = toc;

disp(['>>> Elapsed time: ', num2str(ttime),' seconds. <<<'])

msgbox('Completed temporal data assimilation successfully!');

logmsg(flog,handles,['>>> Elapsed time: ', num2str(ttime),' seconds. <<<']); 
logmsg(flog,handles,'Completed temporal data assimilation successfully!'); 

fclose(flog);
fclose(fstlog);

% --- Executes on button press in pbSpatialDABySCVCM1.
function pbSpatialDABySCVCM1_Callback(hObject, eventdata, handles)
% hObject    handle to pbSpatialDABySCVCM1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% record the start time
tic;

%% get GRACE TWS input data file name
gtws_file = get(handles.txtGraceTwsFile, 'String');
if isempty(gtws_file)
    msgbox('Specify the GRACE TWS Anomaly NetCDF file.');
    return;
end

%% get the output data ENVI header file
hdr_file = get(handles.txtOutputFile, 'String');
out_file = get(handles.txtOutputFile, 'String');
if isempty(out_file)
    msgbox('Specify the output data ENVI header file.');
    return;
end
%% retrive the output data filder
[output_path, outfname, ext] = fileparts(out_file);

%% get EALCO TWS anomaly data file list
daily_files = get(handles.popFileType, 'Value');
if daily_files == 1
    fname = 'daily';
else
    fname = 'monthly';
end
file_list = get(handles.lstTerresTwsTimeSeries,'String');
if iscell(file_list)
    nf = length(file_list);
else
    nf = 1;
end
if nf<1
    msgbox('Add at least one LSM TWS data file to the LSM TWS Time Series list!');
    return;
end
%% create a time and time string list for the EALCO TWS input file list
time_list = zeros(nf,1);
time_str_list = cell(nf,1);
for i = 1:nf
    if ischar(file_list)
        ttwsi_file = file_list;
    else
        ttwsi_file = char(file_list{i});
    end
    if daily_files == 1
        [time_list(i), time_str_list{i}] = daysFromFileName(ttwsi_file,1);
    else
        [time_list(i), time_str_list{i}] = daysFromFileName(ttwsi_file,3);
    end
end
%% note: the time stamps is in days from Jan. 1, 2002

%% get the meta data parameters and processing options
switch get(handles.popDataType, 'Value')
    case 1      % daily data in int16
        data_type =  'int16';
    case 2      % monthly data in float32
        data_type =  'float32';
end
switch get(handles.popNoDataValue, 'Value')
    case 1      
        no_data =  -32760;
    case 2      
        no_data =  9999;
end
%% weighting method for EALCO TWS Anomalies
weight = get(handles.popWeightMethod, 'Value');

row = str2double(get(handles.txtRow, 'String'));    % row number of the GRACE TWS data
col = str2double(get(handles.txtCol, 'String'));    % col number of the GRACE TWS data
band = str2double(get(handles.txtBand, 'String'));  % the number of months 

%% get the precalculated scale factors file name for Canada LCC grid
apply_gain_factors = get(handles.chkApplyGainFactors, 'Value');
scale_file = get(handles.txtGainFactorFile, 'String');

%% add geo reference coordinates to EALCO TWS
out_lcc_latlon_file = [output_path '\Canada_Lat_Lon_lcc_grid.mat'];
if ~exist(out_lcc_latlon_file,'file')
    lamb = Lambert(...          % Define Canada LCC projection parameters
        [49.0, 77.0], ...       % Latitude of the first and second standard parallel
        0.0, -95, ...           % Central latitude and longitude
        0.0, 0.0, ...           % False easting and northing
        6378137.0, 298.257222100882711243 ...     % Ellipsoid: Semi-major axis and inverse flattening
    );

    lcc_ulx_e = -2600000.0;     
    lcc_uly_n = 10500000.0;
    % lcc_lrx_e = 3100000.0;      
    % lcc_lry_n = 5700000.0;
    lcc_dx = 5000.0;            
    lcc_dy = 5000.0;
    lcc_x = zeros(row, col);   
    lcc_y = zeros(row, col);
    lcc_lat = zeros(row, col);   
    lcc_lon = zeros(row, col);

    for i = 1:row
        for j =1:col
            lcc_x(i,j) = lcc_ulx_e + j*lcc_dx - lcc_dx/2;
            lcc_y(i,j) = lcc_uly_n - i*lcc_dy + lcc_dy/2;
            [lcc_lat(i,j), lcc_lon(i,j)] = lamb.cartesian2geographic(lcc_x(i,j), lcc_y(i,j));
        end
    end

    save(out_lcc_latlon_file,'lcc_x', 'lcc_y', 'lcc_lat', 'lcc_lon');
else
    load(out_lcc_latlon_file);
end

%% determine the geolocation area of CANADA
idx = find(lcc_lon < 0);
lcc_lon(idx) = lcc_lon(idx)+360;
minLat = min(lcc_lat(:));   maxLat = max(lcc_lat(:));
minLon = min(lcc_lon(:));   maxLon = max(lcc_lon(:));

%% read the GRACE TWS anomaly data from NETCDF file
% ncdisp(gtws_file)
gtws = ncread(gtws_file,'lwe_thickness');
gtws_std = ncread(gtws_file,'uncertainty');
% lat = ncread(gtws_file,'lat');
% lon = ncread(gtws_file,'lon');
day = ncread(gtws_file,'time');
lat_bounds = ncread(gtws_file,'lat_bounds');
lon_bounds = ncread(gtws_file,'lon_bounds');
day_bounds = ncread(gtws_file,'time_bounds');

%% determine the GRACE TWS anomaly data index range for CANADA area
[idx_lon,idx_lat,b] = size(gtws);
for i = 1:idx_lat
    if minLat >= lat_bounds(1,i) && minLat < lat_bounds(2,i)
        idx_lat_start = i;
        break;
    end
end
for i = 1:idx_lat
    if maxLat >= lat_bounds(1,i) && maxLat < lat_bounds(2,i)
        idx_lat_end = i;
        break;
    end
end
for i = 1:idx_lon
    if minLon >= lon_bounds(1,i) && minLon < lon_bounds(2,i)
        idx_lon_start = i;
        break;
    end
end
for i = 1:idx_lon
    if maxLon >= lon_bounds(1,i) && maxLon < lon_bounds(2,i)
        idx_lon_end = i;
        break;
    end
end
%% retreive GRACE TWS Anomaly and their uncerntainty data for Canada area
gtws = gtws(idx_lon_start:idx_lon_end,idx_lat_start:idx_lat_end,:);
gtws_std = gtws_std(idx_lon_start:idx_lon_end,idx_lat_start:idx_lat_end,:);
% lat = lat(idx_lat_start:idx_lat_end);
% lon = lon(idx_lon_start:idx_lon_end);
lat_bounds = lat_bounds(:,idx_lat_start:idx_lat_end);
lon_bounds = lon_bounds(:,idx_lon_start:idx_lon_end);


%% open log file to record computation progresses
logfile = [output_path '\PM1_Process_' fname '.log'];
if exist(logfile,'file')
    delete logfile;
end
flog = fopen(logfile,'wt');

%% open a statistic log file to record some statistics of results
st_logfile = [output_path '\PM1_Statistics_' fname '.log'];
if exist(st_logfile,'file')
    delete st_logfile;
end
fstlog = fopen(st_logfile,'wt');

logmsg(flog,handles,'Processing started.');
logmsg(fstlog,handles,'Processing Statistics:');

%% initialize the output data container for all monthly datasets
%% the output data is presented in LCC coordinate system
if band > nf
    band = nf;
end
all_atws = ones(row, col, band)*no_data;
all_gtws = ones(row, col, band)*no_data;
all_ttws = ones(row, col, band)*no_data;
all_stws = ones(row, col, band)*no_data;
all_gws = ones(row, col, band)*no_data;

masconid = zeros(row, col);

if daily_files == 1 && apply_gain_factors == 1 && isempty(scale_file)
    all_scal = ones(row, col, band)*no_data;
end

%% determine the size of GRACE TWS anomaly input data
%% note: the size of gtws is for Canada coverage only
[idx_lon,idx_lat,b] = size(gtws);  % b may be >= band

%% compute or read the baseline values 
out_grace_baseline_file = [output_path '\Canada_GRACE_TWS_anomaly_' num2str(band) '_month_baseline.mat'];
if ~exist(out_grace_baseline_file,'file')
    gbase = zeros(idx_lon,idx_lat);
    for i = 1:idx_lon
        for j = 1:idx_lat
            gbase(i,j) = mean(gtws(i,j,1:band));
        end
    end
    save(out_grace_baseline_file, 'gbase');
else
    load(out_grace_baseline_file);
end

%% read the gain factors if selected to consider the leakage errors
if apply_gain_factors == 1 && ~isempty(scale_file)
    fid = fopen(scale_file);
        scale = (fread(fid,[col, row], data_type)');
        scale_std = (fread(fid,[col, row], data_type)');
    fclose(fid);   
end

%% process data band by band for each monscon
%% note: each band represnts a monthly dataset
for i = 1:band         
    %% determine the day ranges for the monthly dataset i
    day_bounds_i = day_bounds(:,i);
    day_start = day_bounds_i(1);
    day_end = day_bounds_i(2);
    dayi = day(i);

    logmsg(flog,handles,' ');
    logmsg(flog,handles,['Processing GRACE TWS anomaly for the day: ' num2str(dayi)]);
    logmsg(flog,handles,' ');    
    
    %% read or compute monthly averaged EALCO TWS anomaly data
    if daily_files == 1
        %% retreive the daily EALCO TWS files for the GRACE TWS monthly dataset   
        idx_days = find(time_list>=day_start & time_list<=day_end);       
        if isempty(idx_days)        
            logmsg(flog,handles,['No matched daily LSM TWS files found for the days from ' num2str(day_start) ' to ' num2str(day_end)]);
            continue;
        else   
            %% continue processing the matched GRACE and EALCO data    
            logmsg(flog,handles,[num2str(length(idx_days)) ' daily LSM TWS files found for the days from ' num2str(day_start) ' to ' num2str(day_end)]);
            mfile_list = file_list(idx_days);
            mtime_str_list = time_str_list(idx_days);
            %% compute monthly average of EALCO TWS from the input file list
            if apply_gain_factors == 1 && isempty(scale_file)
                [ttwsi,tstdi,dtwsi] = averageFromInputFiles(mfile_list, row, col, data_type, no_data);
            else
                [ttwsi,tstdi] = averageFromInputFiles(mfile_list, row, col, data_type, no_data);
            end
            idx_middle = fix(length(idx_days)/2+0.5);
            %% determine output file name according to the time
            time_str = char(mtime_str_list(idx_middle));
        end
    else
        % idx_month = find(time_list>day_start+1 & time_list<day_end-1);
        idx_month = find(time_list>dayi-5 & time_list<dayi+5);
        if isempty(idx_month)        
            logmsg(flog,handles,['No matched monthly LSM TWS files found for the days from ' num2str(day_start) ' to ' num2str(day_end)]);
            continue;
        else   
            %% continue processing the matched GRACE and EALCO data    
            logmsg(flog,handles,[num2str(length(idx_month)) ' monthly LSM TWS file found for the days from ' num2str(day_start) ' to ' num2str(day_end)]);
            %% read in monthly average EALCO TWS data directly from file
            if ischar(file_list)
                ttwsi_file = file_list;
            else
                ttwsi_file = char(file_list{idx_month});
            end
            %% read in LSM TWSA data from the data file
            fid = fopen(ttwsi_file);
                ttwsi = (fread(fid,[col, row], data_type)');
                tstdi = (fread(fid,[col, row], data_type)');
            fclose(fid);
            %% determine output file name according to the time
            time_str = char(time_str_list(idx_month));
        end
    end
   
    %% retreive GRACE TWS data for the same month
    gtwsi = gtws(:,:,i)';
    gtwsi_std = gtws_std(:,:,i)';
    %% deduct the base line value
    gtwsi = (gtwsi - gbase');  
    %% convert the unit from cm to mm
    gtwsi = gtwsi*10;
    gtwsi_std = gtwsi_std*10;
    
    fprintf(fstlog,'%s\n', ' ');
    fprintf(fstlog,'%s\n',['Processing EALCO TWS for the days from ' num2str(day_start) ' to ' num2str(day_end)]);
    fprintf(fstlog,'%s\n','   ID  #OfTTWS    MinTWSA     MaxTWSA    MeanTWSA      GTWSA     MeanSCAL   MeanSTWSA   MeanATWSA   MeanGSWSA      TTStd0      TTStd1      GTStd0      GTStd1       STStd      ATStd       GWStd        MinD        MaxD       MeanD    Iteration');

    %% initialize the output data after assimilation vy copy the ttws
    atwsio = ones(row, col)*no_data;
    gtwsio = ones(row, col)*no_data;
    ttwsio = ones(row, col)*no_data;
    stwsio = ones(row, col)*no_data;
    gwsio = ones(row, col)*no_data;
    
    atwsio_std = ones(row, col)*no_data;
    gtwsio_std = ones(row, col)*no_data;
    ttwsio_std = ones(row, col)*no_data;
    gwsio_std = ones(row, col)*no_data;
        
    if daily_files == 1 && apply_gain_factors == 1 && isempty(scale_file)
        scalio = ones(row, col)*no_data;
        scalio_std = ones(row, col)*no_data;
    end
    
    %% start a loop to downscale and assimilate GRACE TWS to EALCO TWS
    continue_process_mascon = 1; 
    idx_lat1 = 0; idx_lon1 = 0;
    next_lat = 1;
    while continue_process_mascon > 0 && idx_lat1 <= idx_lat-1
        %% look for a GRACE TWS mascon data and its Lat Lon bounds
        if next_lat == 1
            idx_lat0 = idx_lat1+1; idx_lon0 = idx_lon1 + 1;
            mgtws = gtwsi(idx_lat0, idx_lon0);
            mgstd = gtwsi_std(idx_lat0,idx_lon0);
            for j=idx_lat0:idx_lat
                if gtwsi(j,idx_lon0) ~= mgtws
                    idx_lat1 = j-1;
                    break;
                elseif j == idx_lat
                    idx_lat1 = j;
                end
            end
        else
            idx_lon0 = idx_lon1 + 1;
            mgtws = gtwsi(idx_lat0, idx_lon0);
            mgstd = gtwsi_std(idx_lat0,idx_lon0);
        end
        for j=idx_lon0:idx_lon
            if gtwsi(idx_lat0,j) ~= mgtws
                idx_lon1 = j-1;
                break;
            elseif j == idx_lon
                idx_lon1 = j;
            end
        end
        lat0 = lat_bounds(1,idx_lat0);
        lat1 = lat_bounds(2,idx_lat1);
        lon0 = lon_bounds(1,idx_lon0);
        lon1 = lon_bounds(2,idx_lon1);

        %% check all valuse within the mascon are same (for debug test only)
        % mgtws = gtwsi(idx_lat0:idx_lat1,idx_lon0:idx_lon1);
        % mgstd = gtwsi_std(idx_lat0:idx_lat1,idx_lon0:idx_lon1);
        G = unique(mgtws); G_std = unique(mgstd);
        ng = length(G); ng_std = length(G_std);
        logmsg(flog,handles,[num2str(ng) ' GRACE TWS Values found for the mascon id ' num2str(continue_process_mascon)]);
        logmsg(flog,handles,[num2str(ng_std) ' GRACE TWS STD found for the mascon id ' num2str(continue_process_mascon)]);
        if ng > 1 
            countg = zeros(ng,1);
            for k=1:ng
                countg(k) = length(find(G == G(k)));
                logmsg(flog,handles,['GRACE TWS Value: ' num2str(G(k)) ' n = ' num2str(countg(k))]);
            end
            [~, idx2] = max(countg);
            G = G(idx2);
            G_std= G_std(idx2);
            logmsg(flog,handles,['Apply GRACE TWS Value ' num2str(G) ' for the mascon id ' num2str(continue_process_mascon)]);
            logmsg(flog,handles,['Apply GRACE TWS STD ' num2str(G_std) ' for the mascon id ' num2str(continue_process_mascon)]);
        end
        
        %% determine the mascon id        
        m_id = continue_process_mascon; 
        %% retreive EALCO TWS within the mascon
        midx = find((lcc_lat>=lat0 & lcc_lat<lat1) & (lcc_lon>=lon0 & lcc_lon<lon1));
        mttws = ttwsi(midx);
        mtstd = tstdi(midx);

        %% remove invalid EALCO TWS
        inval = find(mttws == -32760 | mttws == 9999 | isnan(mttws)); 
        if ~isempty(inval)
            midx(inval) = [];
            mttws(inval) = [];
            mtstd(inval) = [];
        end          
        if isempty(mttws)
            logmsg(flog,handles,' ');
            logmsg(flog,handles,['No valid EALCO TWS found within the mascon ' num2str(continue_process_mascon) ': lat ' num2str(lat0) '-' num2str(lat1) ' lon ' num2str(lon0) '-' num2str(lon1)]);
            %% contune to next mascon
            continue_process_mascon = continue_process_mascon + 1;
        else
            %if ismember(m_id,[235,236,271,272,273,305,306,307,339,340,370])
            if ismember(m_id,[273,274,308])
                logmsg(flog,handles,' ');
                logmsg(flog,handles,['Found ' num2str(length(mttws)) ' EALCO TWS within the mascon ' num2str(continue_process_mascon) ': lat ' num2str(lat0) '-' num2str(lat1) ' lon ' num2str(lon0) '-' num2str(lon1)]);            

                %% assimilate the matched mascon GRACE TWS with EALCO TWS

                %% calculate or retreive scale factors
                if daily_files == 1 && apply_gain_factors == 1 && isempty(scale_file)
                    [mscal,msstd] = scakeFactorFromEALCO(dtwsi,midx);
                    scalio(midx) = mscal;
                    scalio_std(midx) = msstd;
                elseif apply_gain_factors == 1 && ~isempty(scale_file)
                    mscal = scale(midx);
                    msstd = scale_std(midx);      
                end

                %% call assimilation function 
                % [matws,mgws,mstws,se,sg,sk,iteration] = assimilateMasconToGridBySCVCM5(mgtws,mgstd,mttws,mtstd,mscal,msstd,weight);
                [matws, mastd, mgsws, mgstd, mstws, tstd, gstd, iteration] = assimilateMasconBySCVCM(G,G_std,mttws,mtstd,mscal,msstd,weight);
                logmsg(flog,handles,['SCVCM iteration number: ' num2str(iteration) ]);

                %% put the assimilated TWS to the initialized output dataset
                atwsio(midx) = matws; 
                gtwsio(midx) = ones(length(midx),1)*G;
                ttwsio(midx) = mttws;
                stwsio(midx) = mstws;
                gwsio(midx) = mgsws; 
                
                atwsio_std(midx) = mastd;
                gtwsio_std(midx) = ones(length(midx),1)*G_std;
                ttwsio_std(midx) = mtstd;
                gwsio_std(midx) = mgstd;

                %% check the differences between the TTWS and ATWS
                delta = matws - mttws - mgsws;       
                fprintf(fstlog,'%+5s', num2str(m_id,'% 5d'));
                fprintf(fstlog,'%+8s', num2str(length(mttws),'% 8d'));
                fprintf(fstlog,'%+12s', num2str(min(mttws),'% 12.4f'));
                fprintf(fstlog,'%+12s', num2str(max(mttws),'% 12.4f'));

                fprintf(fstlog,'%+12s', num2str(mean(mttws),'% 12.4f'));
                fprintf(fstlog,'%+12s', num2str(G,'% 12.4f'));
                fprintf(fstlog,'%+12s', num2str(mean(mscal),'% 12.4f'));        
                fprintf(fstlog,'%+12s', num2str(mean(mstws),'% 12.4f'));
                fprintf(fstlog,'%+12s', num2str(mean(matws),'% 12.4f'));
                fprintf(fstlog,'%+12s', num2str(mean(mgsws),'% 12.4f'));

                fprintf(fstlog,'%+12s', num2str(mean(mtstd),'% 12.4f'));
                fprintf(fstlog,'%+12s', num2str(tstd,'% 12.4f'));
                fprintf(fstlog,'%+12s', num2str(G_std,'% 12.4f'));
                fprintf(fstlog,'%+12s', num2str(gstd,'% 12.4f'));
                fprintf(fstlog,'%+12s', num2str(mean(msstd),'% 12.4f'));        
                fprintf(fstlog,'%+12s', num2str(mean(mastd),'% 12.4f'));
                fprintf(fstlog,'%+12s', num2str(mean(mgstd),'% 12.4f'));

                fprintf(fstlog,'%+12s', num2str(min(delta),'% 12.4f'));
                fprintf(fstlog,'%+12s', num2str(max(delta),'% 12.4f'));
                fprintf(fstlog,'%+12s', num2str(mean(delta),'% 12.4f'));

                fprintf(fstlog,'%+8s\n', num2str(iteration,'% 8d'));
            end
            %% assign a valid mascon id 
            masconid(midx) = m_id;
            %% contune to next mascon
            continue_process_mascon = continue_process_mascon + 1;
        end
        if idx_lon1 == idx_lon && idx_lat1 == idx_lat
            continue_process_mascon = 0;
        elseif idx_lon1 == idx_lon && idx_lat1 < idx_lat
            next_lat = 1;
            idx_lon1 = 0;
        else
            idx_lon0 = idx_lon1+1;
            next_lat = 0;
        end                                          
    end
    %% completed processing one dataset (a time stamp of one month)
        
    %% put the monthly data into a channel
    all_atws(:,:,i) = atwsio;
    all_gtws(:,:,i) = gtwsio;
    all_ttws(:,:,i) = ttwsio;
    all_stws(:,:,i) = stwsio;
    all_gws(:,:,i) = gwsio;
        
    %% save the results to the output file  
    out_atwsi_file = [output_path '\DA_MTWSA_' time_str '.dat'];
    fid = fopen(out_atwsi_file,'w');
        fwrite(fid,atwsio','float32');
        fwrite(fid,atwsio_std','float32');
    fclose(fid);
    copyfile(hdr_file,[output_path '\DA_MTWSA_' time_str '.hdr']);
    
    out_gwsi_file = [output_path '\DA_MGWSA_' time_str '.dat'];
    fid = fopen(out_gwsi_file,'w');
        fwrite(fid,gwsio','float32');
        fwrite(fid,gwsio_std','float32');
    fclose(fid);
    copyfile(hdr_file,[output_path '\DA_MGWSA_' time_str '.hdr']);
    
    out_gtwsi_file = [output_path '\GR_MTWSA_' time_str '.dat'];
    fid = fopen(out_gtwsi_file,'w');
        fwrite(fid,gtwsio','float32');
        fwrite(fid,gtwsio_std','float32');
    fclose(fid);
    copyfile(hdr_file,[output_path '\GR_MTWSA_' time_str '.hdr']);
    
    out_ttwsi_file = [output_path '\EA_MTWSA_' time_str '.dat'];
    fid = fopen(out_ttwsi_file,'w');
        fwrite(fid,ttwsio','float32');
        fwrite(fid,ttwsio_std','float32');
    fclose(fid);
    copyfile(hdr_file,[output_path '\EA_MTWSA_' time_str '.hdr']);
       
    if daily_files == 1 && apply_gain_factors == 1 && isempty(scale_file)
        all_scal(:,:,i) = scalio;
        out_scali_file = [output_path '\EA_SCALE_' time_str '.dat'];
        fid = fopen(out_scali_file,'w');
            fwrite(fid,scalio','float32');
            fwrite(fid,scalio_std','float32');
        fclose(fid);
        copyfile(hdr_file,[output_path '\EA_SCALE_' time_str '.hdr']);
    end
    
    %% save the mask id into file
    out_mask_file = [output_path '\Canada_JPL_mascon_id_lcc_grid.dat'];
    if ~exist(out_mask_file,'file')
        fid = fopen(out_mask_file,'w');
        fwrite(fid,masconid','float32');
        fclose(fid);
        copyfile(hdr_file,[output_path '\Canada_JPL_mascon_id_lcc_grid.hdr']);
    end
end

%% Output all assimilated monthly results into one file
out_file = [output_path '\DA_ALL_TWSA.dat'];
fid = fopen(out_file,'w');
for k=1:band
    fwrite(fid,all_atws(:,:,k)','float32');
end
fclose(fid);

out_file = [output_path '\GR_ALL_TWSA.dat'];
fid = fopen(out_file,'w');
for k=1:band
    fwrite(fid,all_gtws(:,:,k)','float32');
end
fclose(fid);

out_file = [output_path '\EA_ALL_TWSA.dat'];
fid = fopen(out_file,'w');
for k=1:band
    fwrite(fid,all_ttws(:,:,k)','float32');
end
fclose(fid);

out_file = [output_path '\DA_All_GWSA.dat'];
fid = fopen(out_file,'w');
for k=1:band
    fwrite(fid,all_gws(:,:,k)','float32');
end
fclose(fid);    
out_file = [output_path '\DA_All_STWSA.dat'];
fid = fopen(out_file,'w');
for k=1:band
    fwrite(fid,all_stws(:,:,k)','float32');
end
fclose(fid);

if daily_files == 1 && apply_gain_factors == 1 && isempty(scale_file)
    out_file = [output_path '\DA_ALL_SCALE.dat'];
    fid = fopen(out_file,'w');
    for k=1:band
        fwrite(fid,all_scal(:,:,k)','float32');
    end
    fclose(fid);
end

ttime = toc;

disp(['>>> Elapsed time: ', num2str(ttime),' seconds. <<<'])

msgbox('Completed adjustment successfully!');

logmsg(flog,handles,['>>> Elapsed time: ', num2str(ttime),' seconds. <<<']); 
logmsg(flog,handles,'Completed adjustment successfully!'); 

fclose(flog);
fclose(fstlog);


% --- Executes on button press in pbSyncGRACEandLSMData.
function pbSyncGRACEandLSMData_Callback(hObject, ~, handles)
% hObject    handle to pbSyncGRACEandLSMData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% record the start time
tic;

%% get GRACE TWS input data file name
gtws_file = get(handles.txtGraceTwsFile, 'String');
if isempty(gtws_file)
    msgbox('Specify the GRACE TWS Anomaly NetCDF file.');
    return;
end

%% get the output data ENVI header file
hdr_file = get(handles.txtOutputFile, 'String');
out_file = get(handles.txtOutputFile, 'String');
if isempty(out_file)
    msgbox('Specify the output data ENVI header file.');
    return;
end
%% retrive the output data filder
[output_path, outfname, ext] = fileparts(out_file);

%% get EALCO TWS anomaly data file list
daily_files = get(handles.popFileType, 'Value');
if daily_files == 1
    fname = 'daily';
else
    fname = 'monthly';
end
file_list = get(handles.lstTerresTwsTimeSeries,'String');
if iscell(file_list)
    nf = length(file_list);
else
    nf = 1;
end
if nf<1
    msgbox('Add at least one LSM TWS data file to the LSM TWS Time Series list!');
    return;
end
%% create a time and time string list for the EALCO TWS input file list
time_list = zeros(nf,1);
time_str_list = cell(nf,1);
for i = 1:nf
    if ischar(file_list)
        ttwsi_file = file_list;
    else
        ttwsi_file = char(file_list{i});
    end
    if daily_files == 1
        [time_list(i), time_str_list{i}] = daysFromFileName(ttwsi_file,1);
    else
        [time_list(i), time_str_list{i}] = daysFromFileName(ttwsi_file,3);
    end
end
%% note: the time stamps is in days from Jan. 1, 2002

%% get the meta data parameters and processing options
switch get(handles.popDataType, 'Value')
    case 1      % daily data in int16
        data_type =  'int16';
    case 2      % monthly data in float32
        data_type =  'float32';
end
switch get(handles.popNoDataValue, 'Value')
    case 1      
        no_data =  -32760;
    case 2      
        no_data =  9999;
end

row = str2double(get(handles.txtRow, 'String'));    % row number of the GRACE TWS data
col = str2double(get(handles.txtCol, 'String'));    % col number of the GRACE TWS data
band = str2double(get(handles.txtBand, 'String'));  % the number of months 

%% get the precalculated scale factors file name for Canada LCC grid
apply_gain_factors = get(handles.chkApplyGainFactors, 'Value');
scale_file = get(handles.txtGainFactorFile, 'String');

%% add geo reference coordinates to EALCO TWS
out_lcc_latlon_file = [output_path '\Canada_Lat_Lon_lcc_grid.mat'];
if ~exist(out_lcc_latlon_file,'file')
    lamb = Lambert(...          % Define Canada LCC projection parameters
        [49.0, 77.0], ...       % Latitude of the first and second standard parallel
        0.0, -95, ...           % Central latitude and longitude
        0.0, 0.0, ...           % False easting and northing
        6378137.0, 298.257222100882711243 ...     % Ellipsoid: Semi-major axis and inverse flattening
    );

    lcc_ulx_e = -2600000.0;     
    lcc_uly_n = 10500000.0;
    % lcc_lrx_e = 3100000.0;      
    % lcc_lry_n = 5700000.0;
    lcc_dx = 5000.0;            
    lcc_dy = 5000.0;
    lcc_x = zeros(row, col);   
    lcc_y = zeros(row, col);
    lcc_lat = zeros(row, col);   
    lcc_lon = zeros(row, col);

    for i = 1:row
        for j =1:col
            lcc_x(i,j) = lcc_ulx_e + j*lcc_dx - lcc_dx/2;
            lcc_y(i,j) = lcc_uly_n - i*lcc_dy + lcc_dy/2;
            [lcc_lat(i,j), lcc_lon(i,j)] = lamb.cartesian2geographic(lcc_x(i,j), lcc_y(i,j));
        end
    end

    save(out_lcc_latlon_file,'lcc_x', 'lcc_y', 'lcc_lat', 'lcc_lon');
else
    load(out_lcc_latlon_file);
end

%% determine the geolocation area of CANADA
idx = find(lcc_lon < 0);
lcc_lon(idx) = lcc_lon(idx)+360;
minLat = min(lcc_lat(:));   maxLat = max(lcc_lat(:));
minLon = min(lcc_lon(:));   maxLon = max(lcc_lon(:));

%% read the GRACE TWS anomaly data from NETCDF file
% ncdisp(gtws_file)
gtws = ncread(gtws_file,'lwe_thickness');
gtws_std = ncread(gtws_file,'uncertainty');
% lat = ncread(gtws_file,'lat');
% lon = ncread(gtws_file,'lon');
day = ncread(gtws_file,'time');
lat_bounds = ncread(gtws_file,'lat_bounds');
lon_bounds = ncread(gtws_file,'lon_bounds');
day_bounds = ncread(gtws_file,'time_bounds');

%% determine the GRACE TWS anomaly data index range for CANADA area
[idx_lon,idx_lat,b] = size(gtws);
for i = 1:idx_lat
    if minLat >= lat_bounds(1,i) && minLat < lat_bounds(2,i)
        idx_lat_start = i;
        break;
    end
end
for i = 1:idx_lat
    if maxLat >= lat_bounds(1,i) && maxLat < lat_bounds(2,i)
        idx_lat_end = i;
        break;
    end
end
for i = 1:idx_lon
    if minLon >= lon_bounds(1,i) && minLon < lon_bounds(2,i)
        idx_lon_start = i;
        break;
    end
end
for i = 1:idx_lon
    if maxLon >= lon_bounds(1,i) && maxLon < lon_bounds(2,i)
        idx_lon_end = i;
        break;
    end
end
%% retreive GRACE TWS Anomaly and their uncerntainty data for Canada area
gtws = gtws(idx_lon_start:idx_lon_end,idx_lat_start:idx_lat_end,:);
gtws_std = gtws_std(idx_lon_start:idx_lon_end,idx_lat_start:idx_lat_end,:);
% lat = lat(idx_lat_start:idx_lat_end);
% lon = lon(idx_lon_start:idx_lon_end);
lat_bounds = lat_bounds(:,idx_lat_start:idx_lat_end);
lon_bounds = lon_bounds(:,idx_lon_start:idx_lon_end);


%% open log file to record computation progresses
logfile = [output_path '\PM1_Process_' fname '.log'];
if exist(logfile,'file')
    delete logfile;
end
flog = fopen(logfile,'wt');

logmsg(flog,handles,'Processing started.');

%% initialize the output data container for all monthly datasets
%% the output data is presented in LCC coordinate system
if band > nf
    band = nf;
end

masconid = zeros(row, col);

%% determine the size of GRACE TWS anomaly input data
%% note: the size of gtws is for Canada coverage only
[idx_lon,idx_lat,b] = size(gtws);  % b may be >= band

%% compute or read the baseline values 
out_grace_baseline_file = [output_path '\Canada_GRACE_TWS_anomaly_' num2str(band) '_month_baseline.mat'];
if ~exist(out_grace_baseline_file,'file')
    gbase = zeros(idx_lon,idx_lat);
    for i = 1:idx_lon
        for j = 1:idx_lat
            gbase(i,j) = mean(gtws(i,j,1:band));
        end
    end
    save(out_grace_baseline_file, 'gbase');
else
    load(out_grace_baseline_file);
end

%% process data band by band for each monscon
%% note: each band represnts a monthly dataset
for i = 1:band         
    %% determine the day ranges for the monthly dataset i
    day_bounds_i = day_bounds(:,i);
    day_start = day_bounds_i(1);
    day_end = day_bounds_i(2);
    dayi = day(i);

    logmsg(flog,handles,' ');
    logmsg(flog,handles,['Processing GRACE TWS anomaly for the day: ' num2str(dayi)]);
    logmsg(flog,handles,' ');    
    
    %% read or compute monthly averaged EALCO TWS anomaly data
    if daily_files == 1
        %% retreive the daily EALCO TWS files for the GRACE TWS monthly dataset   
        idx_days = find(time_list>=day_start & time_list<=day_end);       
        if isempty(idx_days)        
            logmsg(flog,handles,['No matched daily LSM TWS files found for the days from ' num2str(day_start) ' to ' num2str(day_end)]);
            continue;
        else   
            %% continue processing the matched GRACE and EALCO data    
            logmsg(flog,handles,[num2str(length(idx_days)) ' daily LSM TWS files found for the days from ' num2str(day_start) ' to ' num2str(day_end)]);
            mfile_list = file_list(idx_days);
            mtime_str_list = time_str_list(idx_days);
            %% compute monthly average of EALCO TWS from the input file list
            if apply_gain_factors == 1 && isempty(scale_file)
                [ttwsi,tstdi,dtwsi] = averageFromInputFiles(mfile_list, row, col, data_type, no_data);
            else
                [ttwsi,tstdi] = averageFromInputFiles(mfile_list, row, col, data_type, no_data);
            end
            idx_middle = fix(length(idx_days)/2+0.5);
            %% determine output file name according to the time
            time_str = char(mtime_str_list(idx_middle));
        end
    else
        % idx_month = find(time_list>day_start+1 & time_list<day_end-1);
        idx_month = find(time_list>dayi-5 & time_list<dayi+5);
        if isempty(idx_month)        
            logmsg(flog,handles,['No matched monthly LSM TWS files found for the days from ' num2str(day_start) ' to ' num2str(day_end)]);
            continue;
        else   
            %% continue processing the matched GRACE and EALCO data    
            logmsg(flog,handles,[num2str(length(idx_month)) ' monthly LSM TWS file found for the days from ' num2str(day_start) ' to ' num2str(day_end)]);
            %% read in monthly average EALCO TWS data directly from file
            if ischar(file_list)
                ttwsi_file = file_list;
            else
                ttwsi_file = char(file_list{idx_month});
            end
            %% read in LSM TWSA data from the data file
            fid = fopen(ttwsi_file);
                ttwsi = (fread(fid,[col, row], data_type)');
                tstdi = (fread(fid,[col, row], data_type)');
            fclose(fid);
            %% determine output file name according to the time
            time_str = char(time_str_list(idx_month));
        end
    end
   
    %% retreive GRACE TWS data for the same month
    gtwsi = gtws(:,:,i)';
    gtwsi_std = gtws_std(:,:,i)';
    %% deduct the base line value
    gtwsi = (gtwsi - gbase');  
    %% convert the unit from cm to mm
    gtwsi = gtwsi*10;
    gtwsi_std = gtwsi_std*10;

    %% initialize the output data after assimilation vy copy the ttws
    gtwsio = ones(row, col)*no_data;
    ttwsio = ones(row, col)*no_data;
    
    gtwsio_std = ones(row, col)*no_data;
    ttwsio_std = ones(row, col)*no_data;
        
    if daily_files == 1 && apply_gain_factors == 1 && isempty(scale_file)
        scalio = ones(row, col)*no_data;
        scalio_std = ones(row, col)*no_data;
    end
    
    %% start a loop to downscale and assimilate GRACE TWS to EALCO TWS
    continue_process_mascon = 1; 
    idx_lat1 = 0; idx_lon1 = 0;
    next_lat = 1;
    while continue_process_mascon > 0 && idx_lat1 <= idx_lat-1
        %% look for a GRACE TWS mascon data and its Lat Lon bounds
        if next_lat == 1
            idx_lat0 = idx_lat1+1; idx_lon0 = idx_lon1 + 1;
            mgtws = gtwsi(idx_lat0, idx_lon0);
            mgstd = gtwsi_std(idx_lat0,idx_lon0);
            for j=idx_lat0:idx_lat
                if gtwsi(j,idx_lon0) ~= mgtws
                    idx_lat1 = j-1;
                    break;
                elseif j == idx_lat
                    idx_lat1 = j;
                end
            end
        else
            idx_lon0 = idx_lon1 + 1;
            mgtws = gtwsi(idx_lat0, idx_lon0);
            mgstd = gtwsi_std(idx_lat0,idx_lon0);
        end
        for j=idx_lon0:idx_lon
            if gtwsi(idx_lat0,j) ~= mgtws
                idx_lon1 = j-1;
                break;
            elseif j == idx_lon
                idx_lon1 = j;
            end
        end
        lat0 = lat_bounds(1,idx_lat0);
        lat1 = lat_bounds(2,idx_lat1);
        lon0 = lon_bounds(1,idx_lon0);
        lon1 = lon_bounds(2,idx_lon1);

        %% check all valuse within the mascon are same (for debug test only)
        G = unique(mgtws); G_std = unique(mgstd);
        ng = length(G); ng_std = length(G_std);
        logmsg(flog,handles,[num2str(ng) ' GRACE TWS Values found for the mascon id ' num2str(continue_process_mascon)]);
        logmsg(flog,handles,[num2str(ng_std) ' GRACE TWS STD found for the mascon id ' num2str(continue_process_mascon)]);
        if ng > 1 
            countg = zeros(ng,1);
            for k=1:ng
                countg(k) = length(find(G == G(k)));
                logmsg(flog,handles,['GRACE TWS Value: ' num2str(G(k)) ' n = ' num2str(countg(k))]);
            end
            [~, idx2] = max(countg);
            G = G(idx2);
            G_std= G_std(idx2);
            logmsg(flog,handles,['Apply GRACE TWS Value ' num2str(G) ' for the mascon id ' num2str(continue_process_mascon)]);
            logmsg(flog,handles,['Apply GRACE TWS STD ' num2str(G_std) ' for the mascon id ' num2str(continue_process_mascon)]);
        end
        
        %% determine the mascon id        
        m_id = continue_process_mascon; 
        %% retreive EALCO TWS within the mascon
        midx = find((lcc_lat>=lat0 & lcc_lat<lat1) & (lcc_lon>=lon0 & lcc_lon<lon1));
        mttws = ttwsi(midx);
        mtstd = tstdi(midx);

        %% remove invalid EALCO TWS
        inval = find(mttws == -32760 | mttws == 9999 | isnan(mttws)); 
        if ~isempty(inval)
            midx(inval) = [];
            mttws(inval) = [];
            mtstd(inval) = [];
        end          
        if isempty(mttws)
            logmsg(flog,handles,' ');
            logmsg(flog,handles,['No valid EALCO TWS found within the mascon ' num2str(continue_process_mascon) ': lat ' num2str(lat0) '-' num2str(lat1) ' lon ' num2str(lon0) '-' num2str(lon1)]);
            %% contune to next mascon
            continue_process_mascon = continue_process_mascon + 1;
        else
            logmsg(flog,handles,' ');
            logmsg(flog,handles,['Found ' num2str(length(mttws)) ' EALCO TWS within the mascon ' num2str(continue_process_mascon) ': lat ' num2str(lat0) '-' num2str(lat1) ' lon ' num2str(lon0) '-' num2str(lon1)]);            

            %% calculate or retreive scale factors
            if daily_files == 1 && apply_gain_factors == 1 && isempty(scale_file)
                [mscal,msstd] = scaleFactorFromEALCO(dtwsi,midx);
                scalio(midx) = mscal;
                scalio_std(midx) = msstd;    
            end

            gtwsio(midx) = ones(length(midx),1)*G;
            gtwsio_std(midx) = ones(length(midx),1)*G_std;
            ttwsio(midx) = mttws;
            ttwsio_std(midx) = mtstd;
            %% assign a valid mascon id 
            masconid(midx) = m_id;
            %% contune to next mascon
            continue_process_mascon = continue_process_mascon + 1;
        end
        if idx_lon1 == idx_lon && idx_lat1 == idx_lat
            continue_process_mascon = 0;
        elseif idx_lon1 == idx_lon && idx_lat1 < idx_lat
            next_lat = 1;
            idx_lon1 = 0;
        else
            idx_lon0 = idx_lon1+1;
            next_lat = 0;
        end                                          
    end
    %% completed processing one dataset (a time stamp of one month)
          
    out_gtwsi_file = [output_path '\GRC_MTWSA_' time_str '.dat'];
    fid = fopen(out_gtwsi_file,'w');
        fwrite(fid,gtwsio','float32');
        fwrite(fid,gtwsio_std','float32');
    fclose(fid);
    copyfile(hdr_file,[output_path '\GRC_MTWSA_' time_str '.hdr']);
    
    out_ttwsi_file = [output_path '\LSM_MTWSA_' time_str '.dat'];
    fid = fopen(out_ttwsi_file,'w');
        fwrite(fid,ttwsio','float32');
        fwrite(fid,ttwsio_std','float32');
    fclose(fid);
    copyfile(hdr_file,[output_path '\LSM_MTWSA_' time_str '.hdr']);
       
    if daily_files == 1 && apply_gain_factors == 1 && isempty(scale_file)
        out_scali_file = [output_path '\LSM_SCALE_' time_str '.dat'];
        fid = fopen(out_scali_file,'w');
            fwrite(fid,scalio','float32');
            fwrite(fid,scalio_std','float32');
        fclose(fid);
        copyfile(hdr_file,[output_path '\LSM_SCALE_' time_str '.hdr']);
    end
    
    %% save the mask id into file
    out_mask_file = [output_path '\Canada_JPL_mascon_id_lcc_grid.dat'];
    if ~exist(out_mask_file,'file')
        fid = fopen(out_mask_file,'w');
        fwrite(fid,masconid','float32');
        fclose(fid);
        copyfile(hdr_file,[output_path '\Canada_JPL_mascon_id_lcc_grid.hdr']);
    end
end

ttime = toc;

disp(['>>> Elapsed time: ', num2str(ttime),' seconds. <<<'])

msgbox('Sync GRACE and LSM TWSA is done!');

logmsg(flog,handles,['>>> Elapsed time: ', num2str(ttime),' seconds. <<<']); 
logmsg(flog,handles,'Sync GRACE and LSM TWSA is done!'); 

fclose(flog);



function txtEALCOBaselineDataFile_Callback(hObject, eventdata, handles)
% hObject    handle to txtEALCOBaselineDataFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtEALCOBaselineDataFile as text
%        str2double(get(hObject,'String')) returns contents of txtEALCOBaselineDataFile as a double


% --- Executes during object creation, after setting all properties.
function txtEALCOBaselineDataFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtEALCOBaselineDataFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowseBaselineDataFile.
function pbBrowseBaselineDataFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbBrowseBaselineDataFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile( ...
    {  '*Ealco_baseline*.*','Specify EALCO Baseline Data File (*Ealco_baseline*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify the EALCO Baseline Data File');
 if  filename > 0  
    set(handles.txtEALCOBaselineDataFile, 'String', [pathname filename]);
 end 
