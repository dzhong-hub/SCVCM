function varargout = TWSA_Monthly_TS_Analysis(varargin)
% TWSA_Monthly_TS_Analysis MATLAB code for TWSA_Monthly_TS_Analysis.fig
%      TWSA_Monthly_TS_Analysis, by itself, creates a new TWSA_Monthly_TS_Analysis or raises the existing
%      singleton*.
%
%      H = TWSA_Monthly_TS_Analysis returns the handle to a new TWSA_Monthly_TS_Analysis or the handle to
%      the existing singleton*.
%
%      TWSA_Monthly_TS_Analysis('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TWSA_Monthly_TS_Analysis.M with the given input arguments.
%
%      TWSA_Monthly_TS_Analysis('Property','Value',...) creates a new TWSA_Monthly_TS_Analysis or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TWSA_Monthly_TS_Analysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TWSA_Monthly_TS_Analysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TWSA_Monthly_TS_Analysis

% Last Modified by GUIDE v2.5 29-Oct-2020 15:00:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TWSA_Monthly_TS_Analysis_OpeningFcn, ...
                   'gui_OutputFcn',  @TWSA_Monthly_TS_Analysis_OutputFcn, ...
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


% --- Executes just before TWSA_Monthly_TS_Analysis is made visible.
function TWSA_Monthly_TS_Analysis_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TWSA_Monthly_TS_Analysis (see VARARGIN)

% Choose default command line output for TWSA_Monthly_TS_Analysis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Add the path to PACE: Principal Analysis by Conditional Expectation
% addpath(genpath('C:\TWSModel\PACEV217\release2.17'));

% UIWAIT makes TWSA_Monthly_TS_Analysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TWSA_Monthly_TS_Analysis_OutputFcn(hObject, eventdata, handles) 
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
    {  '*.hdr','Specify output data folder (*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify output data folder');
 if  filename > 0  
    %set(handles.txtOutputFile, 'String', [pathname filename]);
    set(handles.txtOutputFile, 'String', pathname);
 end 


function txtOutputFile_Callback(~, eventdata, handles)
% hObject    handle to txtOutputFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtOutputFile as text
%        str2double(get(hObject,'String')) returns contents of txtOutputFile as a double


% --- Executes during object creation, after setting all properties.
function txtOutputFile_CreateFcn(hObject, eventdata, ~)
% hObject    handle to txtOutputFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbClose.
function pbClose_Callback(~, ~, ~)
% hObject    handle to pbClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close;

% --- Executes on button press in pbValidateByGridMGMWObs.
function pbValidateByGridMGMWObs_Callback(hObject, eventdata, handles)
% hObject    handle to pbValidateByGridMGMWObs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% record the start time
tic;

row = str2double(get(handles.txtRow, 'String'));    % row number of the GRACE TWS data
col = str2double(get(handles.txtCol, 'String'));    % col number of the GRACE TWS data
band = str2double(get(handles.txtBand, 'String'));  % the number of months

sy = str2double(get(handles.txtSpecificYiel, 'String'));    % the specific yield value
fl = 12*get(handles.popFilterLength, 'Value')+1;            % The filter length


gtws_file = get(handles.txtGraceTwsFile, 'String');
if isempty(gtws_file)
    msgbox('Specify the GRACE TWS Anomaly NetCDF file.');
    return;
end
%% get JPL mascon id file name for Canada LCC grid
mask_file = get(handles.txtGraceTwsIdMaskFile, 'String');
if isempty(mask_file)
    msgbox('Specify the JPL mascon id file for Canada LCC grid.');
    return;
end
%% get Canada LCC grid Lat Lon position file
lcc_latlon_file = get(handles.txtLCCGridLatLonFile, 'String');
if isempty(lcc_latlon_file)
    msgbox('Specify the Canada LCC grid lat lon position file.');
    return;
end

%% get GWMWL observation input folder
gmwl_in_dir = get(handles.txtGWTSInputFolder, 'String');
if isempty(gmwl_in_dir)
    msgbox('Specify the GMWL observations input folder.');
    return;
end
%% get GWMWL lat lon position file
gmwl_pos_file = get(handles.txtGWLatLonFile, 'String');
if isempty(gmwl_pos_file)
    msgbox('Specify the GMWL Lat Lon position file.');
    return;
end

%% get the output data ENVI header file
output_path = get(handles.txtOutputFile, 'String');
if isempty(output_path)
    msgbox('Specify the output data folder.');
    return;
end

%% read all original GRACE TWSA data
% get all original grace twsa file names
gtwsa_folder = get(handles.txtAllGraceTwsFile, 'String');
if isempty(gtwsa_folder)
    msgbox('Specify the original GRACE TWSA data file folder.');
    return;
end
allfnames = dir(gtwsa_folder);
gtwsa_fnames = {allfnames(~[allfnames.isdir]).name};
n_gtwsa = length(gtwsa_fnames);
gtws = zeros(row,col,n_gtwsa);
time_list = zeros(n_gtwsa,1);
time_str_list = cell(n_gtwsa,1);
for i = 1:n_gtwsa
    gtwsai_file = [gtwsa_folder '\' gtwsa_fnames{i}];
    [time_list(i), time_str_list{i}] = daysFromFileName(gtwsai_file,3);
    fid = fopen(gtwsai_file);
        gtws(:,:,i) = (fread(fid,[col, row],'float32')');
    fclose(fid);
end  

%% read all EALCO TWSA data
% get all EALCO twsa file names
ttwsa_folder = get(handles.txtAllEalcoTwsFile, 'String');
if isempty(ttwsa_folder)
    msgbox('Specify the EALCO TWSA data file folder.');
    return;
end
allfnames = dir(ttwsa_folder);
ttwsa_fnames = {allfnames(~[allfnames.isdir]).name};
n_ttwsa = length(ttwsa_fnames);
ttws = ones(row,col,n_gtwsa)*-32760;
file_read_count = 0;
for i = 1:n_ttwsa
    ttwsai_file = [ttwsa_folder '\' ttwsa_fnames{i}];
    % sync the time stamps
    [time_day, ttws_time_str] = daysFromFileName(ttwsai_file,3);
    [tf, idx] = ismember(ttws_time_str,time_str_list);
    if tf
        file_read_count = file_read_count + 1;
        fid = fopen(ttwsai_file);
        ttws(:,:,idx) = (fread(fid,[col, row],'float32')');
        fclose(fid);
    end
end
if file_read_count ~= n_gtwsa
    msgbox('Error: the GTWSA file number ~= the read TTWSA file number! Program stopped!');
    return;
end

%% read all simply scaled GRACE TWS data
stwsa_file = get(handles.txtAllScaledTwsFile, 'String');
if isempty(stwsa_file)
    msgbox('Specify the simply scaled TWSA data file.');
    return;
end
stws = zeros(row,col,n_gtwsa);
fid = fopen(stwsa_file);
for i = 1:n_gtwsa
    stws(:,:,i) = (fread(fid,[col, row],'float32')');  
end 
fclose(fid);

% stwsa_folder = get(handles.txtAllScaledTwsFile, 'String');
% if isempty(stwsa_folder)
%     msgbox('Specify the simply scaled TWSA data file folder.');
%     return;
% end
% allfnames = dir(stwsa_folder);
% stwsa_fnames = {allfnames(~[allfnames.isdir]).name};
% n_stwsa = length(stwsa_fnames);
% stws = ones(row,col,n_gtwsa)*-32760;
% file_read_count = 0;
% for i = 1:n_stwsa
%     stwsai_file = [stwsa_folder '\' stwsa_fnames{i}];
%     % sync the time stamps
%     [time_day, stws_time_str] = daysFromFileName(stwsai_file,3);
%     [tf, idx] = ismember(stws_time_str,time_str_list);
%     if tf
%         file_read_count = file_read_count + 1;
%         fid = fopen(stwsai_file);
%         stws(:,:,idx) = (fread(fid,[col, row],'float32')');
%         fclose(fid);
%     end
% end
% if file_read_count ~= n_gtwsa
%     msgbox('Error: the GTWSA file number ~= the read STWSA file number! Program stopped!');
%     return;
% end

%% read all spatially downscaled TWSA data 1
% get all saptially downscaled and assimilated TWSA file names
sdtwsa_folder = get(handles.txtAllFinalTwsFile1, 'String');
if isempty(sdtwsa_folder)
    msgbox('Specify the saptially downscaled TWSA file folder.');
    return;
end
allfnames = dir(sdtwsa_folder);
sdtwsa_fnames = {allfnames(~[allfnames.isdir]).name};
% idx_hdr = contains('.hdr',sdtwsa_fnames);
% sdtwsa_fnames(idx_hdr)=[];
n_sdtwsa = length(sdtwsa_fnames);
atws1 = ones(row,col,n_gtwsa)*-32760;
file_read_count = 0;
for i = 1:n_sdtwsa
    sdtwsai_file = [sdtwsa_folder '\' sdtwsa_fnames{i}];
    % sync the time stamps
    [time_day, atws_time_str] = daysFromFileName(sdtwsai_file,3);
    [tf, idx] = ismember(atws_time_str,time_str_list);
    if tf
        file_read_count = file_read_count + 1;
        fid = fopen(sdtwsai_file);
        atws1(:,:,idx) = (fread(fid,[col, row],'float32')');
        fclose(fid);
    end
end
if file_read_count ~= n_gtwsa
    msgbox('Error: the GTWSA file number ~= the read SDTWSA file number! Program stopped!');
    return;
end


%% read all estimated GSWSA data 1
gswsa_folder = get(handles.txtAllGswsFile1, 'String');
if isempty(sdtwsa_folder)
    msgbox('Specify the saptially downscaled GSWSA file folder.');
    return;
end
allfnames = dir(gswsa_folder);
gswsa_fnames = {allfnames(~[allfnames.isdir]).name};
n_gswsa = length(gswsa_fnames);
gsws1 = ones(row,col,n_gtwsa)*-32760;
file_read_count = 0;
for i = 1:n_gswsa
    gswsai_file = [gswsa_folder '\' gswsa_fnames{i}];
    % sync the time stamps
    [time_day, gswsa_time_str] = daysFromFileName(gswsai_file,3);
    [tf, idx] = ismember(gswsa_time_str,time_str_list);
    if tf
        file_read_count = file_read_count + 1;
        fid = fopen(gswsai_file);
        gsws1(:,:,idx) = (fread(fid,[col, row],'float32')');
        fclose(fid);
    end
end
if file_read_count ~= n_gtwsa
    msgbox('Error: the GTWSA file number ~= the read GSWSA1 file number! Program stopped!');
    return;
end

%% read all spatially downscaled TWSA data 2
% get all saptially downscaled and assimilated TWSA file names
sdtwsa_folder = get(handles.txtAllFinalTwsFile2, 'String');
if isempty(sdtwsa_folder)
    msgbox('Specify the saptially downscaled TWSA file folder 2.');
    return;
end
allfnames = dir(sdtwsa_folder);
sdtwsa_fnames = {allfnames(~[allfnames.isdir]).name};
n_sdtwsa = length(sdtwsa_fnames);
atws2 = ones(row,col,n_gtwsa)*-32760;
file_read_count = 0;
for i = 1:n_sdtwsa
    sdtwsai_file = [sdtwsa_folder '\' sdtwsa_fnames{i}];
    % sync the time stamps
    [time_day, atws_time_str] = daysFromFileName(sdtwsai_file,3);
    [tf, idx] = ismember(atws_time_str,time_str_list);
    if tf
        file_read_count = file_read_count + 1;
        fid = fopen(sdtwsai_file);
        atws2(:,:,idx) = (fread(fid,[col, row],'float32')');
        fclose(fid);
    end
end
if file_read_count ~= n_gtwsa
    msgbox('Error: the GTWSA file number ~= the read SDTWSA 2 file number! Program stopped!');
    return;
end


%% read all estimated GSWSA data 2
gswsa_folder = get(handles.txtAllGswsFile2, 'String');
if isempty(gswsa_folder)
    msgbox('Specify the saptially downscaled GSWSA file folder.');
    return;
end
allfnames = dir(gswsa_folder);
gswsa_fnames = {allfnames(~[allfnames.isdir]).name};
n_gswsa = length(gswsa_fnames);
gsws2 = ones(row,col,n_gtwsa)*-32760;
file_read_count = 0;
for i = 1:n_gswsa
    gswsai_file = [gswsa_folder '\' gswsa_fnames{i}];
    % sync the time stamps
    [time_day, gswsa_time_str] = daysFromFileName(gswsai_file,3);
    [tf, idx] = ismember(gswsa_time_str,time_str_list);
    if tf
        file_read_count = file_read_count + 1;
        fid = fopen(gswsai_file);
        gsws2(:,:,idx) = (fread(fid,[col, row],'float32')');
        fclose(fid);
    end
end
if file_read_count ~= n_gtwsa
    msgbox('Error: the GTWSA file number ~= the read GSWSA 2 file number! Program stopped!');
    return;
end

%% read in the pixel id mask input data
fid = fopen(mask_file);
    mask = (fread(fid,[col, row],'float32')');
fclose(fid);

%% load lcc grid lat lon positions
% save(out_lcc_latlon_file,'lcc_x', 'lcc_y', 'lcc_lat', 'lcc_lon');
load(lcc_latlon_file);
lamb = Lambert(...          % Define Canada LCC projection parameters
    [49.0, 77.0], ...       % Latitude of the first and second standard parallel
    0.0, -95, ...           % Central latitude and longitude
    0.0, 0.0, ...           % False easting and northing
    6378137.0, 298.257222100882711243 ...     % Ellipsoid: Semi-major axis and inverse flattening
);
lcc_dx = 5000.0;            
lcc_dy = 5000.0;

%% specify mascon ids to be processed
%maskvalues = [273,274,307,308,341]';
maskvalues = [273,274,308]';
%maskvalues = [162,163,203,204,205,206,273,274,307,308,341]';
%maskvalues = [234,235,236,271,272,305,306,307,339,340]';

st_logfile = [output_path '\Grid_Statistics.log'];
if exist(st_logfile,'file')
    delete st_logfile;
end
fstlog = fopen(st_logfile,'wt');

%% get groundwater well level time series file names 
allfnames = dir(gmwl_in_dir);
fnames = {allfnames(~[allfnames.isdir]).name};
nf = length(fnames);

%% read in lat lon positions of groundwater wells
%[wid,lat,lon,y1,y2,ym] = readvars(gmwl_pos_file);
[wid,lat,lon] = readvars(gmwl_pos_file);
nwid = length(wid);    

%% convert the well id from double to string
if isa(wid,'double')
    wid = num2str(wid);
end

fprintf(fstlog,'%s\n', ' ');
fprintf(fstlog,'%s\n','Temporal correlation comparison results for the estimated and measured GSWS');
fprintf(fstlog,'%s\n','           WellID           Lat         Lon          Mascon     Pixeli      Pixelj       N(month)  b1(yield)    RMSE1       RMSE2       RMSE3       RMSE4         R1          R2          R3          R4        RMSE1       RMSE2       RMSE3       RMSE4         R1          R2          R3          R4');

%% process data well by well
QI = zeros(nwid,17);
gws0 = ones(nwid, band)*-32760;
gws1 = ones(nwid, band)*-32760;
gws2 = ones(nwid, band)*-32760;
gws3 = ones(nwid, band)*-32760;
gws4 = ones(nwid, band)*-32760;

gws0_hp = ones(nwid, band)*-32760;
gws1_hp = ones(nwid, band)*-32760;
gws2_hp = ones(nwid, band)*-32760;
gws3_hp = ones(nwid, band)*-32760;
gws4_hp = ones(nwid, band)*-32760;

for i =1:nwid
    find_wid = 0;
    for j = 1:nf
        if contains(fnames{j},wid(i))
            find_wid = j;
            break;
        end
    end
    if find_wid > 0       
        [we, wn] = lamb.geographic2cartesian(lat(i), lon(i));
        pixi = 0; pixj = 0;
        %% find the pixel location
        for ii = 1: row
            for jj =1:col
                x1 = lcc_x(ii,jj) - lcc_dx/2; x2 = lcc_x(ii,jj) + lcc_dx/2;
                y1 = lcc_y(ii,jj) - lcc_dy/2; y2 = lcc_y(ii,jj) + lcc_dy/2;
                if (we >= x1 && we < x2) && (wn >= y1 && wn < y2)
                    pixi = ii; pixj = jj;
                    break;
                end
            end
            if pixi > 0 && pixj > 0
                break;
            end
        end
        if pixi>0 && pixj>0
            %% get estimated gsws
            mascon_id = mask(pixi,pixj);            
            if ismember(mascon_id,maskvalues)
                g = gsws1(pixi,pixj,:); % estimated gsws with option 1
                y = gsws2(pixi,pixj,:); % estimated gsws with option 2
                z = gtws(pixi,pixj,:);  % grace tws 
                x = ttws(pixi,pixj,:);  % ealco tws
                s = stws(pixi,pixj,:);  % scaled tws           
                g = g(:); z = z(:); x = x(:); y = y(:); s = s(:);

                csvfn = [gmwl_in_dir '\' fnames{find_wid}];
                [t0,w] = readvars(csvfn);
                if length(g) == length(w)                              
                    %% calculate the temporal correlations
                    idxg = find(g ~= -32760);
                    idxw = find(w ~= -32760);
                    idx = intersect(idxg, idxw);                
                    if ~isempty(idx) 
                        g = g(idx); % gsws 1
                        y = y(idx); % gsws 2
                        z = z(idx); % grace tws
                        x = x(idx); % ealco tws
                        s = s(idx); % scaled tws
                        
                        w = w(idx)*1000; % observed gmw height changes
                        
                        t = t0(idx);
                        n = length(idx);
                                                
                        %% deduct the mean values
                        g = g - mean(g);
                        y = y - mean(y);
                        z = z - mean(z);
                        x = x - mean(x);
                        s = s - mean(s);
                        
                        %% calculate GWSA
                        g0 = w - mean(w); % observed gmw height changes
                        g1 = z - x;       % no downscaled gwsa
                        g2 = s - x;       % simple downscaled gwsa
                        g3 = g;           % scvcm estimated gwsa 1
                        g4 = y;           % scvcm estimated gwsa 2
                        
                        %% apply HP filter
%                         g0_hp = hpfilter2(g0,14400);
%                         g1_hp = hpfilter2(g1,14400);
%                         g2_hp = hpfilter2(g2,14400);
%                         g3_hp = hpfilter2(g3,14400);
%                         g4_hp = hpfilter2(g4,14400);

                        B = 1/fl*ones(fl,1);
                        g0_hp = filter(B,1,g0);
                        g1_hp = filter(B,1,g1);
                        g2_hp = filter(B,1,g2);
                        g3_hp = filter(B,1,g3);
                        g4_hp = filter(B,1,g4);

                        %% convert gw water level height to GWS
                        switch sy
                            case 1
                               sv = g0_hp\g1_hp; 
                            case 2
                               sv = g0_hp\g2_hp;
                            case 3
                                %sv = g0_hp\g3_hp;
                                sv = g0\g3;
                            case 4
                                sv = g0_hp\g4_hp;
                            otherwise
                                % use specified value
                                sv = sy;
                        end 
                        g0 = sv*g0;
                        g0_hp = sv*g0_hp;

                        %% plot out the unfiltered and filtered data for comparison                                               
                        widstr = wid(i); 
                        if  contains(widstr,'234745')||contains(widstr,'9649')||contains(widstr,'19128')
                        % if  contains(widstr,'31803')||contains(widstr,'234747')||contains(widstr,'32791') 
% %                             figure;
% %                             hold on;
% %                             plot(t,g0,'r-o',...
% %                                 'LineWidth',1,...
% %                                 'MarkerSize',3,...
% %                                 'MarkerEdgeColor','r',...
% %                                 'MarkerFaceColor',[0.5,0.5,0.5]);
% %                             plot(t,g1,'c-o',...
% %                                 'LineWidth',1,...
% %                                 'MarkerSize',3,...
% %                                 'MarkerEdgeColor','c',...
% %                                 'MarkerFaceColor',[0.5,0.5,0.5]);
% %                             plot(t,g2,'b-*',...
% %                                 'LineWidth',1,...
% %                                 'MarkerSize',3,...
% %                                 'MarkerEdgeColor','b',...
% %                                 'MarkerFaceColor',[0.5,0.5,0.5]);
% %                             plot(t,g3,'g-*',...
% %                                 'LineWidth',1,...
% %                                 'MarkerSize',3,...
% %                                 'MarkerEdgeColor','g',...
% %                                 'MarkerFaceColor',[0.5,0.5,0.5]);
% % %                             plot(t,g4,'b-*',...
% % %                                 'LineWidth',1,...
% % %                                 'MarkerSize',3,...
% % %                                 'MarkerEdgeColor','b',...
% % %                                 'MarkerFaceColor',[0.5,0.5,0.5]);
% %                             %title([wid(i) ' GWSA Comparison before HP Filter']);
% %                             title(['Well ID: ' wid(i)]);
% %                             xlabel('Time in days from 2002-01-01')
% %                             ylabel('GWSA (mm EWH)') 
% %                             legend('OGWSA','EGWSA1','EGWSA2','EGWSA3')
% %                             
% %                             figure;
% %                             hold on;
% %                             plot(t,g0_hp,'r-o',...
% %                                 'LineWidth',1,...
% %                                 'MarkerSize',3,...
% %                                 'MarkerEdgeColor','r',...
% %                                 'MarkerFaceColor',[0.5,0.5,0.5]);
% %                             plot(t,g1_hp,'c-o',...
% %                                 'LineWidth',1,...
% %                                 'MarkerSize',3,...
% %                                 'MarkerEdgeColor','c',...
% %                                 'MarkerFaceColor',[0.5,0.5,0.5]);
% %                             plot(t,g2_hp,'b-*',...
% %                                 'LineWidth',1,...
% %                                 'MarkerSize',3,...
% %                                 'MarkerEdgeColor','b',...
% %                                 'MarkerFaceColor',[0.5,0.5,0.5]);
% %                             plot(t,g3_hp,'g-*',...
% %                                 'LineWidth',1,...
% %                                 'MarkerSize',3,...
% %                                 'MarkerEdgeColor','g',...
% %                                 'MarkerFaceColor',[0.5,0.5,0.5]);
% % %                             plot(t,g4_hp,'b-*',...
% % %                                 'LineWidth',1,...
% % %                                 'MarkerSize',3,...
% % %                                 'MarkerEdgeColor','b',...
% % %                                 'MarkerFaceColor',[0.5,0.5,0.5]);
% %                             title(['Well ID: ' wid(i)]);
% %                             xlabel('Time in days from 2002-01-01')
% %                             ylabel('Smoothed GWSA (mm EWH)') 
% %                             legend('OGWSA','EGWSA1','EGWSA2','EGWSA3')
                            
                            %% save GMW data for daily trend changes plots
                            if  contains(widstr,'234745')
                                well_num = '234745';
                            elseif contains(widstr,'9649')
                                well_num = '9649';
                            elseif contains(widstr,'19128')
                                well_num = '19128';
                            end    
                            used_well_file = [output_path '\Well_' well_num '_in_mascon_' num2str(mascon_id) '_monthly_trend_data.txt'];
                            if exist(used_well_file,'file')
                                delete used_well_file;
                            end
                            fwell = fopen(used_well_file,'wt');
                            fprintf(fwell,'%s\n','             Id       Date        TimeDays      EGWSA1      EGWSA2      EGWSA3       OGWSA      EGWSA1      EGWSA2      EGWSA3      OGWSA');
                            
                            for k = 1:length(g0)
                                fprintf(fwell,'%+16s', char(wid(i)));
                                fprintf(fwell,'%+14s', convertDaysToDateTime(t(k)));
                                fprintf(fwell,'%+12s', num2str(t(k),'% 12.1f'));                                
                                fprintf(fwell,'%+12s', num2str(g1(k),'% 12.6f'));
                                fprintf(fwell,'%+12s', num2str(g2(k),'% 12.6f'));
                                fprintf(fwell,'%+12s', num2str(g3(k),'% 12.6f'));
                                fprintf(fwell,'%+12s', num2str(g0(k),'% 12.6f'));                                
                                fprintf(fwell,'%+12s', num2str(g1_hp(k),'% 12.6f'));
                                fprintf(fwell,'%+12s', num2str(g2_hp(k),'% 12.6f'));
                                fprintf(fwell,'%+12s', num2str(g3_hp(k),'% 12.6f'));
                                fprintf(fwell,'%+12s\n', num2str(g0_hp(k),'% 12.6f'));
                            end
                            fclose(fwell);
                        end
                        

                        %% keep a copy of the time series for spatial correlation comparison
                        gws0(i,idx) = g0';
                        gws1(i,idx) = g1';
                        gws2(i,idx) = g2';
                        gws3(i,idx) = g3';
                        gws4(i,idx) = g4';
                       
                        gws0_hp(i,idx) = g0_hp';
                        gws1_hp(i,idx) = g1_hp';
                        gws2_hp(i,idx) = g2_hp';
                        gws3_hp(i,idx) = g3_hp';
                        gws4_hp(i,idx) = g4_hp';
                        
                        QI(i,1) = sv;                        
                        QI(i,2) = sqrt(mean((g1-g0).^2));
                        QI(i,3) = sqrt(mean((g2-g0).^2));
                        QI(i,4) = sqrt(mean((g3-g0).^2));
                        QI(i,5) = sqrt(mean((g4-g0).^2));                      
                        QI(i,6) = corr(g1,g0); 
                        QI(i,7) = corr(g2,g0);
                        QI(i,8) = corr(g3,g0);
                        QI(i,9) = corr(g4,g0);
                        QI(i,10) = sqrt(mean((g1_hp-g0_hp).^2));
                        QI(i,11) = sqrt(mean((g2_hp-g0_hp).^2));
                        QI(i,12) = sqrt(mean((g3_hp-g0_hp).^2));
                        QI(i,13) = sqrt(mean((g4_hp-g0_hp).^2));                      
                        QI(i,14) = corr(g1_hp,g0_hp); 
                        QI(i,15) = corr(g2_hp,g0_hp);
                        QI(i,16) = corr(g3_hp,g0_hp);
                        QI(i,17) = corr(g4_hp,g0_hp);
                        
                        fprintf(fstlog,'%+6s', num2str(i,'% 6d'));
                        fprintf(fstlog,'%+16s', char(wid(i)));

                        fprintf(fstlog,'%+12s', num2str(lat(i),'% 12.6f'));
                        fprintf(fstlog,'%+12s', num2str(lon(i),'% 12.6f'));
                        fprintf(fstlog,'%+12s', num2str(mascon_id,'% 8d'));
                        fprintf(fstlog,'%+12s', num2str(pixi,'% 8d'));
                        fprintf(fstlog,'%+12s', num2str(pixj,'% 8d'));
                        fprintf(fstlog,'%+12s', num2str(n,'% 12d'));
                        for j=1:17
                            fprintf(fstlog,'%+12s', num2str(QI(i,j),'% 12.4f'));
                        end
                        fprintf(fstlog,'%s\n', ' ');
                        
                    end
                end
            end
        end
    end
end

%% print out mean values
fprintf(fstlog,'%s\n', ' ');
fprintf(fstlog,'%+94s', '   Mean Value  ');
for j=1:17
    fprintf(fstlog,'%+12s', num2str(mean(QI(:,j)),'% 12.4f'));
end
fprintf(fstlog,'%s\n', ' ');

%% plot temporal RMSE comparison
figure;
hold on;
plot(QI(:,2),'r-o',...
    'LineWidth',1,...
    'MarkerSize',3,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
plot(QI(:,3),'c-o',...
    'LineWidth',1,...
    'MarkerSize',3,...
    'MarkerEdgeColor','c',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
plot(QI(:,4),'g-o',...
    'LineWidth',1,...
    'MarkerSize',3,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
% plot(QI(:,5),'b-o',...
%     'LineWidth',1,...
%     'MarkerSize',3,...
%     'MarkerEdgeColor','b',...
%     'MarkerFaceColor',[0.5,0.5,0.5]);
title('Temporal RMSE Comparison');
xlabel('Well ID')
ylabel('RMSE in EWH [mm]') 
legend('GWSA1','GWSA2','GWSA3','GWSA4')

%% plot temporal Peason's correlation comparison
figure;
hold on;
plot(QI(:,6),'r-o',...
    'LineWidth',1,...
    'MarkerSize',3,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
plot(QI(:,7),'c-o',...
    'LineWidth',1,...
    'MarkerSize',3,...
    'MarkerEdgeColor','c',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
plot(QI(:,8),'g-o',...
    'LineWidth',1,...
    'MarkerSize',3,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
% plot(QI(:,9),'b-o',...
%     'LineWidth',1,...
%     'MarkerSize',3,...
%     'MarkerEdgeColor','b',...
%     'MarkerFaceColor',[0.5,0.5,0.5]);
title('Temporal Pearson Correlation Comparison');
xlabel('Well ID')
ylabel('Pearson Correlation') 
legend('GWSA1','GWSA2','GWSA3','GWSA4')

%% scatter plot
x = gws0(:); 
y1 = gws1(:); 
y2 = gws2(:); 
y3 = gws3(:);  
y4 = gws4(:); 

% x_hp = gws0_hp(:);
% y1_hp = gws1_hp(:);
% y2_hp = gws2_hp(:);
% y3_hp = gws3_hp(:); 
% y4_hp = gws4_hp(:); 
    
idxg = find(x ~= -32760);
idxw = find(y1 ~= -32760);
idx = intersect(idxg, idxw);                
if ~isempty(idx)
    x = x(idx); 
    y1 = y1(idx);
    y2 = y2(idx); 
    y3 = y3(idx); 
%     y4 = y4(idx);
% 
%     x_hp = x_hp(idx); 
%     y1_hp = y1_hp(idx);
%     y2_hp = y2_hp(idx); 
%     y3_hp = y3_hp(idx);
%     y4_hp = y4_hp(idx);
end

h1 = figure;
set(h1,'Position',[300 300 1200 320]);
y = y1;
subplot(1,3,1)
plotScaterPoints(x, y, 'OGWSA (mm EWH)', 'EGWSA1 (mm EWH)');
y = y2;
subplot(1,3,2)
plotScaterPoints(x, y, 'OGWSA (mm EWH)', 'EGWSA2 (mm EWH)');
y = y3;
subplot(1,3,3)
plotScaterPoints(x, y, 'OGWSA (mm EWH)', 'EGWSA3 (mm EWH)');

% h1 = figure;
% set(h1,'Position',[300 300 1200 320]);
% x = x_hp;
% y = y1_hp;
% subplot(1,3,1)
% plotScaterPoints(x, y, 'OGWSA (mm EWH)', 'EGWSA1 (mm EWH)');
% y = y2_hp;
% subplot(1,3,2)
% plotScaterPoints(x, y, 'OGWSA (mm EWH)', 'EGWSA2 (mm EWH)');
% y = y3_hp;
% subplot(1,3,3)
% plotScaterPoints(x, y, 'OGWSA (mm EWH)', 'EGWSA3 (mm EWH)');

%% calculate the spatial correlations

fprintf(fstlog,'%s\n', ' ');
fprintf(fstlog,'%s\n','Spatial correlation comparison results for the estimated and measured GSWS');
fprintf(fstlog,'%s\n','    Id       Date        TimeDays     N(wells)    RMSE1       RMSE2       RMSE3       RMSE4         R1          R2          R3          R4       RMSE1       RMSE2       RMSE3       RMSE4         R1          R2          R3          R4');

QIs = zeros(band,17);
for i = 1:band
    g0 = gws0(:,i); 
    g1 = gws1(:,i); 
    g2 = gws2(:,i); 
    g3 = gws3(:,i);  
    g4 = gws4(:,i); 
        
    g0_hp = gws0_hp(:,i);
    g1_hp = gws1_hp(:,i);
    g2_hp = gws2_hp(:,i);
    g3_hp = gws3_hp(:,i); 
    g4_hp = gws4_hp(:,i); 
    
    idxg = find(g0 ~= -32760);
    idxw = find(g1 ~= -32760);
    idx = intersect(idxg, idxw);                
    if ~isempty(idx)
        g0 = g0(idx); 
        g1 = g1(idx);
        g2 = g2(idx); 
        g3 = g3(idx); 
        g4 = g4(idx);
        
        g0_hp = g0_hp(idx); 
        g1_hp = g1_hp(idx);
        g2_hp = g2_hp(idx); 
        g3_hp = g3_hp(idx);
        g4_hp = g4_hp(idx);
        
        QIs(i,1) = 1;                        
        QIs(i,2) = sqrt(mean((g1-g0).^2));
        QIs(i,3) = sqrt(mean((g2-g0).^2));
        QIs(i,4) = sqrt(mean((g3-g0).^2));
        QIs(i,5) = sqrt(mean((g4-g0).^2));
        QIs(i,6) = corr(g1,g0); 
        QIs(i,7) = corr(g2,g0);
        QIs(i,8) = corr(g3,g0);
        QIs(i,9) = corr(g4,g0);
        QIs(i,10) = sqrt(mean((g1_hp-g0_hp).^2));
        QIs(i,11) = sqrt(mean((g2_hp-g0_hp).^2));
        QIs(i,12) = sqrt(mean((g3_hp-g0_hp).^2));
        QIs(i,13) = sqrt(mean((g4_hp-g0_hp).^2));
        QIs(i,14) = corr(g1_hp,g0_hp); 
        QIs(i,15) = corr(g2_hp,g0_hp);
        QIs(i,16) = corr(g3_hp,g0_hp);
        QIs(i,17) = corr(g4_hp,g0_hp);
        
        n = length(idx);
    else
        n = 0;
    end
    
    fprintf(fstlog,'%+6s', num2str(i,'% 6d'));
    fprintf(fstlog,'%+14s', convertDaysToDateTime(t0(i)));
    fprintf(fstlog,'%+12s', num2str(t0(i),'% 12.1f'));
    fprintf(fstlog,'%+12s', num2str(n,'% 12d'));
    for j=2:17
        fprintf(fstlog,'%+12s', num2str(QIs(i,j),'% 12.4f'));
    end
    fprintf(fstlog,'%s\n', ' ');
end

%% print out mean values
fprintf(fstlog,'%s\n', ' ');
fprintf(fstlog,'%+44s', '   Mean Value  ');
for j=2:17
    fprintf(fstlog,'%+12s', num2str(mean(QIs(:,j)),'% 12.4f'));
end
fprintf(fstlog,'%s\n', ' ');

%% Plot Spatial RMSE Comparison
figure;
hold on;
plot(QIs(:,2),'r-o',...
    'LineWidth',1,...
    'MarkerSize',3,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
plot(QIs(:,3),'c-o',...
    'LineWidth',1,...
    'MarkerSize',3,...
    'MarkerEdgeColor','c',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
plot(QIs(:,4),'g-o',...
    'LineWidth',1,...
    'MarkerSize',3,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
% plot(QIs(:,5),'b-o',...
%     'LineWidth',1,...
%     'MarkerSize',3,...
%     'MarkerEdgeColor','b',...
%     'MarkerFaceColor',[0.5,0.5,0.5]);
title('Spatial RMSE Comparison');
xlabel('Time Epoch in Month')
ylabel('RMSE in EWH [mm]') 
legend('GWSA1','GWSA2','GWSA3','GWSA4')

%% Plot Spatial Pearson Correlation Comparison
figure;
hold on;
plot(QIs(:,6),'r-o',...
    'LineWidth',1,...
    'MarkerSize',3,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
plot(QIs(:,7),'c-o',...
    'LineWidth',1,...
    'MarkerSize',3,...
    'MarkerEdgeColor','c',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
plot(QIs(:,8),'g-o',...
    'LineWidth',1,...
    'MarkerSize',3,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
% plot(QIs(:,9),'b-o',...
%     'LineWidth',1,...
%     'MarkerSize',3,...
%     'MarkerEdgeColor','b',...
%     'MarkerFaceColor',[0.5,0.5,0.5]);
title('Spatial Pearson Correlation Comparison');
xlabel('Time Epoch in Month')
ylabel('Pearson Correlation') 
legend('GWSA1','GWSA2','GWSA3','GWSA4')

ttime = toc;

disp(['>>> Elapsed time: ', num2str(ttime),' seconds. <<<'])

msgbox('Completed TWS data analysis successfully!');

fclose(fstlog);



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



function txtTerresTwsStd_Callback(hObject, eventdata, handles)
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
function popUnitTypeGraceTwsStd_CreateFcn(hObject, ~, handles)
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
function pbBrowseGraceTwsFile_Callback(~, eventdata, handles)
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


function txtGraceTwsIdMaskFile_Callback(hObject, ~, handles)
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
    {  '*tws*.bin','Input image files (*tws*.bin)'; ...
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
function lstTerresTwsTimeSeries_Callback(hObject, eventdata, handles)
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
    {  '*tws*.bin','Input GRACE TWS files (*tws*.bin)'; ...
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




function txtLCCGridLatLonFile_Callback(hObject, eventdata, handles)
% hObject    handle to txtLCCGridLatLonFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtLCCGridLatLonFile as text
%        str2double(get(hObject,'String')) returns contents of txtLCCGridLatLonFile as a double


% --- Executes during object creation, after setting all properties.
function txtLCCGridLatLonFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtLCCGridLatLonFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowseLCCGridLatLonFile.
function pbBrowseLCCGridLatLonFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbBrowseLCCGridLatLonFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {  '*Lat_Lon*.*','Specify Canada LCC Lat Lon file (*Lat_Lon*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify Canada LCC Lat Lon file');
 if  filename > 0  
    set(handles.txtLCCGridLatLonFile, 'String', [pathname filename]);
 end 


function txtAllGraceTwsFile_Callback(hObject, eventdata, handles)
% hObject    handle to txtAllGraceTwsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtAllGraceTwsFile as text
%        str2double(get(hObject,'String')) returns contents of txtAllGraceTwsFile as a double


% --- Executes during object creation, after setting all properties.
function txtAllGraceTwsFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtAllGraceTwsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowseAllGraceTwsFile.
function pbBrowseAllGraceTwsFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbBrowseAllGraceTwsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {  'GR_All*.*','Specify a GRACE TWSA LCC Grid Time Series file (*GRACE_TWSA*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify a GRACE TWSA LCC Grid Time Series file');
 if  filename > 0  
    set(handles.txtAllGraceTwsFile, 'String', pathname);
    %set(handles.txtAllGraceTwsFile, 'String', [pathname filename]);
 end 


function txtAllEalcoTwsFile_Callback(hObject, eventdata, handles)
% hObject    handle to txtAllEalcoTwsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtAllEalcoTwsFile as text
%        str2double(get(hObject,'String')) returns contents of txtAllEalcoTwsFile as a double


% --- Executes during object creation, after setting all properties.
function txtAllEalcoTwsFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtAllEalcoTwsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowseAllEalcoTwsFile.
function pbBrowseAllEalcoTwsFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbBrowseAllEalcoTwsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile( ...
    {  'EA_All*.*','Specify an EALCO TWSA LCC Grid Time Series file (*EALCO_TWSA*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify EALCO TWSA LCC Grid Time Series file');
 if  filename > 0  
    set(handles.txtAllEalcoTwsFile, 'String', pathname);
    %set(handles.txtAllEalcoTwsFile, 'String', [pathname filename]);
 end 


function txtAllScaledTwsFile_Callback(hObject, eventdata, handles)
% hObject    handle to txtAllScaledTwsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtAllScaledTwsFile as text
%        str2double(get(hObject,'String')) returns contents of txtAllScaledTwsFile as a double


% --- Executes during object creation, after setting all properties.
function txtAllScaledTwsFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtAllScaledTwsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowseAllScaledTwsFile.
function pbBrowseAllScaledTwsFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbBrowseAllScaledTwsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {  'DA_All*.*','Specify a Scaled TWSA LCC Grid Time Series file (*All_STWS*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify a Scaled TWSA LCC Grid Time Series file');
 if  filename > 0  
    % set(handles.txtAllScaledTwsFile, 'String', pathname);
    set(handles.txtAllScaledTwsFile, 'String', [pathname filename]);
 end 


function txtAllFinalTwsFile1_Callback(hObject, eventdata, handles)
% hObject    handle to txtAllFinalTwsFile1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtAllFinalTwsFile1 as text
%        str2double(get(hObject,'String')) returns contents of txtAllFinalTwsFile1 as a double


% --- Executes during object creation, after setting all properties.
function txtAllFinalTwsFile1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtAllFinalTwsFile1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowseAllFinalTwsFile1.
function pbBrowseAllFinalTwsFile1_Callback(hObject, eventdata, handles)
% hObject    handle to pbBrowseAllFinalTwsFile1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile( ...
    {  'DA_All*.*','Specify a Downscaled TWSA LCC Grid Time Series file 1 (*All_TWSA*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify a Downscaled TWSA LCC Grid Time Series file 1');
 if  filename > 0  
    set(handles.txtAllFinalTwsFile1, 'String', pathname);
    %set(handles.txtAllFinalTwsFile1, 'String', [pathname filename]);
 end 

function txtAllGswsFile1_Callback(hObject, eventdata, handles)
% hObject    handle to txtAllGswsFile1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtAllGswsFile1 as text
%        str2double(get(hObject,'String')) returns contents of txtAllGswsFile1 as a double


% --- Executes during object creation, after setting all properties.
function txtAllGswsFile1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtAllGswsFile1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowseAllGswsFile1.
function pbBrowseAllGswsFile1_Callback(hObject, eventdata, handles)
% hObject    handle to pbBrowseAllGswsFile1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {  'DA_All*.*','Specify a GSWSA LCC Grid Time Series file 1 (*All_GSWSA*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify a GSWSA LCC Grid Time Series file 1');
 if  filename > 0  
    set(handles.txtAllGswsFile1, 'String', pathname);
    %set(handles.txtAllGswsFile1, 'String', [pathname filename]);
 end 


function txtGWTSInputFolder_Callback(hObject, eventdata, handles)
% hObject    handle to txtGWTSInputFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtGWTSInputFolder as text
%        str2double(get(hObject,'String')) returns contents of txtGWTSInputFolder as a double


% --- Executes during object creation, after setting all properties.
function txtGWTSInputFolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtGWTSInputFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowseGWTSInputFolder.
function pbBrowseGWTSInputFolder_Callback(hObject, ~, handles)
% hObject    handle to pbBrowseGWTSInputFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile( ...
    {  '*.*','Specify Groundwater Well Level Observation data input folder (*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify Groundwater Well Level Observation data input folder');
 if  filename > 0  
    set(handles.txtGWTSInputFolder, 'String', pathname);
    %set(handles.txtGWTSInputFolder, 'String', [pathname filename]);
 end 


function txtGWLatLonFile_Callback(hObject, ~, handles)
% hObject    handle to txtGWLatLonFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtGWLatLonFile as text
%        str2double(get(hObject,'String')) returns contents of txtGWLatLonFile as a double


% --- Executes during object creation, after setting all properties.
function txtGWLatLonFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtGWLatLonFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowseGWLatLonFile.
function pbBrowseGWLatLonFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbBrowseGWLatLonFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {  '*.*','Specify Groundwater Well Lat Lon file (*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify Groundwater Well Lat Lon file');
 if  filename > 0  
    set(handles.txtGWLatLonFile, 'String', [pathname filename]);
 end 


% --- Executes on button press in pbValidateByBasinMGMWObs.
function pbValidateByBasinMGMWObs_Callback(hObject, eventdata, handles)
% hObject    handle to pbValidateByBasinMGMWObs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% record the start time
tic;

%% get GRACE TWS NetCDF input data file name
gtws_file = get(handles.txtGraceTwsFile, 'String');
if isempty(gtws_file)
    msgbox('Specify the GRACE TWS Anomaly NetCDF file.');
    return;
end
%% get JPL mascon id file name for Canada LCC grid
mask_file = get(handles.txtGraceTwsIdMaskFile, 'String');
if isempty(mask_file)
    msgbox('Specify the JPL mascon id file for Canada LCC grid.');
    return;
end
%% get Canada LCC grid Lat Lon position file
lcc_latlon_file = get(handles.txtLCCGridLatLonFile, 'String');
if isempty(lcc_latlon_file)
    msgbox('Specify the Canada LCC grid lat lon position file.');
    return;
end
%% get all grace tws file
all_gtws_file = get(handles.txtAllGraceTwsFile, 'String');
if isempty(all_gtws_file)
    msgbox('Specify the all GRACE TWS file.');
    return;
end
%% get all EALCO tws file
all_ttws_file = get(handles.txtAllEalcoTwsFile, 'String');
if isempty(all_ttws_file)
    msgbox('Specify the all EALCO TWS file.');
    return;
end
%% get all scaled GRACE tws file
all_stws_file = get(handles.txtAllScaledTwsFile, 'String');
if isempty(all_stws_file)
    msgbox('Specify the all scaled GRACE TWS file.');
    return;
end
%% get all assimilated GRACE tws file 1
all_atws_file1 = get(handles.txtAllFinalTwsFile1, 'String');
if isempty(all_atws_file1)
    msgbox('Specify the all assimilated GRACE TWS file 1.');
    return;
end
%% get all estimated GSWS file 1
all_gsws_file1 = get(handles.txtAllGswsFile1, 'String');
if isempty(all_gsws_file1)
    msgbox('Specify the all estimated GSWS file 1.');
    return;
end
%% get all assimilated GRACE tws file 2
all_atws_file2 = get(handles.txtAllFinalTwsFile2, 'String');
if isempty(all_atws_file2)
    msgbox('Specify the all assimilated GRACE TWS file 2.');
    return;
end
%% get all estimated GSWS file 2
all_gsws_file2 = get(handles.txtAllGswsFile2, 'String');
if isempty(all_gsws_file2)
    msgbox('Specify the all estimated GSWS file 2.');
    return;
end

%% get Basin average GWS observation file
gws_obs_file = get(handles.txtBasinGWSDataFile, 'String');
if isempty(gws_obs_file)
    msgbox('Specify the Basin average GWS observation file.');
    return;
end
%% get Basin average SWS observation file
sws_obs_file = get(handles.txtBasinSWSDataFile, 'String');
if isempty(sws_obs_file)
    msgbox('Specify the Basin average SWS observation file.');
    return;
end
%% get basin mask file
basin_mask_file = get(handles.txtBasinMaskFile, 'String');
if isempty(basin_mask_file)
    msgbox('Specify Basin Mask file.');
    return;
end
%% get time epoch file
time_epoch_file = get(handles.txtTimeEpochFile, 'String');
if isempty(time_epoch_file)
    msgbox('Specify the time epoch file.');
    return;
end
%% get the output data ENVI header file
out_file = get(handles.txtOutputFile, 'String');
if isempty(out_file)
    msgbox('Specify the output data ENVI header file.');
    return;
end
%% retrive the output data filder
[output_path, outfname, ext] = fileparts(out_file);

row = str2double(get(handles.txtRow, 'String'));    % row number of the GRACE TWS data
col = str2double(get(handles.txtCol, 'String'));    % col number of the GRACE TWS data
band = str2double(get(handles.txtBand, 'String'));  % the number of months 
bsn = str2double(get(handles.txtBasinNum, 'String'));  % the number of basins

%% read all GRACE TWS data
gtws = zeros(row,col,band);
fid = fopen(all_gtws_file);
for b=1:band
    gtws(:,:,b) = (fread(fid,[col, row],'float32')');
end
fclose(fid);
%% read all EALCO TWS data
ttws = zeros(row,col,band);
fid = fopen(all_ttws_file);
for b=1:band
    ttws(:,:,b) = (fread(fid,[col, row],'float32')');
end
fclose(fid);
%% read all simply scaled GRACE TWS data
stws = zeros(row,col,band);
fid = fopen(all_stws_file);
for b=1:band
    stws(:,:,b) = (fread(fid,[col, row],'float32')');
end
fclose(fid);
%% read all assimilated GRACE TWS data 1
atws1 = zeros(row,col,band);
fid = fopen(all_atws_file1);
for b=1:band
    atws1(:,:,b) = (fread(fid,[col, row],'float32')');
end
fclose(fid);
%% read all estimated GSWS data 1
gsws1 = zeros(row,col,band);
fid = fopen(all_gsws_file1);
for b=1:band
    gsws1(:,:,b) = (fread(fid,[col, row],'float32')');
end
fclose(fid);
%% read all assimilated GRACE TWS data 2
atws2 = zeros(row,col,band);
fid = fopen(all_atws_file2);
for b=1:band
    atws2(:,:,b) = (fread(fid,[col, row],'float32')');
end
fclose(fid);
%% read all estimated GSWS data 2
gsws2 = zeros(row,col,band);
fid = fopen(all_gsws_file2);
for b=1:band
    gsws2(:,:,b) = (fread(fid,[col, row],'float32')');
end
fclose(fid);

% %% read in the pixel id mask input data
% fid = fopen(mask_file);
%     mask = (fread(fid,[col, row],'float32')');
% fclose(fid);

%% read in the basin mask data
fid = fopen(basin_mask_file);
    bs_mask = (fread(fid,[col, row],'int8')');
fclose(fid);

%% read in the time epchs
t0 = readvars(time_epoch_file);

%% read in the observed GWS time series
obs_gsws = ones(band,bsn)*-32760;        % observed GSWS
grc_gsws0 = ones(band,bsn)*-32760;      % GRACE derived GWS+SWS without downscaling
grc_gsws1 = ones(band,bsn)*-32760;      % GRACE derived GWS+SWS with simple scale factor for downscaling
grc_gsws2 = ones(band,bsn)*-32760;      % GRACE derived GWS+SWS with iterative adjustment for downscaling

%% calculate GRACE derived GSWS
for i=1:bsn
    idx = find(bs_mask == i);
    if ~isempty(idx)
        for b = 1:band
            tmp_gsws0 = gtws(:,:,b) - ttws(:,:,b);
            tmp_gsws0 = tmp_gsws0(idx);
            tmp_gsws0(tmp_gsws0==-32760)=[];
            grc_gsws0(b,i) = mean(tmp_gsws0);
            
            tmp_gsws1 = stws(:,:,b) - ttws(:,:,b);
            tmp_gsws1 = tmp_gsws1(idx);
            tmp_gsws1(tmp_gsws1==-32760)=[];
            grc_gsws1(b,i) = mean(tmp_gsws1);
                       
            tmp_gsws2 = gsws1(:,:,b);
            tmp_gsws2 = tmp_gsws2(idx);
            tmp_gsws2(tmp_gsws2==-32760)=[];
            grc_gsws2(b,i) = mean(tmp_gsws2);
        end
    end
end

%% read in the observed GWS and SWS
gws = readmatrix(gws_obs_file);
sws = readmatrix(sws_obs_file);
gws = gws + sws;
gws(isnan(gws)) = -32760;

%% sync the time epoch
%% manually to match the time epcochs
obs_gsws(8:12,:) = gws(1:5,:);
obs_gsws(13:102,:) = gws(7:96,:);
obs_gsws(103:106,:) = gws(98:101,:);
obs_gsws(107:116,:) = gws(103:112,:);
obs_gsws(117:120,:) = gws(105:108,:);
obs_gsws(121:124,:) = gws(110:113,:);
obs_gsws(125:128,:) = gws(115:118,:);
obs_gsws(129:132,:) = gws(120:123,:);
obs_gsws(133:136,:) = gws(125:128,:);
obs_gsws(137:140,:) = gws(130:133,:);

st_logfile = [output_path '\Basin_Statistics.log'];
if exist(st_logfile,'file')
    delete st_logfile;
end

fstlog = fopen(st_logfile,'wt');

fprintf(fstlog,'%s\n', ' ');
fprintf(fstlog,'%s\n','Basin based model evaluation results of the SCVCM');
fprintf(fstlog,'%s\n','BasinID  NofObs  RMSE0       PR0     RMSE1       PR1     RMSE2       PR2   RMSE0_HP    PR0_HP  RMSE1_HP    PR1_HP  RMSE2_HP    PR2_HP');

%% perform evaluation by comparing RMSE and correlations
rmse = zeros(bsn,6);
r = zeros(bsn,6);
for i = 1:bsn
    x = obs_gsws(:,i);
    y0 = grc_gsws0(:,i);
    y1 = grc_gsws1(:,i);
    y2 = grc_gsws2(:,i);

    idx = find(x ~= -32760);
    idy = find(y2 ~= -32760);
    idxy = intersect(idx,idy);
    if ~isempty(idxy)
        t = t0(idxy);
        x = x(idxy)*10;   %x = x - mean(x);
        y0 = y0(idxy); %y0 = y0 - mean(y0);
        y1 = y1(idxy); %y1 = y1 - mean(y1);
        y2 = y2(idxy); %y2 = y2 - mean(y2);
        rmse(i,1) = sqrt(mean((y0-x).^2));
        rmse(i,2) = sqrt(mean((y1-x).^2));
        rmse(i,3) = sqrt(mean((y2-x).^2));
        r(i,1) = corr(x,y0);
        r(i,2) = corr(x,y1);
        r(i,3) = corr(x,y2);
    end 
    
    plot_data = 1;
    if plot_data >= 1
        figure;
        hold on;
        plot(t,x,'r-o',...
            'LineWidth',1,...
            'MarkerSize',3,...
            'MarkerEdgeColor','r',...
            'MarkerFaceColor',[0.5,0.5,0.5]);
        plot(t,y0,'c-o',...
            'LineWidth',1,...
            'MarkerSize',3,...
            'MarkerEdgeColor','c',...
            'MarkerFaceColor',[0.5,0.5,0.5]);
        plot(t,y1,'g-o',...
            'LineWidth',1,...
            'MarkerSize',3,...
            'MarkerEdgeColor','c',...
            'MarkerFaceColor',[0.5,0.5,0.5]);
        plot(t,y2,'b-*',...
            'LineWidth',1,...
            'MarkerSize',3,...
            'MarkerEdgeColor','b',...
            'MarkerFaceColor',[0.5,0.5,0.5]);
        %title([wid(i) ' rt = ' num2str(r(i,2))]);
        %title(wid(i));
        xlabel('Time in days from 2002-01-01')
        ylabel('GSWS Anomaly in EWH [mm]')
    end
    
    %% apply the HP filter
    x = hpfilter2(x,14400);
    y0 = hpfilter2(y0,14400);
    y1 = hpfilter2(y1,14400);
    y2 = hpfilter2(y2,14400);
    rmse(i,4) = sqrt(mean((y0-x).^2));
    rmse(i,5) = sqrt(mean((y1-x).^2));
    rmse(i,6) = sqrt(mean((y2-x).^2));

    r(i,4) = corr(x,y0);
    r(i,5) = corr(x,y1);
    r(i,6) = corr(x,y2);
    
    plot_data = 1;
    if plot_data >= 1
        figure;
        hold on;
        plot(t,x,'r-o',...
            'LineWidth',1,...
            'MarkerSize',3,...
            'MarkerEdgeColor','r',...
            'MarkerFaceColor',[0.5,0.5,0.5]);
        plot(t,y0,'c-o',...
            'LineWidth',1,...
            'MarkerSize',3,...
            'MarkerEdgeColor','c',...
            'MarkerFaceColor',[0.5,0.5,0.5]);
        plot(t,y1,'g-o',...
            'LineWidth',1,...
            'MarkerSize',3,...
            'MarkerEdgeColor','c',...
            'MarkerFaceColor',[0.5,0.5,0.5]);
        plot(t,y2,'b-*',...
            'LineWidth',1,...
            'MarkerSize',3,...
            'MarkerEdgeColor','b',...
            'MarkerFaceColor',[0.5,0.5,0.5]);
        %title([wid(i) ' rt = ' num2str(r(i,2))]);
        %title(wid(i));
        xlabel('Time in days from 2002-01-01')
        ylabel('GSWS Anomaly in EWH [mm]')
    end
    
    fprintf(fstlog,'%+5s', num2str(i,'% 5d'));
    fprintf(fstlog,'%+8s', num2str(length(idxy),'% 8d'));
    for j = 1:6
        fprintf(fstlog,'%+10s', num2str(rmse(i,j),'% 10.4f'));
        fprintf(fstlog,'%+10s', num2str(r(i,j),'% 10.4f'));
    end
    fprintf(fstlog,'%s\n', ' ');      
end


ttime = toc;

disp(['>>> Elapsed time: ', num2str(ttime),' seconds. <<<'])

msgbox('Completed TWS data analysis successfully!');

fclose(fstlog);


function txtAllFinalTwsFile2_Callback(hObject, eventdata, handles)
% hObject    handle to txtAllFinalTwsFile2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtAllFinalTwsFile2 as text
%        str2double(get(hObject,'String')) returns contents of txtAllFinalTwsFile2 as a double


% --- Executes during object creation, after setting all properties.
function txtAllFinalTwsFile2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtAllFinalTwsFile2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowseAllFinalTwsFile2.
function pbBrowseAllFinalTwsFile2_Callback(hObject, eventdata, handles)
% hObject    handle to pbBrowseAllFinalTwsFile2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile( ...
    {  'DA_All*.*','Specify a Downscaled TWSA LCC Grid Time Series file 2 (*All_TWSA*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify a Downscaled TWSA LCC Grid Time Series file 2');
 if  filename > 0  
    set(handles.txtAllFinalTwsFile2, 'String', pathname);
    %set(handles.txtAllFinalTwsFile2, 'String', [pathname filename]);
 end 


function txtAllGswsFile2_Callback(hObject, eventdata, handles)
% hObject    handle to txtAllGswsFile2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtAllGswsFile2 as text
%        str2double(get(hObject,'String')) returns contents of txtAllGswsFile2 as a double


% --- Executes during object creation, after setting all properties.
function txtAllGswsFile2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtAllGswsFile2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowseAllGswsFile2.
function pbBrowseAllGswsFile2_Callback(hObject, eventdata, handles)
% hObject    handle to pbBrowseAllGswsFile2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile( ...
    {  'DA_All*.*','Specify a GSWSA LCC Grid Time Series file 2 (*All_GSWSA*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify a GSWSA LCC Grid Time Series file 2');
 if  filename > 0  
    set(handles.txtAllGswsFile2, 'String', pathname);
    %set(handles.txtAllGswsFile2, 'String', [pathname filename]);
 end 

function txtBasinMaskFile_Callback(hObject, eventdata, handles)
% hObject    handle to txtBasinMaskFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtBasinMaskFile as text
%        str2double(get(hObject,'String')) returns contents of txtBasinMaskFile as a double


% --- Executes during object creation, after setting all properties.
function txtBasinMaskFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtBasinMaskFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowseBasinMaskFile.
function pbBrowseBasinMaskFile_Callback(~, eventdata, handles)
% hObject    handle to pbBrowseBasinMaskFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {  '*.*','Specify Water Basin Mask File (*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify Surface Water Basin Mask File');
 if  filename > 0  
    set(handles.txtBasinMaskFile, 'String', [pathname filename]);
 end 



function txtBasinNum_Callback(hObject, ~, handles)
% hObject    handle to txtBasinNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtBasinNum as text
%        str2double(get(hObject,'String')) returns contents of txtBasinNum as a double


% --- Executes during object creation, after setting all properties.
function txtBasinNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtBasinNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkAlberta.
function chkAlberta_Callback(hObject, eventdata, handles)
% hObject    handle to chkAlberta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkAlberta



function txtTimeEpochFile_Callback(hObject, ~, handles)
% hObject    handle to txtTimeEpochFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtTimeEpochFile as text
%        str2double(get(hObject,'String')) returns contents of txtTimeEpochFile as a double


% --- Executes during object creation, after setting all properties.
function txtTimeEpochFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtTimeEpochFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowseTimeEpochFile.
function pbBrowseTimeEpochFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbBrowseTimeEpochFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {  '*.*','Specify Time Epoch File (*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify Time Epoch File');
 if  filename > 0  
    set(handles.txtTimeEpochFile, 'String', [pathname filename]);
 end 


% --- Executes on button press in pbPlotImages.
function pbPlotImages_Callback(hObject, eventdata, handles)
% hObject    handle to pbPlotImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

row = str2double(get(handles.txtRow, 'String'));    % row number of the GRACE TWS data
col = str2double(get(handles.txtCol, 'String'));    % col number of the GRACE TWS data
band = str2double(get(handles.txtBand, 'String'));  % the number of months
total_month = band;

years = 2002:2016;
% idx = [1,7,18,30,42,54,66,78,90,102,112,122,131,140,149,158];

%% set the pixel range for the test area in Saschatchewan
% xbound = 623:770;
% ybound = 298:442;
xbound = 1:row;
ybound = 1:col;

%% remove the gray background color for paper figures
set(0,'DefaultTextInterpreter','none');
set(gcf, 'InvertHardCopy', 'off');

%% get all downsvaled GRACE tws file 1
all_atws_file1 = get(handles.txtAllFinalTwsFile1, 'String');
if isempty(all_atws_file1)
    msgbox('Specify the all assimilated GRACE TWS file (option 1).');
    return;
end
%% get all estimated GSWS file 1
all_gsws_file1 = get(handles.txtAllGswsFile1, 'String');
if isempty(all_gsws_file1)
    msgbox('Specify the all estimated GSWS file (option 1).');
    return;
end

%% get JPL mascon id file name for Canada LCC grid
mask_file = get(handles.txtGraceTwsIdMaskFile, 'String');
if isempty(mask_file)
    msgbox('Specify the JPL mascon id file for Canada LCC grid.');
    return;
end
%% read in the pixel id mask input data
fid = fopen(mask_file);
    mask = (fread(fid,[col, row],'float32')');
fclose(fid);

%% get the output data ENVI header file
output_path = get(handles.txtOutputFile, 'String');
if isempty(output_path)
    msgbox('Specify the output data folder.');
    return;
end


%% read all original GRACE TWSA data
% get all original grace twsa file names
gtwsa_folder = get(handles.txtAllGraceTwsFile, 'String');
if isempty(gtwsa_folder)
    msgbox('Specify the original GRACE TWSA data file folder.');
    return;
end
allfnames = dir(gtwsa_folder);
gtwsa_fnames = {allfnames(~[allfnames.isdir]).name};
n_gtwsa = length(gtwsa_fnames);
gtws1 = ones(row,col,n_gtwsa)*-32760;
gtws1_std = ones(row,col,n_gtwsa)*-32760;
time_list = zeros(n_gtwsa,1);
time_str_list = cell(n_gtwsa,1);
for i = 1:n_gtwsa
    gtwsai_file = [gtwsa_folder '\' gtwsa_fnames{i}];
    [time_list(i), time_str_list{i}] = daysFromFileName(gtwsai_file,3);
    fid = fopen(gtwsai_file);
        gtws1(:,:,i) = (fread(fid,[col, row],'float32')');
        gtws1_std(:,:,i) = (fread(fid,[col, row],'float32')');
    fclose(fid);
end  

%% read all EALCO TWSA data
% get all EALCO twsa file names
ttwsa_folder = get(handles.txtAllEalcoTwsFile, 'String');
if isempty(ttwsa_folder)
    msgbox('Specify the EALCO TWSA data file folder.');
    return;
end
allfnames = dir(ttwsa_folder);
ttwsa_fnames = {allfnames(~[allfnames.isdir]).name};
n_ttwsa = length(ttwsa_fnames);
ttws1 = ones(row,col,n_gtwsa)*-32760;
file_read_count = 0;
for i = 1:n_ttwsa
    ttwsai_file = [ttwsa_folder '\' ttwsa_fnames{i}];
    % sync the time stamps
    [time_day, ttws_time_str] = daysFromFileName(ttwsai_file,3);
    [tf, idx] = ismember(ttws_time_str,time_str_list);
    if tf
        file_read_count = file_read_count + 1;
        fid = fopen(ttwsai_file);
        ttws1(:,:,idx) = (fread(fid,[col, row],'float32')');
        fclose(fid);
    end
end
if file_read_count ~= n_gtwsa
    msgbox('Error: the GTWSA file number ~= the read TTWSA file number! Program stopped!');
    return;
end

%% read all spatially downscaled TWSA data 1
% get all saptially downscaled and assimilated TWSA file names
sdtwsa_folder = get(handles.txtAllFinalTwsFile1, 'String');
if isempty(sdtwsa_folder)
    msgbox('Specify the saptially downscaled TWSA file folder.');
    return;
end
allfnames = dir(sdtwsa_folder);
sdtwsa_fnames = {allfnames(~[allfnames.isdir]).name};
n_sdtwsa = length(sdtwsa_fnames);
atws1 = ones(row,col,n_gtwsa)*-32760;
atws1_std = ones(row,col,n_gtwsa)*-32760;
file_read_count = 0;
for i = 1:n_sdtwsa
    if ~contains(sdtwsa_fnames{i}, '.hdr')
        sdtwsai_file = [sdtwsa_folder '\' sdtwsa_fnames{i}];
        % sync the time stamps
        [time_day, atws_time_str] = daysFromFileName(sdtwsai_file,3);
        [tf, idx] = ismember(atws_time_str,time_str_list);
        if tf
            file_read_count = file_read_count + 1;
            fid = fopen(sdtwsai_file);
            atws1(:,:,idx) = (fread(fid,[col, row],'float32')');
            atws1_std(:,:,idx) = (fread(fid,[col, row],'float32')');
            fclose(fid);
        end
    end
end
if file_read_count ~= n_gtwsa
    msgbox('Error: the GTWSA file number ~= the read SDTWSA file number! Program stopped!');
    return;
end


%% read all estimated GSWSA data 1
gswsa_folder = get(handles.txtAllGswsFile1, 'String');
if isempty(sdtwsa_folder)
    msgbox('Specify the saptially downscaled GSWSA file folder.');
    return;
end
allfnames = dir(gswsa_folder);
gswsa_fnames = {allfnames(~[allfnames.isdir]).name};
n_gswsa = length(gswsa_fnames);
gsws1 = ones(row,col,n_gtwsa)*-32760;
gsws1_std = ones(row,col,n_gtwsa)*-32760;
file_read_count = 0;
for i = 1:n_gswsa
    if ~contains(gswsa_fnames{i}, '.hdr')
        gswsai_file = [gswsa_folder '\' gswsa_fnames{i}];
        % sync the time stamps
        [time_day, gswsa_time_str] = daysFromFileName(gswsai_file,3);
        [tf, idx] = ismember(gswsa_time_str,time_str_list);
        if tf
            file_read_count = file_read_count + 1;
            fid = fopen(gswsai_file);
            gsws1(:,:,idx) = (fread(fid,[col, row],'float32')');
            gsws1_std(:,:,idx) = (fread(fid,[col, row],'float32')');
            fclose(fid);
        end
    end
end
if file_read_count ~= n_gtwsa
    msgbox('Error: the GTWSA file number ~= the read GSWSA1 file number! Program stopped!');
    return;
end

%% set a unified data range for easy comparison
cmin = -300;
cmax = 300;
clims = [cmin,cmax];
cmin_std = 0;
cmax_std = 100;
clims_std = [cmin_std,cmax_std];

gtwsa1 = gtws1(xbound,ybound,1);
ttwsa1 = ttws1(xbound,ybound,1);
atwsa1 = atws1(xbound,ybound,1);
gswsa1 = gsws1(xbound,ybound,1);
nodata_idx = find(atwsa1 == -32760);
% gtwsa1(nodata_idx) = 0;
% ttwsa1(nodata_idx) = 0;
% atwsa1(nodata_idx) = 0;
% gswsa1(nodata_idx) = 0;
% cmin = min([min(gtwsa1(:)),min(ttwsa1(:)),min(atwsa1(:)),min(gswsa1(:))]);
% cmax = max([max(gtwsa1(:)),max(ttwsa1(:)),max(atwsa1(:)),max(gswsa1(:))]);
gtwsa1(nodata_idx)= cmin;
ttwsa1(nodata_idx)= cmin;
atwsa1(nodata_idx)= cmin;
gswsa1(nodata_idx)= cmin;

%% set the no data area as transparent
[r,c] = size(atwsa1);
A = ones(r,c); A(nodata_idx)= 0;

figure(1);
subplot(2,2,1)
imagesc(gtwsa1,clims);
alpha(A);
set(gca,'XTick',[], 'YTick', []);
colorbar
subplot(2,2,2)
imagesc(ttwsa1,clims);
alpha(A);
set(gca,'XTick',[], 'YTick', []);
colorbar
subplot(2,2,3)
imagesc(atwsa1,clims);
alpha(A);
set(gca,'XTick',[], 'YTick', []);
colorbar
subplot(2,2,4)
imagesc(gswsa1,clims);
alpha(A);
set(gca,'XTick',[], 'YTick', []);
colorbar
saveas(figure(1),[output_path 'Figure2a.png']);

gtwsa1 = gtws1(xbound,ybound,end);
ttwsa1 = ttws1(xbound,ybound,end);
atwsa1 = atws1(xbound,ybound,end);
gswsa1 = gsws1(xbound,ybound,end);
nodata_idx = find(atwsa1 == -32760);
% gtwsa1(nodata_idx) = 0;
% ttwsa1(nodata_idx) = 0;
% atwsa1(nodata_idx) = 0;
% gswsa1(nodata_idx) = 0;
% cmin = min([min(gtwsa1(:)),min(ttwsa1(:)),min(atwsa1(:)),min(gswsa1(:))]);
% cmax = max([max(gtwsa1(:)),max(ttwsa1(:)),max(atwsa1(:)),max(gswsa1(:))]);
gtwsa1(nodata_idx)= cmin;
ttwsa1(nodata_idx)= cmin;
atwsa1(nodata_idx)= cmin;
gswsa1(nodata_idx)= cmin;

%% set the no data area as transparent
[r,c] = size(gtwsa1);
A = ones(r,c); A(nodata_idx)= 0;

figure(2);
subplot(2,2,1)
imagesc(gtwsa1,clims);
alpha(A);
set(gca,'XTick',[], 'YTick', []);
colorbar
subplot(2,2,2)
imagesc(ttwsa1,clims);
alpha(A);
set(gca,'XTick',[], 'YTick', []);
colorbar
subplot(2,2,3)
imagesc(atwsa1,clims);
alpha(A);
set(gca,'XTick',[], 'YTick', []);
colorbar
subplot(2,2,4)
imagesc(gswsa1,clims);
alpha(A);
set(gca,'XTick',[], 'YTick', []);
colorbar
saveas(figure(2),[output_path 'Figure2b.png']);

%% plot the animated monthly and annual data
%% plot monthly image in an animated gif file
filename1 = [output_path 'ttwsa_monthly.gif'];
filename2 = [output_path 'sdtwsa_monthly.gif'];
filename3 = [output_path 'sdgwsa_monthly.gif'];
filename4 = [output_path 'grtwsa_monthly.gif'];
filename5 = [output_path 'grtwsa_std_monthly.gif'];
filename6 = [output_path 'sdtwsa_std_monthly.gif'];
filename7 = [output_path 'sdgswsa_std_monthly.gif'];
for i = 1:total_month
    ttwsa1 = ttws1(xbound,ybound,i);
    sdtwsa1 = atws1(xbound,ybound,i);
    sdgwsa1 = gsws1(xbound,ybound,i);
    grtwsa1 = gtws1(xbound,ybound,i);
    grtwsa_std1 = gtws1_std(xbound,ybound,i);
    sdtwsa_std1 = atws1_std(xbound,ybound,i);
    gswsa_std1 = gsws1_std(xbound,ybound,i);

    nodata_idx = find(sdtwsa1 == -32760);
    % ttwsa1(nodata_idx) = 0;
    % stdtwsa1(nodata_idx) = 0;
    % stdgwsa1(nodata_idx) = 0;
    % cmin = min([min(gtwsa1(:)),min(ttwsa1(:)),min(stdtwsa1(:)),min(stdgwsa1(:))]);
    % cmax = max([max(gtwsa1(:)),max(ttwsa1(:)),max(stdtwsa1(:)),max(stdgwsa1(:))]);
    ttwsa1(nodata_idx)= cmin;
    sdtwsa1(nodata_idx)= cmin;
    sdgwsa1(nodata_idx)= cmin;
    grtwsa1(nodata_idx)= cmin;
    grtwsa_std1(nodata_idx)= cmax_std/2;
    clims = [cmin,cmax];
    clims_std = [cmin_std,cmax_std];

    [r,c] = size(sdtwsa1);
    A = ones(r,c); A(nodata_idx)= 0;

    time_str = convertDaysToDateTime(time_list(i));
    time_str = datestr(time_str, 'yyyy-mm-dd');

    figure(3);
    set(gcf,'color','w');
    imagesc(ttwsa1,clims)
    alpha(A);
    set(gca,'XTick',[], 'YTick', []);
    set(gcf,'color','w');
    title(['EALCO TWSA: ' time_str]);
    colorbar
    drawnow
    frame = getframe(3);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
          imwrite(imind,cm,filename1,'gif', 'Loopcount',inf);
    else
          imwrite(imind,cm,filename1,'gif','WriteMode','append');
    end

    figure(4);
    set(gcf,'color','w');
    imagesc(sdtwsa1,clims)
    alpha(A);
    set(gca,'XTick',[], 'YTick', []);
    title(['SD TWSA: ' time_str]);
    colorbar
    drawnow
    frame = getframe(4);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
          imwrite(imind,cm,filename2,'gif', 'Loopcount',inf);
    else
          imwrite(imind,cm,filename2,'gif','WriteMode','append');
    end

    figure(5);
    set(gcf,'color','w');
    imagesc(sdgwsa1,clims)
    alpha(A);
    set(gca,'XTick',[], 'YTick', []);
    title(['EGSWSA: ' time_str]);
    colorbar
    drawnow
    frame = getframe(5);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
          imwrite(imind,cm,filename3,'gif', 'Loopcount',inf);
    else
          imwrite(imind,cm,filename3,'gif','WriteMode','append');
    end

    figure(6);
    set(gcf,'color','w');
    imagesc(grtwsa1,clims);
    alpha(A);
    set(gca,'XTick',[], 'YTick', []);
    title(['GRACE TWSA: ' time_str]);
    colorbar
    drawnow
    frame = getframe(6);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
          imwrite(imind,cm,filename4,'gif', 'Loopcount',inf);
    else
          imwrite(imind,cm,filename4,'gif','WriteMode','append');
    end

    figure(7);
    set(gcf,'color','w');
    imagesc(grtwsa_std1,clims_std)
    alpha(A);
    set(gca,'XTick',[], 'YTick', []);
    title(['GRACE TWSA Std: ' time_str]);
    colorbar
    drawnow
    frame = getframe(7);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
          imwrite(imind,cm,filename5,'gif', 'Loopcount',inf);
    else
          imwrite(imind,cm,filename5,'gif','WriteMode','append');
    end  
    
    figure(8);
    set(gcf,'color','w');
    imagesc(sdtwsa_std1,clims_std)
    alpha(A);
    set(gca,'XTick',[], 'YTick', []);
    title(['SD TWSA Std: ' time_str]);
    colorbar
    drawnow
    frame = getframe(8);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
          imwrite(imind,cm,filename6,'gif', 'Loopcount',inf);
    else
          imwrite(imind,cm,filename6,'gif','WriteMode','append');
    end  
    
    figure(9);
    set(gcf,'color','w');
    imagesc(gswsa_std1,clims_std)
    alpha(A);
    set(gca,'XTick',[], 'YTick', []);
    title(['EGSWSA Std: ' time_str]);
    colorbar
    drawnow
    frame = getframe(9);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
          imwrite(imind,cm,filename7,'gif', 'Loopcount',inf);
    else
          imwrite(imind,cm,filename7,'gif','WriteMode','append');
    end  
end

%% plot annual changes
filename1 = [output_path 'ttwsa_yearly.gif'];
filename2 = [output_path 'sdtwsa_yearly.gif'];
filename3 = [output_path 'sdgwsa_yearly.gif'];
filename4 = [output_path 'grtwsa_yearly.gif'];
idx = [1,7,18,30,42,54,66,78,90,102,112,122,131,140,149,158];
for i=1:length(years)
    gtws_annual = ones(length(xbound), length(ybound))*-32760;
    ttws_annual = ones(length(xbound), length(ybound))*-32760;
    atws_annual = ones(length(xbound), length(ybound))*-32760;
    gsws_annual = ones(length(xbound), length(ybound))*-32760;
    
    gtws_monthly = gtws1(xbound,ybound,idx(i):idx(i+1));
    ttws_monthly = ttws1(xbound,ybound,idx(i):idx(i+1));
    atws_monthly = atws1(xbound,ybound,idx(i):idx(i+1));
    gsws_monthly = gsws1(xbound,ybound,idx(i):idx(i+1));

    for ii=1:length(xbound)
        for jj=1:length(ybound)
            gtws = gtws_monthly(ii,jj,:);
            ttws = ttws_monthly(ii,jj,:);
            atws = atws_monthly(ii,jj,:);
            gsws = gsws_monthly(ii,jj,:);
            % exclude the invalid data
            gtws(gtws==-32760)=[];
            ttws(ttws==-32760)=[];
            atws(atws==-32760)=[];
            gsws(gsws==-32760)=[];
            if ~isempty(gtws)
                gtws_annual(ii,jj)=mean(gtws);
            else
                gtws_annual(ii,jj)=-32760;
            end
            if ~isempty(ttws)
                ttws_annual(ii,jj)=mean(ttws);
            else
                ttws_annual(ii,jj)=-32760;
            end
            if ~isempty(atws)
                atws_annual(ii,jj)=mean(atws);
            else
                atws_annual(ii,jj)=-32760;
            end
            if ~isempty(gsws)
                gsws_annual(ii,jj)=mean(gsws);
            else
                gsws_annual(ii,jj)=-32760;
            end
        end
    end
    nodata_idx = find(atws_annual == -32760);
%     gtws_annual(nodata_idx) = 0;
%     ttws_annual(nodata_idx) = 0;
%     atws_annual(nodata_idx) = 0;
%     gsws_annual(nodata_idx) = 0;
%     cmin = min([min(gtws_annual(:)),min(ttws_annual(:)),min(atws_annual(:)),min(gsws_annual(:))]);
%     cmax = max([max(gtws_annual(:)),max(ttws_annual(:)),max(atws_annual(:)),max(gsws_annual(:))]);
    gtws_annual(nodata_idx)= cmin;
    ttws_annual(nodata_idx)= cmin;
    atws_annual(nodata_idx)= cmin;
    gsws_annual(nodata_idx)= cmin;
    clims = [cmin,cmax];

    [r,c] = size(gtws_annual);
    A = ones(r,c); A(nodata_idx)= 0;      

    figure(10);
    set(gcf,'color','w');
    imagesc(ttws_annual,clims);
    alpha(A);
    set(gca,'XTick',[], 'YTick', []);
    title(['EALCO TWSA: ' num2str(years(i))]);
    colorbar
    drawnow
    frame = getframe(10);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
          imwrite(imind,cm,filename1,'gif', 'Loopcount',inf);
    else
          imwrite(imind,cm,filename1,'gif','WriteMode','append');
    end

    figure(11);
    set(gcf,'color','w');
    imagesc(atws_annual,clims)
    alpha(A);
    set(gca,'XTick',[], 'YTick', []);
    title(['SD TWSA: ' num2str(years(i))]);
    colorbar
    drawnow
    frame = getframe(11);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
          imwrite(imind,cm,filename2,'gif', 'Loopcount',inf);
    else
          imwrite(imind,cm,filename2,'gif','WriteMode','append');
    end

    figure(12);
    set(gcf,'color','w');
    imagesc(gsws_annual,clims)
    alpha(A);
    set(gca,'XTick',[], 'YTick', []);
    title(['EGSWSA: ' num2str(years(i))]);
    colorbar
    drawnow
    frame = getframe(12);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
          imwrite(imind,cm,filename3,'gif', 'Loopcount',inf);
    else
          imwrite(imind,cm,filename3,'gif','WriteMode','append');
    end

    figure(13)
    set(gcf,'color','w');
    imagesc(gtws_annual,clims);
    alpha(A);
    title(['GRACE TWSA: ' num2str(years(i))]);
    set(gca,'XTick',[], 'YTick', []);
    colorbar
    drawnow
    frame = getframe(13);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
          imwrite(imind,cm,filename4,'gif', 'Loopcount',inf);
    else
          imwrite(imind,cm,filename4,'gif','WriteMode','append');
    end
end


function txtSpecificYiel_Callback(hObject, eventdata, handles)
% hObject    handle to txtSpecificYiel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtSpecificYiel as text
%        str2double(get(hObject,'String')) returns contents of txtSpecificYiel as a double


% --- Executes during object creation, after setting all properties.
function txtSpecificYiel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtSpecificYiel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popFilterLength.
function popFilterLength_Callback(~, eventdata, handles)
% hObject    handle to popFilterLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popFilterLength contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popFilterLength


% --- Executes during object creation, after setting all properties.
function popFilterLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popFilterLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtBasinGWSDataFile_Callback(hObject, eventdata, handles)
% hObject    handle to txtBasinGWSDataFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtBasinGWSDataFile as text
%        str2double(get(hObject,'String')) returns contents of txtBasinGWSDataFile as a double


% --- Executes during object creation, after setting all properties.
function txtBasinGWSDataFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtBasinGWSDataFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowseBasinGWSDataFile.
function pbBrowseBasinGWSDataFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbBrowseBasinGWSDataFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {  '*.*','Specify Groundwater Storage Time Series (*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify Groundwater Storage Time Series  file');
 if  filename > 0  
    set(handles.txtBasinGWSDataFile, 'String', [pathname filename]);
 end 


function txtBasinSWSDataFile_Callback(hObject, eventdata, handles)
% hObject    handle to txtBasinSWSDataFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtBasinSWSDataFile as text
%        str2double(get(hObject,'String')) returns contents of txtBasinSWSDataFile as a double


% --- Executes during object creation, after setting all properties.
function txtBasinSWSDataFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtBasinSWSDataFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowseBasinSWSDataFile.
function pbBrowseBasinSWSDataFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbBrowseBasinSWSDataFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile( ...
    {  '*.*','Specify Surface Water Storage Time Series (*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify Surface Water Storage Time Series file');
 if  filename > 0  
    set(handles.txtBasinSWSDataFile, 'String', [pathname filename]);
 end 
