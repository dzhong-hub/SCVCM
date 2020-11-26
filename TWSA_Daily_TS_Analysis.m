function varargout = TWSA_Daily_TS_Analysis(varargin)
% TWSA_Daily_TS_Analysis MATLAB code for TWSA_Daily_TS_Analysis.fig
%      TWSA_Daily_TS_Analysis, by itself, creates a new TWSA_Daily_TS_Analysis or raises the existing
%      singleton*.
%
%      H = TWSA_Daily_TS_Analysis returns the handle to a new TWSA_Daily_TS_Analysis or the handle to
%      the existing singleton*.
%
%      TWSA_Daily_TS_Analysis('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TWSA_Daily_TS_Analysis.M with the given input arguments.
%
%      TWSA_Daily_TS_Analysis('Property','Value',...) creates a new TWSA_Daily_TS_Analysis or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TWSA_Daily_TS_Analysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TWSA_Daily_TS_Analysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TWSA_Daily_TS_Analysis

% Last Modified by GUIDE v2.5 25-Nov-2020 21:55:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TWSA_Daily_TS_Analysis_OpeningFcn, ...
                   'gui_OutputFcn',  @TWSA_Daily_TS_Analysis_OutputFcn, ...
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


% --- Executes just before TWSA_Daily_TS_Analysis is made visible.
function TWSA_Daily_TS_Analysis_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TWSA_Daily_TS_Analysis (see VARARGIN)

% Choose default command line output for TWSA_Daily_TS_Analysis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Add the path to PACE: Principal Analysis by Conditional Expectation
% addpath(genpath('C:\TWSModel\PACEV217\release2.17'));

% UIWAIT makes TWSA_Daily_TS_Analysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TWSA_Daily_TS_Analysis_OutputFcn(hObject, eventdata, handles) 
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
function pbClose_Callback(hObject, ~, ~)
% hObject    handle to pbClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close;

function txtGraceTwsStd_Callback(hObject, eventdata, handles)
% hObject    handle to txtGraceTwsStd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtGraceTwsStd as text
%        str2double(get(hObject,'String')) returns contents of txtGraceTwsStd as a double


% --- Executes during object creation, after setting all properties.
function txtGraceTwsStd_CreateFcn(hObject, ~, handles)
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
function popUnitTypeTerresTwsStd_CreateFcn(hObject, eventdata, ~)
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
function lstTerresTwsTimeSeries_Callback(hObject, eventdata, ~)
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
function popWeightMethod_Callback(hObject, ~, handles)
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




function txtLCCGridLatLonFile_Callback(~, eventdata, handles)
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
    {  'GR_All*.*','Specify a GRACE TWSA LCC Grid Monthly Time Series file (*GRACE_TWSA*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify a GRACE TWSA LCC Grid Monthly Time Series file');
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
    {  '*.*','Specify a EALCO TWSA LCC Grid Daily or Monthly Time Series data file (*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify a EALCO TWSA LCC Grid Daily or Monthly Time Series data file');
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
    {  'DA_All*.*','Specify Scaled TWS LCC Grid Time Series file (*All_STWS*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify Scaled TWS LCC Grid Time Series file');
 if  filename > 0  
    set(handles.txtAllScaledTwsFile, 'String', [pathname filename]);
 end 


function txtAllSpatialDownscaledTwsFile_Callback(hObject, eventdata, handles)
% hObject    handle to txtAllSpatialDownscaledTwsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtAllSpatialDownscaledTwsFile as text
%        str2double(get(hObject,'String')) returns contents of txtAllSpatialDownscaledTwsFile as a double


% --- Executes during object creation, after setting all properties.
function txtAllSpatialDownscaledTwsFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtAllSpatialDownscaledTwsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowseAllFinalTwsFile.
function pbBrowseAllFinalTwsFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbBrowseAllFinalTwsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile( ...
    {  'DA_All*.*','Specify a Spatially Downscaled TWSA LCC Grid Monthly Time Series file (*All_TWSA*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify a Spatially Downscaled TWSA LCC Grid Monthly Time Series file');
 if  filename > 0
    set(handles.txtAllSpatialDownscaledTwsFile, 'String', pathname);
    % set(handles.txtAllSpatialDownscaledTwsFile, 'String', [pathname filename]);
 end 

function txtAllSpatialDownscaledGswsFile_Callback(hObject, eventdata, handles)
% hObject    handle to txtAllSpatialDownscaledGswsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtAllSpatialDownscaledGswsFile as text
%        str2double(get(hObject,'String')) returns contents of txtAllSpatialDownscaledGswsFile as a double


% --- Executes during object creation, after setting all properties.
function txtAllSpatialDownscaledGswsFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtAllSpatialDownscaledGswsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowseAllGswsFile.
function pbBrowseAllGswsFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbBrowseAllGswsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {  'DA_All*.*','Specify a GSWSA LCC Grid Monthly Time Series file (*All_GSWSA*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify a GSWSA LCC Grid Monthly Time Series file');
 if  filename > 0  
    set(handles.txtAllSpatialDownscaledGswsFile, 'String', pathname);
    %set(handles.txtAllSpatialDownscaledGswsFile, 'String', [pathname filename]);
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
    {  '*.*','Specify a Groundwater Well Level Observation data file (*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify a Groundwater Well Level Observation data file');
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

function txtAllTemporalDownscaledTwsFile_Callback(hObject, eventdata, handles)
% hObject    handle to txtAllTemporalDownscaledTwsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtAllTemporalDownscaledTwsFile as text
%        str2double(get(hObject,'String')) returns contents of txtAllTemporalDownscaledTwsFile as a double


% --- Executes during object creation, after setting all properties.
function txtAllTemporalDownscaledTwsFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtAllTemporalDownscaledTwsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowseGWSObsFile.
function pbBrowseGWSObsFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbBrowseGWSObsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {  '*.*','Specify Groundwater Storage Time Series (*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify Groundwater Storage Time Series  file');
 if  filename > 0  
    set(handles.txtAllTemporalDownscaledTwsFile, 'String', pathname);
    %set(handles.txtAllTemporalDownscaledTwsFile, 'String', [pathname filename]);
 end 


function txtAllTemporalDownscaledGwsFile_Callback(hObject, eventdata, handles)
% hObject    handle to txtAllTemporalDownscaledGwsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtAllTemporalDownscaledGwsFile as text
%        str2double(get(hObject,'String')) returns contents of txtAllTemporalDownscaledGwsFile as a double


% --- Executes during object creation, after setting all properties.
function txtAllTemporalDownscaledGwsFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtAllTemporalDownscaledGwsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowseSWSObsFile.
function pbBrowseSWSObsFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbBrowseSWSObsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile( ...
    {  '*.*','Specify a STD GSWSA Daily Time Series data file (*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify a STD GSWSA Daily Time Series data file');
 if  filename > 0  
    set(handles.txtAllTemporalDownscaledGwsFile, 'String', pathname);
    %set(handles.txtAllTemporalDownscaledGwsFile, 'String', [pathname filename]);
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
function pbBrowseBasinMaskFile_Callback(hObject, eventdata, handles)
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
function popFilterLength_Callback(hObject, eventdata, handles)
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


% --- Executes on button press in pbValidateByGridDGMWObs.
function pbValidateByGridDGMWObs_Callback(hObject, eventdata, handles)
% hObject    handle to pbValidateByGridDGMWObs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% record the start time
tic;

%% get meta data
row = str2double(get(handles.txtRow, 'String'));    % row number of the GRACE TWS data
col = str2double(get(handles.txtCol, 'String'));    % col number of the GRACE TWS data
band = str2double(get(handles.txtBand, 'String'));  % the number of months
total_month = band;

sy = str2double(get(handles.txtSpecificYiel, 'String'));        % the specific yield value
daily_files = get(handles.popFileType, 'Value');
if daily_files == 1
    fl = 360*get(handles.popFilterLength, 'Value')+1;           % The filter length
else
    fl = 12*get(handles.popFilterLength, 'Value')+1;            % The filter length
end

%% specify mascon ids to be processed
% maskvalues = [273,274,308]'; % for Saschatchewan
% maskvalues = [162,163,203,204,205,206,273,274,307,308,341]'; % for Ontario
maskvalues = [234,235,236,271,272,305,306,307,339,340]'; % for alberta

%% get GRACE TWS NetCDF input data file name
gtws_file = get(handles.txtGraceTwsFile, 'String');
if isempty(gtws_file)
    msgbox('Specify the GRACE TWS Anomaly NetCDF file.');
    return;
end
%% read the GRACE TWS data from NETCDF file to get the time stamps
%ncdisp(gtws_file)
%day = ncread(gtws_file,'time');
day_bounds = ncread(gtws_file,'time_bounds');

%% get all EALCO daily or monthly twsa file names
ttws_folder = get(handles.txtAllEalcoTwsFile, 'String');
if isempty(ttws_folder)
    msgbox('Specify the EALCO Daily TWSA data file folder.');
    return;
end

%% get groundwater well level time series file names 
gmwl_obs_type = get(handles.popGMWObsType, 'Value'); % 1->water level height, 2->water level depth
% get GMWL daily observation input file folder
gmwl_in_dir = get(handles.txtGWTSInputFolder, 'String');
if isempty(gmwl_in_dir)
    msgbox('Specify the GMWL observations input folder.');
    return;
end
allfnames = dir(gmwl_in_dir);
gmwl_fnames = {allfnames(~[allfnames.isdir]).name};
n_gmwl = length(gmwl_fnames);

%% get the output file folder
output_path = get(handles.txtOutputFile, 'String');
if isempty(output_path)
    msgbox('Specify the output data folder.');
    return;
end

%% read in the JPL mascon id mask data
% get JPL mascon id file name for Canada LCC grid
mask_file = get(handles.txtGraceTwsIdMaskFile, 'String');
if isempty(mask_file)
    msgbox('Specify the JPL mascon id file for Canada LCC grid.');
    return;
end
fid = fopen(mask_file);
    mask = (fread(fid,[col, row],'float32')');
fclose(fid);

%% load lcc grid lat lon positions
% get Canada LCC grid Lat Lon position file
lcc_latlon_file = get(handles.txtLCCGridLatLonFile, 'String');
if isempty(lcc_latlon_file)
    msgbox('Specify the Canada LCC grid lat lon position file.');
    return;
end
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
gr_twsa = zeros(row,col,n_gtwsa);
% time_list = zeros(n_gtwsa,1);
% time_str_list = cell(n_gtwsa,1);
for i = 1:n_gtwsa
    gtwsai_file = [gtwsa_folder '\' gtwsa_fnames{i}];
    % [time_list(i), time_str_list{i}] = daysFromFileName(gtwsai_file,1);
    fid = fopen(gtwsai_file);
        gr_twsa(:,:,i) = (fread(fid,[col, row],'float32')');
    fclose(fid);
end  

%% read all spatially downscaled TWSA data
% get all saptially downscaled and assimilated TWSA file names
sdtwsa_folder = get(handles.txtAllSpatialDownscaledTwsFile, 'String');
if isempty(sdtwsa_folder)
    msgbox('Specify the saptially downscaled TWSA file folder.');
    return;
end
allfnames = dir(sdtwsa_folder);
sdtwsa_fnames = {allfnames(~[allfnames.isdir]).name};
n_sdtwsa = length(sdtwsa_fnames);
sd_twsa = zeros(row,col,n_sdtwsa);
% time_list = zeros(n_sdtwsa,1);
% time_str_list = cell(n_sdtwsa,1);
for i = 1:n_sdtwsa
    sdtwsai_file = [sdtwsa_folder '\' sdtwsa_fnames{i}];
    % [time_list(i), time_str_list{i}] = daysFromFileName(sdtwsai_file,1);
    fid = fopen(sdtwsai_file);
        sd_twsa(:,:,i) = (fread(fid,[col, row],'float32')');
    fclose(fid);
end  

% %% read all spatiotemporally downscaled TWSA data (no use)
% stdtwsa_folder = get(handles.txtAllTemporalDownscaledTwsFile, 'String');
% if isempty(stdtwsa_folder)
%     msgbox('Specify the saptiotemporally downscaled Daily TWSA data file folder.');
%     return;
% end
% allfnames = dir(stdtwsa_folder);
% stdtwsa_fnames = {allfnames(~[allfnames.isdir]).name};
% n_stdtwsa = length(stdtwsa_fnames);
% std_twsa = zeros(row,col,n_stdtwsa);
% time_list = zeros(n_stdtwsa,1);
% time_str_list = cell(n_stdtwsa,1);
% for i = 1:n_stdtwsa
%     std_twsi_file = [stdtwsa_folder '\' stdtwsa_fnames{i}];
%     [time_list(i), time_str_list{i}] = daysFromFileName(std_twsi_file,1);
%     fid = fopen(std_twsi_file);
%     std_twsa(:,:,i) = (fread(fid,[col, row],'float32')');
%     fclose(fid);
% end  

%% read all spatiotemporally downscaled and estimated GSWSA data
gwsa_folder = get(handles.txtAllTemporalDownscaledGwsFile, 'String');
if isempty(gwsa_folder)
    msgbox('Specify the Estimated Daily GWSA data file folder.');
    return;
end
allfnames = dir(gwsa_folder);
gwsa_fnames = {allfnames(~[allfnames.isdir]).name};
n_gwsa = length(gwsa_fnames);
gsws = zeros(row,col,n_gwsa);
time_list = zeros(n_gwsa,1);
time_str_list = cell(n_gwsa,1);
for i = 1:n_gwsa
    gwsi_file = [gwsa_folder '\' gwsa_fnames{i}];
    [time_list(i), time_str_list{i}] = daysFromFileName(gwsi_file,1);
    fid = fopen(gwsi_file);
        gsws(:,:,i) = (fread(fid,[col, row],'float32')');
    fclose(fid);
end  

%% read all daily matched EALCO TWSA data
allfnames = dir(ttws_folder);
ttws_fnames = {allfnames(~[allfnames.isdir]).name};
n_ttws = length(ttws_fnames);
ttws = ones(row,col,n_gwsa)*-32760;
file_read_count = 0;
for i = 1:n_ttws
    ttwsi_file = [ttws_folder '\' ttws_fnames{i}];
    % sync the time stamps
    [time_day, ttws_time_str] = daysFromFileName(ttwsi_file,1);
    [tf, idx] = ismember(ttws_time_str,time_str_list);
    if tf
        file_read_count = file_read_count + 1;
        fid = fopen(ttwsi_file);
        ttws(:,:,idx) = (fread(fid,[col, row],'float32')');
        fclose(fid);
    end
end
if file_read_count ~= n_gwsa
    msgbox('Error: the GSWSA file number ~= the read TTWSA file number! Program stopped!');
    return;
end
%% reset the n_ttws to the n_gwsa
n_ttws = n_gwsa;
    
%% open a log file for statistic information
st_logfile = [output_path '\Grid_Statistics.log'];
if exist(st_logfile,'file')
    delete st_logfile;
end
fstlog = fopen(st_logfile,'wt');

%% read in lat lon positions of groundwater wells
% get GMWL lat lon position file
gmwl_pos_file = get(handles.txtGWLatLonFile, 'String');
if isempty(gmwl_pos_file)
    msgbox('Specify the GMWL Lat Lon position file.');
    return;
end
%[wid,lat,lon,y1,y2,ym] = readvars(gmwl_pos_file);
[wid,lat,lon] = readvars(gmwl_pos_file);
nwid = length(wid);    

%% convert the well id from double to string
if isa(wid,'double')
    wid = num2str(wid);
end

fprintf(fstlog,'%s\n', ' ');
fprintf(fstlog,'%s\n','Temporal correlation comparison results for the estimated and measured GSWS');
fprintf(fstlog,'%s\n','           WellID           Lat         Lon          Mascon     Pixeli      Pixelj       N(days)   b1(yield)    RMSE1       RMSE2       RMSE3         R1          R2          R3        RMSE1       RMSE2       RMSE3         R1          R2          R3');

%% process data well by well
QI = zeros(nwid,13);
gws0 = ones(nwid, n_ttws)*-32760;
gws1 = ones(nwid, n_ttws)*-32760;
gws2 = ones(nwid, n_ttws)*-32760;
gws3 = ones(nwid, n_ttws)*-32760;

gws0_hp = ones(nwid, n_ttws)*-32760;
gws1_hp = ones(nwid, n_ttws)*-32760;
gws2_hp = ones(nwid, n_ttws)*-32760;
gws3_hp = ones(nwid, n_ttws)*-32760;

for i =1:nwid
    find_wid = 0;
    gtws1 = ones(n_ttws,1)*-32760;
    gtws2 = ones(n_ttws,1)*-32760;

    for j = 1:n_gmwl
        if contains(gmwl_fnames{j},wid(i,:))
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
                gsws_i = gsws(pixi,pixj,:);     % spatiotemporally dowscaled and estimated gsws
                twsa_i = ttws(pixi,pixj,:);     % EALCO TWSA daily data at mascon_id
                gtws_i = gr_twsa(pixi,pixj,:);  % original grace twsa
                atws_i = sd_twsa(pixi,pixj,:);  % spatially dowscaled and estimated twsa                
                
                gsws_i = gsws_i(:); 
                twsa_i = twsa_i(:); 
                gtws_i = gtws_i(:); 
                atws_i = atws_i(:); 
               
                %% calculate daily GSWS from the original GRACE TWSA and the spatially downscaled TWSA
                for j=1:total_month         
                    %% determine the day ranges for the monthly dataset j
                    day_bounds_j = day_bounds(:,j);
                    day_start = day_bounds_j(1);
                    day_end = day_bounds_j(2);

                    idx_days = find(time_list>=day_start & time_list<=day_end); % the start and end days should be included.
                    if ~isempty(idx_days)
                        temp_ones = ones(length(idx_days),1);
                        gtws1(idx_days)= temp_ones*gtws_i(j);
                        gtws2(idx_days)= temp_ones*atws_i(j);
                    end
                end                               
                
                %% read in GMWL observations
                csvfn = [gmwl_in_dir '\' gmwl_fnames{find_wid}];
                [tw,gmwl_i] = readvars(csvfn);
                
                tg = time_list; 
                
                %% remove NaN data in the time series
                idxg1 = find(gtws1 == -32760);
                idxg2 = find(gtws2 == -32760);
                idxg3 = find(gsws_i == -32760); %% Note: idxg1, idxg2 and idxg3 should be same.
                if ~isempty(idxg3)
                    gtws1(idxg3) = [];
                    gtws2(idxg3) = [];
                    gsws_i(idxg3) = [];
                    twsa_i(idxg3) = [];
                    tg(idxg3) = [];
                end
                
                idxw = find(gmwl_i == -32760);
                if ~isempty(idxw)
                    gmwl_i(idxw) = [];
                    tw(idxw) = [];
                end
                
                %% compare the three solutions by the temporal correlations
                [t,idxw,idxg] = intersect(tw,tg); 
                if ~isempty(t)
                    if gmwl_obs_type == 1
                        gsws0 = gmwl_i(idxw)*1000;   % for water height observations
                    else
                        gsws0 = -gmwl_i(idxw)*1000;   % for water depth observations
                    end
                    gsws1 = gtws1(idxg) - twsa_i(idxg); % gsws 1
                    gsws2 = gtws2(idxg) - twsa_i(idxg); % gsws 2
                    gsws3 = gsws_i(idxg);               % gsws 3
                    
                    n = length(t);

                    %% deduct the mean values
                    g0 = gsws0 - mean(gsws0);
                    g1 = gsws1 - mean(gsws1);
                    g2 = gsws2 - mean(gsws2);
                    g3 = gsws3 - mean(gsws3);
                    
                    %% apply HP filter
                    %% issue: how to avoid the discontinuty problem
                    % g0_hp = hpfilter2(g0,14400);
                    % g1_hp = hpfilter2(g1,14400);
                    % g2_hp = hpfilter2(g2,14400);
                    % g3_hp = hpfilter2(g3,14400);

                    B = 1/fl*ones(fl,1);
                    g0_hp = filter(B,1,g0);
                    g1_hp = filter(B,1,g1);
                    g2_hp = filter(B,1,g2);
                    g3_hp = filter(B,1,g3);

                    %% convert gw water level height to GWS
                    switch sy
                        case 1
                           sv = g0_hp\g1_hp; 
                        case 2
                           sv = g0_hp\g2_hp;
                        case 3
                           sv = g0_hp\g3_hp;
                        otherwise
                           sv = sy;
                    end 
                    g0 = sv*g0;
                    g0_hp = sv*g0_hp;

                    %% plot out the unfiltered and filtered data fora comparison
                    widstr = wid(i,:); 
%                     if  contains(widstr,'931')||contains(widstr,'229')||contains(widstr,'105') || contains(widstr,'260') || contains(widstr,'252')||contains(widstr,'251')||contains(widstr,'369') ||contains(widstr,'960')
%                     if contains(wid(i),'369') || contains(wid(i),'263') || contains(wid(i),'252')
%                         figure;
%                         hold on;
%                         plot(t,g0,'r-o',...
%                             'LineWidth',1,...
%                             'MarkerSize',3,...
%                             'MarkerEdgeColor','r',...
%                             'MarkerFaceColor',[0.5,0.5,0.5]);
%                         plot(t,g1,'c-o',...
%                             'LineWidth',1,...
%                             'MarkerSize',3,...
%                             'MarkerEdgeColor','c',...
%                             'MarkerFaceColor',[0.5,0.5,0.5]);
%                         plot(t,g2,'b-*',...
%                             'LineWidth',1,...
%                             'MarkerSize',3,...
%                             'MarkerEdgeColor','b',...
%                             'MarkerFaceColor',[0.5,0.5,0.5]);
%                         plot(t,g3,'g-*',...
%                             'LineWidth',1,...
%                             'MarkerSize',3,...
%                             'MarkerEdgeColor','g',...
%                             'MarkerFaceColor',[0.5,0.5,0.5]);                        
%                         title(['Well #: ' widstr]);
%                         xlabel('Time in days from 2002-01-01')
%                         ylabel('GWSA (mm EWH)') 
%                         legend('OGWSA','EGWSA1','EGWSA2','EGWSA3')
% 
%                         figure;
%                         hold on;
%                         plot(t,g0_hp,'r-o',...
%                             'LineWidth',1,...
%                             'MarkerSize',3,...
%                             'MarkerEdgeColor','r',...
%                             'MarkerFaceColor',[0.5,0.5,0.5]);
%                         plot(t,g1_hp,'c-o',...
%                             'LineWidth',1,...
%                             'MarkerSize',3,...
%                             'MarkerEdgeColor','c',...
%                             'MarkerFaceColor',[0.5,0.5,0.5]);
%                         plot(t,g2_hp,'b-*',...
%                             'LineWidth',1,...
%                             'MarkerSize',3,...
%                             'MarkerEdgeColor','b',...
%                             'MarkerFaceColor',[0.5,0.5,0.5]);
%                         plot(t,g3_hp,'g-*',...
%                             'LineWidth',1,...
%                             'MarkerSize',3,...
%                             'MarkerEdgeColor','g',...
%                             'MarkerFaceColor',[0.5,0.5,0.5]);
%                         title(['Well #: ' widstr]);
%                         xlabel('Time in days from 2002-01-01')
%                         ylabel('Smoothed GWSA (mm EWH)') 
%                         legend('OGWSA','EGWSA1','EGWSA2','EGWSA3')
%                     end 
                    
                    %% save GMW data for daily trend changes plots
                    used_well_file = [output_path '\Well_' widstr '_in_mascon_' num2str(mascon_id) '_daily_trend_data.txt'];
                    if exist(used_well_file,'file')
                        delete used_well_file;
                    end
                    fwell = fopen(used_well_file,'wt');
                    fprintf(fwell,'%s\n','             Id       Date        TimeDays      EGWSA1      EGWSA2      EGWSA3       OGWSA      EGWSA1      EGWSA2      EGWSA3      OGWSA');
                    for k = 1:length(g0)
                        fprintf(fwell,'%+16s', widstr);
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
                                            

                    %% calculate RMSE and temporal correlation coefficients
                    QI(i,1) = sv;                        
                    QI(i,2) = sqrt(mean((g1-g0).^2));
                    QI(i,3) = sqrt(mean((g2-g0).^2));
                    QI(i,4) = sqrt(mean((g3-g0).^2));
                    QI(i,5) = corr(g1,g0); 
                    QI(i,6) = corr(g2,g0);
                    QI(i,7) = corr(g3,g0);
                    QI(i,8) = sqrt(mean((g1_hp-g0_hp).^2));
                    QI(i,9) = sqrt(mean((g2_hp-g0_hp).^2));
                    QI(i,10) = sqrt(mean((g3_hp-g0_hp).^2));
                    QI(i,11) = corr(g1_hp,g0_hp); 
                    QI(i,12) = corr(g2_hp,g0_hp);
                    QI(i,13) = corr(g3_hp,g0_hp);
                    
                    fprintf(fstlog,'%+6s', num2str(i,'% 6d'));
                    fprintf(fstlog,'%+16s', wid(i,:));

                    fprintf(fstlog,'%+12s', num2str(lat(i),'% 12.6f'));
                    fprintf(fstlog,'%+12s', num2str(lon(i),'% 12.6f'));
                    fprintf(fstlog,'%+12s', num2str(mascon_id,'% 8d'));
                    fprintf(fstlog,'%+12s', num2str(pixi,'% 8d'));
                    fprintf(fstlog,'%+12s', num2str(pixj,'% 8d'));
                    fprintf(fstlog,'%+12s', num2str(n,'% 12d'));
                    for j=1:13
                        fprintf(fstlog,'%+12s', num2str(QI(i,j),'% 12.4f'));
                    end
                    fprintf(fstlog,'%s\n', ' ');
                    
                    %% keep a copy of the time series for spatial correlation comparison
                    % determine the idx of t in time_list
                    [t,idxw,idx] = intersect(t,time_list); % suppose that t should be all included.
                    % copy the data to the corresponding idx in time list
                    gws0(i,idx) = g0';
                    gws1(i,idx) = g1';
                    gws2(i,idx) = g2';
                    gws3(i,idx) = g3';

                    gws0_hp(i,idx) = g0_hp';
                    gws1_hp(i,idx) = g1_hp';
                    gws2_hp(i,idx) = g2_hp';
                    gws3_hp(i,idx) = g3_hp';
                end
            end
        end
    end
end

%% print out mean values
fprintf(fstlog,'%s\n', ' ');
fprintf(fstlog,'%+94s', '   Mean Value  ');
for j=1:13
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
title('Temporal RMSE Comparison');
xlabel('Well ID')
ylabel('RMSD in EWH [mm]') 
legend('EGWSA1','EGWSA2','EGWSA3')

%% plot temporal Peason's correlation comparison
figure;
hold on;
plot(QI(:,5),'r-o',...
    'LineWidth',1,...
    'MarkerSize',3,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
plot(QI(:,6),'c-o',...
    'LineWidth',1,...
    'MarkerSize',3,...
    'MarkerEdgeColor','c',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
plot(QI(:,7),'g-o',...
    'LineWidth',1,...
    'MarkerSize',3,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
title('Temporal Pearson Correlation Comparison');
xlabel('Well ID')
ylabel('Pearson Correlation') 
legend('EGWSA1','EGWSA2','EGWSA3')

%% scatter plot
x = gws0(:); 
y1 = gws1(:); 
y2 = gws2(:); 
y3 = gws3(:);  

% x_hp = gws0_hp(:);
% y1_hp = gws1_hp(:);
% y2_hp = gws2_hp(:);
% y3_hp = gws3_hp(:); 
    
idxg = find(x ~= -32760);
idxw = find(y1 ~= -32760);
idx = intersect(idxg, idxw);                
if ~isempty(idx)
    x = x(idx); 
    y1 = y1(idx);
    y2 = y2(idx); 
    y3 = y3(idx); 

%     x_hp = x_hp(idx); 
%     y1_hp = y1_hp(idx);
%     y2_hp = y2_hp(idx); 
%     y3_hp = y3_hp(idx);

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
fprintf(fstlog,'%s\n','    Id       Date        TimeDays     N(wells)    RMSE1       RMSE2       RMSE3         R1          R2          R3        RMSE1       RMSE2       RMSE3         R1          R2          R3');

QIs = zeros(n_ttws,13);
for i = 1:n_ttws
    g0 = gws0(:,i); 
    g1 = gws1(:,i); 
    g2 = gws2(:,i); 
    g3 = gws3(:,i);  
        
    g0_hp = gws0_hp(:,i);
    g1_hp = gws1_hp(:,i);
    g2_hp = gws2_hp(:,i);
    g3_hp = gws3_hp(:,i); 
    
    idxg = find(g0 ~= -32760);
    idxw = find(g1 ~= -32760);
    idx = intersect(idxg, idxw);                
    if ~isempty(idx)
        g0 = g0(idx); 
        g1 = g1(idx);
        g2 = g2(idx); 
        g3 = g3(idx); 
        
        g0_hp = g0_hp(idx); 
        g1_hp = g1_hp(idx);
        g2_hp = g2_hp(idx); 
        g3_hp = g3_hp(idx);
        
        QIs(i,1) = 1;                        
        QIs(i,2) = sqrt(mean((g1-g0).^2));
        QIs(i,3) = sqrt(mean((g2-g0).^2));
        QIs(i,4) = sqrt(mean((g3-g0).^2));
        QIs(i,5) = corr(g1,g0); 
        QIs(i,6) = corr(g2,g0);
        QIs(i,7) = corr(g3,g0);
        QIs(i,8) = sqrt(mean((g1_hp-g0_hp).^2));
        QIs(i,9) = sqrt(mean((g2_hp-g0_hp).^2));
        QIs(i,10) = sqrt(mean((g3_hp-g0_hp).^2));
        QIs(i,11) = corr(g1_hp,g0_hp); 
        QIs(i,12) = corr(g2_hp,g0_hp);
        QIs(i,13) = corr(g3_hp,g0_hp);
        
        n = length(idx);
    else
        n = 0;
    end
    
    fprintf(fstlog,'%+6s', num2str(i,'% 6d'));
    fprintf(fstlog,'%+14s', convertDaysToDateTime(time_list(i)));
    fprintf(fstlog,'%+12s', num2str(time_list(i),'% 12.1f'));
    fprintf(fstlog,'%+12s', num2str(n,'% 12d'));
    for j=2:13
        fprintf(fstlog,'%+12s', num2str(QIs(i,j),'% 12.4f'));
    end
    fprintf(fstlog,'%s\n', ' ');
end

%% print out mean values
fprintf(fstlog,'%s\n', ' ');
fprintf(fstlog,'%+44s', '   Mean Value  ');
for j=2:13
    fprintf(fstlog,'%+12s', num2str(mean(QIs(:,j)),'% 12.4f'));
end
fprintf(fstlog,'%s\n', ' ');

%% Plot Spatial RMSE Comparison
figure;
hold on;
plot(QIs(:,2),'r*',...
    'LineWidth',1,...
    'MarkerSize',3,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
plot(QIs(:,3),'c*',...
    'LineWidth',1,...
    'MarkerSize',3,...
    'MarkerEdgeColor','c',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
plot(QIs(:,4),'g*',...
    'LineWidth',1,...
    'MarkerSize',3,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
title('Spatial RMSE Comparison');
xlabel('Time Epoch in Month')
ylabel('RMSD in EWH [mm]') 
legend('EGWSA1','EGWSA2','EGWSA3')

%% Plot Spatial Pearson Correlation Comparison
figure;
hold on;
plot(QIs(:,5),'r*',...
    'LineWidth',1,...
    'MarkerSize',3,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
plot(QIs(:,6),'c*',...
    'LineWidth',1,...
    'MarkerSize',3,...
    'MarkerEdgeColor','c',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
plot(QIs(:,7),'g*',...
    'LineWidth',1,...
    'MarkerSize',3,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
title('Spatial Pearson Correlation Comparison');
xlabel('Time Epoch in Month')
ylabel('Pearson Correlation') 
legend('EGWSA1','EGWSA2','EGWSA3')

ttime = toc;

disp(['>>> Elapsed time: ', num2str(ttime),' seconds. <<<'])

msgbox('Completed TWS data analysis successfully!');

fclose(fstlog);


% --- Executes on button press in pbValidateByBasinDGMWObs.
function pbValidateByBasinDGMWObs_Callback(hObject, ~, handles)
% hObject    handle to pbValidateByBasinDGMWObs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pbPlotImages.
function pbPlotImages_Callback(hObject, eventdata, handles)
% hObject    handle to pbPlotImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% record the start time
tic;

%% get meta data
row = str2double(get(handles.txtRow, 'String'));    % row number of the GRACE TWS data
col = str2double(get(handles.txtCol, 'String'));    % col number of the GRACE TWS data
band = str2double(get(handles.txtBand, 'String'));  % the number of months
total_month = band;

years = 2002:2016;
idx = [1,7,18,30,42,54,66,78,90,102,112,122,131,140,149,158];

daily_files = get(handles.popFileType, 'Value');

% %% read in the JPL mascon id mask data
% % get JPL mascon id file name for Canada LCC grid
% mask_file = get(handles.txtGraceTwsIdMaskFile, 'String');
% if isempty(mask_file)
%     msgbox('Specify the JPL mascon id file for Canada LCC grid.');
%     return;
% end
% fid = fopen(mask_file);
%     mask = (fread(fid,[col, row],'float32')');
% fclose(fid);

%% get all EALCO daily or monthly twsa file names
ttws_folder = get(handles.txtAllEalcoTwsFile, 'String');
if isempty(ttws_folder)
    msgbox('Specify the EALCO Daily TWSA data file folder.');
    return;
end
allfnames = dir(ttws_folder);
ttws_fnames = {allfnames(~[allfnames.isdir]).name};
n_ttws = length(ttws_fnames);
%% read in all EALCO TWSA daily or monthly data
ttwsa = zeros(row,col,n_ttws);
time_list = zeros(n_ttws,1);
time_str_list = cell(n_ttws,1);
for i = 1:n_ttws
    ttwsi_file = [ttws_folder '\' ttws_fnames{i}];
    fid = fopen([ttws_folder '\' ttws_fnames{i}]);
    ttwsa(:,:,i) = (fread(fid,[col, row],'float32')');
    fclose(fid);
    if daily_files == 1
        [time_list(i), time_str_list{i}] = daysFromFileName(ttwsi_file,1);
    else
        [time_list(i), time_str_list{i}] = daysFromFileName(ttwsi_file,3);
    end
end

%% get all original GRACE TWSA data file names
grtwsa_folder = get(handles.txtAllGraceTwsFile, 'String');
if isempty(grtwsa_folder)
    msgbox('Specify the monthly GRACE TWSA data file folder.');
    return;
end
allfnames = dir(grtwsa_folder);
grtwsa_fnames = {allfnames(~[allfnames.isdir]).name}; 

%% get all spatially downscaled TWSA data file names
sdtwsa_folder = get(handles.txtAllSpatialDownscaledTwsFile, 'String');
if isempty(sdtwsa_folder)
    msgbox('Specify the saptially downscaled monthly TWSA data file folder.');
    return;
end
allfnames = dir(sdtwsa_folder);
sdtwsa_fnames = {allfnames(~[allfnames.isdir]).name}; 

%% get all spatiotemporally downscaled TWSA data file names
stdtwsa_folder = get(handles.txtAllTemporalDownscaledTwsFile, 'String');
if isempty(stdtwsa_folder)
    msgbox('Specify the saptiotemporally downscaled Daily TWSA data file folder.');
    return;
end
allfnames = dir(stdtwsa_folder);
stdtwsa_fnames = {allfnames(~[allfnames.isdir]).name}; 

%% get all spatiotemporally estimated GSWSA data file names
gwsa_folder = get(handles.txtAllTemporalDownscaledGwsFile, 'String');
if isempty(gwsa_folder)
    msgbox('Specify the Estimated Daily GWSA data file folder.');
    return;
end
allfnames = dir(gwsa_folder);
gwsa_fnames = {allfnames(~[allfnames.isdir]).name};
n_gwsa = length(gwsa_fnames);

%% get the output file folder
output_path = get(handles.txtOutputFile, 'String');
if isempty(output_path)
    msgbox('Specify the output data folder.');
    return;
end
   
%% read all original monthly GRACE TWSA data
n_grtwsa = length(grtwsa_fnames);
grtwsa = zeros(row,col,n_grtwsa);
grtwsa_std = zeros(row,col,n_grtwsa);
for i = 1:n_grtwsa
    grtwsa_file = [grtwsa_folder '\' grtwsa_fnames{i}];
    fid = fopen(grtwsa_file);
        grtwsa(:,:,i) = (fread(fid,[col, row],'float32')');
        grtwsa_std(:,:,i) = (fread(fid,[col, row],'float32')');
    fclose(fid);
    clear temp;
end 

%% read all spatially downscaled TWSA data
%     % get all saptially downscaled and assimilated TWSA file name
%     all_sdtwsa_file = get(handles.txtAllSpatialDownscaledTwsFile, 'String');
%     if isempty(all_sdtwsa_file)
%         msgbox('Specify the all saptially downscaled and assimilated TWSA file.');
%         return;
%     end
%     sdtwsa = zeros(row,col,band);
%     fid = fopen(all_sdtwsa_file);
%     for b=1:band
%         sdtwsa(:,:,b) = (fread(fid,[col, row],'float32')');
%     end
%     fclose(fid);

%     %% read all spatially estimated GSWSA data
%     % get all saptially downscaled and assimilated TWSA file name
%     all_sdgwsa_file = get(handles.txtAllSpatialDownscaledGswsFile, 'String');
%     if isempty(all_sdgwsa_file)
%         msgbox('Specify the all saptially downscaled and assimilated TWSA file.');
%         return;
%     end
%     sdgwsa = zeros(row,col,band);
%     fid = fopen(all_sdgwsa_file);
%     for b=1:band
%         sdgwsa(:,:,b) = (fread(fid,[col, row],'float32')');
%     end
%     fclose(fid);
    

%% specify the pixels range for the selected test region Alberta
% xbound = 500:780;
% ybound = 140:370;
xbound = 1:row;
ybound = 1:col;

%% specify the comparable scale range for all plots of TWSA
cmin = -250;
cmax = 250;
%% specify the comparable scale range for all plots of TWSA Std
cmin_std = 0;
cmax_std = 100;

%% remove the gray background color for paper figures
set(0,'DefaultTextInterpreter','none');
set(gcf, 'InvertHardCopy', 'off');

% %% plot Figure 1a of the paper
% gtwsa1 = grtwsa(xbound,ybound,1);
% ttwsa1 = ttwsa(xbound,ybound,1);
% % sdtwsa1 = sdtwsa(xbound,ybound,1);
% % sdgwsa1 = sdgwsa(xbound,ybound,1);
% % read the twsa and gswsa from the data files
% gwsi_file = [gwsa_folder '\' gwsa_fnames{1}];
% std_twsi_file = [stdtwsa_folder '\' stdtwsa_fnames{1}];        
% fid = fopen(std_twsi_file);
%     sdtwsa1 = (fread(fid,[col, row],'float32')');
% fclose(fid);
% fid = fopen(gwsi_file);
%     sdgwsa1 = (fread(fid,[col, row],'float32')');
% fclose(fid);
% sdtwsa1 = sdtwsa1(xbound,ybound);
% sdgwsa1 = sdgwsa1(xbound,ybound);
% 
% nodata_idx = find(sdtwsa1 == -32760);
% % gtwsa1(nodata_idx) = 0;
% % ttwsa1(nodata_idx) = 0;
% % sdtwsa1(nodata_idx) = 0;
% % sdgwsa1(nodata_idx) = 0;
% % cmin = min([min(grtwsa1(:)),min(ttwsa1(:)),min(sdtwsa1(:)),min(sdgwsa1(:))]);
% % cmax = max([max(grtwsa1(:)),max(ttwsa1(:)),max(sdtwsa1(:)),max(sdgwsa1(:))]);
% gtwsa1(nodata_idx)= cmin;
% ttwsa1(nodata_idx)= cmin;
% sdtwsa1(nodata_idx)= cmin;
% sdgwsa1(nodata_idx)= cmin;
% clims = [cmin,cmax];
% 
% %% set the no data area as transparent
% [r,c] = size(gtwsa1);
% A = ones(r,c); A(nodata_idx)= 0;
% 
% figure(1);
% set(gcf,'color','w');
% subplot(2,2,1)
% imagesc(gtwsa1,clims)
% alpha(A);
% set(gca,'XTick',[], 'YTick', []);
% colorbar
% subplot(2,2,2)
% imagesc(ttwsa1,clims)
% alpha(A);
% set(gca,'XTick',[], 'YTick', []);
% colorbar
% subplot(2,2,3)
% imagesc(sdtwsa1,clims)
% alpha(A);
% set(gca,'XTick',[], 'YTick', []);
% colorbar
% subplot(2,2,4)
% imagesc(sdgwsa1,clims)
% alpha(A);
% set(gca,'XTick',[], 'YTick', []);
% colorbar
% saveas(figure(1),[output_path 'Figure1a.png']);
% 
% %% plot Figure 1b of the paper
% gtwsa1 = grtwsa(xbound,ybound,end);
% ttwsa1 = ttwsa(xbound,ybound,end);
% % sdtwsa1 = sdtwsa(xbound,ybound,end);
% % sdgwsa1 = sdgwsa(xbound,ybound,end);
% % read the twsa and gswsa from files
% gwsi_file = [gwsa_folder '\' gwsa_fnames{end}];
% std_twsi_file = [stdtwsa_folder '\' stdtwsa_fnames{end}];        
% fid = fopen(std_twsi_file);
%     sdtwsa1 = (fread(fid,[col, row],'float32')');
% fclose(fid);
% fid = fopen(gwsi_file);
%     sdgwsa1 = (fread(fid,[col, row],'float32')');
% fclose(fid);
% sdtwsa1 = sdtwsa1(xbound,ybound);
% sdgwsa1 = sdgwsa1(xbound,ybound);
% 
% nodata_idx = find(sdtwsa1 == -32760);
% % gtwsa1(nodata_idx) = 0;
% % ttwsa1(nodata_idx) = 0;
% % sdtwsa1(nodata_idx) = 0;
% % sdgwsa1(nodata_idx) = 0;
% % cmin = min([min(grtwsa1(:)),min(ttwsa1(:)),min(sdtwsa1(:)),min(sdgwsa1(:))]);
% % cmax = max([max(grtwsa1(:)),max(ttwsa1(:)),max(sdtwsa1(:)),max(sdgwsa1(:))]);
% gtwsa1(nodata_idx)= cmin;
% ttwsa1(nodata_idx)= cmin;
% sdtwsa1(nodata_idx)= cmin;
% sdgwsa1(nodata_idx)= cmin;
% clims = [cmin,cmax];
% 
% [r,c] = size(gtwsa1);
% A = ones(r,c); A(nodata_idx)= 0;
% 
% figure(2);
% set(gcf,'color','w');
% subplot(2,2,1)
% imagesc(gtwsa1,clims);
% alpha(A);
% set(gca,'XTick',[], 'YTick', []);
% colorbar
% subplot(2,2,2)
% imagesc(ttwsa1,clims);
% alpha(A);
% set(gca,'XTick',[], 'YTick', []);
% colorbar
% subplot(2,2,3)
% imagesc(sdtwsa1,clims);
% alpha(A);
% set(gca,'XTick',[], 'YTick', []);
% colorbar
% subplot(2,2,4)
% imagesc(sdgwsa1,clims);
% alpha(A);
% set(gca,'XTick',[], 'YTick', []);
% colorbar
% saveas(figure(2),[output_path 'Figure1b.png']);

%% plot the animated monthly and annual data
if daily_files == 2  
    %% plot monthly image in an animated gif file
    filename1 = [output_path 'ttwsa_monthly.gif'];
    filename2 = [output_path 'sdtwsa_monthly.gif'];
    filename3 = [output_path 'sdgwsa_monthly.gif'];
    filename4 = [output_path 'grtwsa_monthly.gif'];
    filename5 = [output_path 'grtwsa_std_monthly.gif'];
    for i = 1:total_month
        ttwsa1 = ttwsa(xbound,ybound,i);
        sdtwsa1 = sdtwsa(xbound,ybound,i);
        sdgwsa1 = sdgwsa(xbound,ybound,i);
        grtwsa1 = grtwsa(xbound,ybound,i);
        grtwsa_std1 = grtwsa_std(xbound,ybound,i);
        
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
    end
    
    %% plot annual changes
    filename1 = [output_path 'ttwsa_yearly.gif'];
    filename2 = [output_path 'sdtwsa_yearly.gif'];
    filename3 = [output_path 'sdgwsa_yearly.gif'];
    filename4 = [output_path 'grtwsa_yearly.gif'];
    for i=1:length(years)
        gtws_annual = ones(length(xbound), length(ybound))*-32760;
        ttws_annual = ones(length(xbound), length(ybound))*-32760;
        atws_annual = ones(length(xbound), length(ybound))*-32760;
        gsws_annual = ones(length(xbound), length(ybound))*-32760;

        gtws_monthly = grtwsa(xbound,ybound,idx(i):idx(i+1));
        ttws_monthly = ttwsa(xbound,ybound,idx(i):idx(i+1));
        atws_monthly = sdtwsa(xbound,ybound,idx(i):idx(i+1));
        gsws_monthly = sdgwsa(xbound,ybound,idx(i):idx(i+1));
        
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

        figure(6);
        set(gcf,'color','w');
        imagesc(ttws_annual,clims);
        alpha(A);
        set(gca,'XTick',[], 'YTick', []);
        title(['EALCO TWSA: ' num2str(years(i))]);
        colorbar
        drawnow
        frame = getframe(6);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if i == 1
              imwrite(imind,cm,filename1,'gif', 'Loopcount',inf);
        else
              imwrite(imind,cm,filename1,'gif','WriteMode','append');
        end

        figure(7);
        set(gcf,'color','w');
        imagesc(atws_annual,clims)
        alpha(A);
        set(gca,'XTick',[], 'YTick', []);
        title(['SD TWSA: ' num2str(years(i))]);
        colorbar
        drawnow
        frame = getframe(7);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if i == 1
              imwrite(imind,cm,filename2,'gif', 'Loopcount',inf);
        else
              imwrite(imind,cm,filename2,'gif','WriteMode','append');
        end

        figure(8);
        set(gcf,'color','w');
        imagesc(gsws_annual,clims)
        alpha(A);
        set(gca,'XTick',[], 'YTick', []);
        title(['EGSWSA: ' num2str(years(i))]);
        colorbar
        drawnow
        frame = getframe(8);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if i == 1
              imwrite(imind,cm,filename3,'gif', 'Loopcount',inf);
        else
              imwrite(imind,cm,filename3,'gif','WriteMode','append');
        end
        
        figure(9)
        set(gcf,'color','w');
        imagesc(gtws_annual,clims);
        alpha(A);
        title(['GRACE TWSA: ' num2str(years(i))]);
        set(gca,'XTick',[], 'YTick', []);
        colorbar
        drawnow
        frame = getframe(9);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if i == 1
              imwrite(imind,cm,filename4,'gif', 'Loopcount',inf);
        else
              imwrite(imind,cm,filename4,'gif','WriteMode','append');
        end
    end
    
else %% Plot animated daily images
    filename1 = [output_path 'ttwsa_daily.gif'];
    filename2 = [output_path 'stdtwsa_daily.gif'];
    filename3 = [output_path 'stdgwsa_daily.gif'];
    filename4 = [output_path 'stdtwsa_std_daily.gif'];
    filename5 = [output_path 'stdgwsa_std_daily.gif'];
    for i = 1:n_gwsa
        gwsi_file = [gwsa_folder '\' gwsa_fnames{i}];
        std_twsi_file = [stdtwsa_folder '\' stdtwsa_fnames{i}];        
        [time_day, time_str] = daysFromFileName(gwsi_file,1);
        [tf, idx] = ismember(time_str,time_str_list);
        if tf
            time_str = convertDaysToDateTime(time_day);
            time_str = datestr(time_str, 'yyyy-mm-dd');
            % read the std twsa and uncertainties
            fid = fopen(std_twsi_file);
                stdtwsa = (fread(fid,[col, row],'float32')');
                stdtwsa_std = (fread(fid,[col, row],'float32')');
            fclose(fid);
            % read the gswsa and uncertainties
            fid = fopen(gwsi_file);
                stdgwsa = (fread(fid,[col, row],'float32')');
                stdgwsa_std = (fread(fid,[col, row],'float32')');
            fclose(fid);

            ttwsa1 = ttwsa(xbound,ybound,idx);
            stdtwsa1 = stdtwsa(xbound,ybound);
            stdgwsa1 = stdgwsa(xbound,ybound);
            stdtwsa_std1 = stdtwsa_std(xbound,ybound);
            stdgwsa_std1 = stdgwsa_std(xbound,ybound);

            nodata_idx = find(stdtwsa1 == -32760);
            % ttwsa1(nodata_idx) = 0;
            % stdtwsa1(nodata_idx) = 0;
            % stdgwsa1(nodata_idx) = 0;
            % cmin = min([min(gtwsa1(:)),min(ttwsa1(:)),min(stdtwsa1(:)),min(stdgwsa1(:))]);
            % cmax = max([max(gtwsa1(:)),max(ttwsa1(:)),max(stdtwsa1(:)),max(stdgwsa1(:))]);
            ttwsa1(nodata_idx)= cmin;
            stdtwsa1(nodata_idx)= cmin;
            stdgwsa1(nodata_idx)= cmin;
            clims = [cmin,cmax];
            stdtwsa_std1(nodata_idx)= cmax_std/2;
            stdgwsa_std1(nodata_idx)= cmax_std/2;
            clims_std = [cmin_std,cmax_std];
            
            [r,c] = size(stdtwsa1);
            A = ones(r,c); A(nodata_idx)= 0;

            if ~isempty(find(stdtwsa1 > 0, 1))
                figure(1);
                set(gcf,'color','w');
                imagesc(ttwsa1,clims);
                alpha(A);
                set(gca,'XTick',[], 'YTick', []);
                title(['EALCO TWSA: ' time_str]);
                colorbar
                drawnow
                frame = getframe(1);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                if i == 1
                      imwrite(imind,cm,filename1,'gif', 'Loopcount',inf);
                else
                      imwrite(imind,cm,filename1,'gif','WriteMode','append');
                end

                figure(2);
                set(gcf,'color','w');
                imagesc(stdtwsa1,clims);
                alpha(A);
                set(gca,'XTick',[], 'YTick', []);
                title(['STD TWSA: ' time_str]);
                colorbar
                drawnow
                frame = getframe(2);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                if i == 1
                      imwrite(imind,cm,filename2,'gif', 'Loopcount',inf);
                else
                      imwrite(imind,cm,filename2,'gif','WriteMode','append');
                end

                figure(3);
                set(gcf,'color','w');
                imagesc(stdgwsa1,clims);
                alpha(A);
                set(gca,'XTick',[], 'YTick', []);
                title(['EGSWSA: ' time_str]);
                colorbar
                drawnow
                frame = getframe(3);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                if i == 1
                      imwrite(imind,cm,filename3,'gif', 'Loopcount',inf);
                else
                      imwrite(imind,cm,filename3,'gif','WriteMode','append');
                end
                
                figure(4);
                set(gcf,'color','w');
                imagesc(stdtwsa_std1,clims_std);
                alpha(A);
                set(gca,'XTick',[], 'YTick', []);
                title(['STD TWSA Std: ' time_str]);
                colorbar
                drawnow
                frame = getframe(4);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                if i == 1
                      imwrite(imind,cm,filename4,'gif', 'Loopcount',inf);
                else
                      imwrite(imind,cm,filename4,'gif','WriteMode','append');
                end

                figure(5);
                set(gcf,'color','w');
                imagesc(stdgwsa_std1,clims_std);
                alpha(A);
                set(gca,'XTick',[], 'YTick', []);
                title(['EGSWSA Std: ' time_str]);
                colorbar
                drawnow
                frame = getframe(5);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                if i == 1
                      imwrite(imind,cm,filename5,'gif', 'Loopcount',inf);
                else
                      imwrite(imind,cm,filename5,'gif','WriteMode','append');
                end
            end
        end
    end
end



function edit29_Callback(hObject, eventdata, ~)
% hObject    handle to txtGWLatLonFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtGWLatLonFile as text
%        str2double(get(hObject,'String')) returns contents of txtGWLatLonFile as a double


% --- Executes during object creation, after setting all properties.
function edit29_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtGWLatLonFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowseGMWPostionFile.
function pbBrowseGMWPostionFile_Callback(hObject, eventdata, ~)
% hObject    handle to pbBrowseGMWPostionFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile( ...
    {  '*.*','Specify a Groundwater Well Lat Lon position file (*.*)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Specify a Groundwater Well Lat Lon position file');
 if  filename > 0  
    set(handles.txtGWLatLonFile, 'String', [pathname filename]);
 end 

% --- Executes on selection change in popGMWObsType.
function popGMWObsType_Callback(hObject, eventdata, handles)
% hObject    handle to popGMWObsType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popGMWObsType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popGMWObsType


% --- Executes during object creation, after setting all properties.
function popGMWObsType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popGMWObsType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
