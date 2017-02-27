function varargout = dprstable(varargin)
%DPRSTABLE - dprstable table GUI
% dprstable M-file for dprstable.fig
%      dprstable, by itself, creates ds new dprstable or raises the existing
%      singleton*.
%
%      H = dprstable returns the handle to ds new dprstable or the handle to
%      the existing singleton*.
%
%      dprstable('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in dprstable.M with the given input arguments.
%
%      dprstable('Property','Value',...) creates ds new dprstable or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dprstable_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dprstable_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

% Edit the above text to modify the response to help dprstable

% Last Modified by GUIDE v2.5 25-Jun-2010 18:35:41
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dprstable_OpeningFcn, ...
                   'gui_OutputFcn',  @dprstable_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before dprstable is made visible.
function dprstable_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in ds future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dprstable (see VARARGIN)

% Choose default command line output for dprstable
handles.output = hObject;

scrsz=get(0,'ScreenSize');
pos_act=get(gcf,'Position');
xr=scrsz(3)-pos_act(3);
xp=round(xr/2);
yr=scrsz(4)-pos_act(4);
yp=round(yr/2);
set(gcf,'position',[xp yp pos_act(3) pos_act(4)]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dprstable wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dprstable_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in ds future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function ds_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ds (see GCBO)
% eventdata  reserved - to be defined in ds future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have ds white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function ds_Callback(hObject, eventdata, handles)
% hObject    handle to ds (see GCBO)
% eventdata  reserved - to be defined in ds future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ds as text
%        str2double(get(hObject,'String')) returns contents of ds as ds double
updatenum(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function ps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ps (see GCBO)
% eventdata  reserved - to be defined in ds future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have ds white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function ps_Callback(hObject, eventdata, handles)
% hObject    handle to ps (see GCBO)
% eventdata  reserved - to be defined in ds future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ps as text
%        str2double(get(hObject,'String')) returns contents of ps as ds double
updatenum(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function dr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dr (see GCBO)
% eventdata  reserved - to be defined in ds future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have ds white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function dr_Callback(hObject, eventdata, handles)
% hObject    handle to dr (see GCBO)
% eventdata  reserved - to be defined in ds future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dr as text
%        str2double(get(hObject,'String')) returns contents of dr as ds double
updatenum(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function pr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pr (see GCBO)
% eventdata  reserved - to be defined in ds future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have ds white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function pr_Callback(hObject, eventdata, handles)
% hObject    handle to pr (see GCBO)
% eventdata  reserved - to be defined in ds future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pr as text
%        str2double(get(hObject,'String')) returns contents of pr as ds double
updatenum(hObject, eventdata, handles);

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in ds future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Close.
function Close_Callback(hObject, eventdata, handles)
% hObject    handle to Close (see GCBO)
% eventdata  reserved - to be defined in ds future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close;



% --- Executes during object creation, after setting all properties.
function dsps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dsps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
%if ispc
%    set(hObject,'BackgroundColor','gray');
%else
%    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
%end

function dsps_Callback(hObject, eventdata, handles)
% hObject    handle to dsps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dsps as text
%        str2double(get(hObject,'String')) returns contents of dsps as a double



% --- Executes on button press in fisher.
function fisher_Callback(hObject, eventdata, handles)
% hObject    handle to fisher (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[ds,ps,dr,pr]=updatenum(hObject, eventdata, handles);
[P]=fisherextest(round(ds),round(ps),round(dr),round(pr));
i_dispheader('Fisher''s Exact Test')
fprintf(['   P-value = %g %s\n'], P, sigtag(P));
i_dispfooter


% --- Executes on button press in chisquare.
function chisquare_Callback(hObject, eventdata, handles)
% hObject    handle to chisquare (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[ds,ps,dr,pr]=updatenum(hObject, eventdata, handles);
[P,P2,X2,X2C]=chi2test(ds,ps,dr,pr);
i_dispheader('Chi-square Test')
fprintf(['Chi-square value: %10.3f\n'], X2);
fprintf(['      P-value = %g %s\n'], P, sigtag(P));
fprintf(['Chi-square with Yates'' correction: %10.3f\n'], X2C);
fprintf(['      P-value = %g %s\n'], P2,sigtag(P2));
%dr/ds
%pr/ps
fprintf(['      alpha = 1-(ds*pr)/(dr*ps) = %10.5f\n'], 1-(ds*pr)/(dr*ps));
i_dispfooter




% --- Executes on button press in gtest.
function gtest_Callback(hObject, eventdata, handles)
% hObject    handle to gtest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[ds,ps,dr,pr]=updatenum(hObject, eventdata, handles);
[P,G]=gtest(ds,ps,dr,pr);
i_dispheader('G Test')
fprintf(['G value: %g\n'], G);
fprintf(['   P-value = %g %s\n'], P, sigtag(P));
i_dispfooter




% --- Executes during object creation, after setting all properties.
function drpr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to drpr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
%if ispc
%    set(hObject,'BackgroundColor','white');
%else
%    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
%end



function drpr_Callback(hObject, eventdata, handles)
% hObject    handle to drpr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of drpr as text
%        str2double(get(hObject,'String')) returns contents of drpr as a double


% --- Executes during object creation, after setting all properties.
function total_CreateFcn(hObject, eventdata, handles)
% hObject    handle to total (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
%if ispc
%    set(hObject,'BackgroundColor','white');
%else
%    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
%end



function total_Callback(hObject, eventdata, handles)
% hObject    handle to total (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of total as text
%        str2double(get(hObject,'String')) returns contents of total as a double


% --- Executes during object creation, after setting all properties.
function dsdr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dsdr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
%if ispc
%    set(hObject,'BackgroundColor','white');
%else
%    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
%end



function dsdr_Callback(hObject, eventdata, handles)
% hObject    handle to dsdr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dsdr as text
%        str2double(get(hObject,'String')) returns contents of dsdr as a double


% --- Executes during object creation, after setting all properties.
function pspr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pspr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
%if ispc
%    set(hObject,'BackgroundColor','white');
%else
%    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
%end

function figure1_CreateFcn(hObject, eventdata, handles)



function pspr_Callback(hObject, eventdata, handles)
% hObject    handle to pspr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pspr as text
%        str2double(get(hObject,'String')) returns contents of pspr as a double



function [ds,ps,dr,pr]=updatenum(hObject, eventdata, handles)
ps=str2double(get(handles.ps, 'String'));
pr=str2double(get(handles.pr, 'String'));
ds=str2double(get(handles.ds, 'String'));
dr=str2double(get(handles.dr, 'String'));

set(handles.dsps,'String',num2str(ds+ps));
set(handles.drpr,'String',num2str(dr+pr));
set(handles.dsdr,'String',num2str(ds+dr));
set(handles.pspr,'String',num2str(ps+pr));
set(handles.total,'String',num2str(ps+pr+ds+dr));
