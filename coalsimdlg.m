function varargout = coalsimdlg(varargin)
%COALSIMDLG M-file for coalsimdlg.fig
%      COALSIMDLG, by itself, creates a new COALSIMDLG or raises the existing
%      singleton*.
%
%      H = COALSIMDLG returns the handle to a new COALSIMDLG or the handle to
%      the existing singleton*.
%
%      COALSIMDLG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COALSIMDLG.M with the given input arguments.
%
%      COALSIMDLG('Property','Value',...) creates a new COALSIMDLG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before coalsimdlg_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to coalsimdlg_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help coalsimdlg

% Last Modified by GUIDE v2.5 23-Feb-2007 18:37:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @coalsimdlg_OpeningFcn, ...
                   'gui_OutputFcn',  @coalsimdlg_OutputFcn, ...
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


% --- Executes just before coalsimdlg is made visible.
function coalsimdlg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to coalsimdlg (see VARARGIN)

% Choose default command line output for coalsimdlg
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

% UIWAIT makes coalsimdlg wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = coalsimdlg_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in thetacheck.
function thetacheck_Callback(hObject, eventdata, handles)
% hObject    handle to thetacheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of thetacheck
if(get(handles.segsitescheck,'Value'))
    set(handles.segsitescheck,'Value',0);
    set(handles.segsites,'Enable','off');
end
set(handles.theta,'Enable','on');


% --- Executes on button press in segsitescheck.
function segsitescheck_Callback(hObject, eventdata, handles)
% hObject    handle to segsitescheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of segsitescheck
if(get(handles.thetacheck,'Value'))
    set(handles.thetacheck,'Value',0);
    set(handles.theta,'Enable','off');
end
set(handles.segsites,'Enable','on');





% --- Executes during object creation, after setting all properties.
function theta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function theta_Callback(hObject, eventdata, handles)
% hObject    handle to theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of theta as text
%        str2double(get(hObject,'String')) returns contents of theta as a double


% --- Executes during object creation, after setting all properties.
function segsites_CreateFcn(hObject, eventdata, handles)
% hObject    handle to segsites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function segsites_Callback(hObject, eventdata, handles)
% hObject    handle to segsites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of segsites as text
%        str2double(get(hObject,'String')) returns contents of segsites as a double


% --- Executes during object creation, after setting all properties.
function samplesize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to samplesize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function samplesize_Callback(hObject, eventdata, handles)
% hObject    handle to samplesize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of samplesize as text
%        str2double(get(hObject,'String')) returns contents of samplesize as a double


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nsam=str2num(get(handles.samplesize,'String'));
nrep=str2num(get(handles.nreplicate,'String'));

if (isempty(nsam)), warndlg('Sample size must be a number.'); return; end
if (isempty(nrep)), warndlg('Number of replicates must be a number.'); return; end

if (nsam==0), warndlg('Sample size cannot be zero.'); return; end
if (nrep==0), warndlg('Number of replicates cannot be zero.'); return; end

if(get(handles.thetacheck,'Value')),
    theta=str2num(get(handles.theta,'String'));
    segs=0;
else
    segs=str2num(get(handles.segsites,'String'));
    theta=0;
end


if(get(handles.norecombinationcheck,'Value')),
     rho=0.0;
     nsites=2;
     %warndlg('No recomibination.');
else
    rho=str2num(get(handles.recombinationrho,'String'));
    nsites=str2num(get(handles.nsites,'String'));
    if (rho<0 | rho>500), warndlg('4Nr should between 0 - 500.'); return; end
    if (nsites<2), warndlg('nsites should >=2.'); return; end
end



if ~((segs>0)||(theta>0)),
    warndlg('Either THETA or SEGS has to be specified.');
    return;
end


switch (get(handles.statistic,'Value'))
    case (1)
     	[x]=thetapi_simu(nsam,nrep,theta,segs,rho,nsites);
        i_showreport(x,'Nucleotide Diversity, Pi')
    case (2)   %
        [x]=thetaw_simu(nsam,nrep,theta,segs,rho,nsites);
        i_showreport(x,'Thet-W (from S)')
    case (3)   % tajima's d
        [x]=tajima89d_simu(nsam,nrep,theta,segs,rho,nsites);
        i_showreport(x,'Tajima''s D');
    case (4)
        [x]=faywu00h_simu(nsam,nrep,theta,segs,rho,nsites);
        i_showreport(x,'Fay Wu''s H');
    case (5)
     	[x]=thetah_simu(nsam,nrep,theta,segs,rho,nsites);
        i_showreport(x,'Fay''s thetaH');
    case (6)
        [x,y]=fuli93dsfs_simu(nsam,nrep,theta,segs,rho,nsites);
        i_showreport(x,'Fu and Li''s D*');
        i_showreport(y,'Fu and Li''s F*');
    case (7)
     	[x]=hapdiv_simu(nsam,nrep,theta,segs,rho,nsites);
        i_showreport(x,'Haplotype diversity');

end



function i_showreport(x,txt)
        i_dispheader(txt)
        [lowci,upci]=cireport(x);
        i_dispfooter
        figure;
        hist(x,50);
        %h = findobj(gca,'Type','patch');
        %set(h,'FaceColor',[.3 .4 .4])
        xlabel(txt)
        ylabel('#')
        title('Coalescent Simulation')
        vline([lowci,upci])





% --- Executes on button press in Cancel.
function Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close

% --- Executes during object creation, after setting all properties.
function nreplicate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nreplicate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function nreplicate_Callback(hObject, eventdata, handles)
% hObject    handle to nreplicate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nreplicate as text
%        str2double(get(hObject,'String')) returns contents of nreplicate as a double


% --- Executes during object creation, after setting all properties.
function statistic_CreateFcn(hObject, eventdata, handles)
% hObject    handle to statistic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in statistic.
function statistic_Callback(hObject, eventdata, handles)
% hObject    handle to statistic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns statistic contents as cell array
%        contents{get(hObject,'Value')} returns selected item from statistic




% --- Executes on button press in NorecombinationCheck.
%function NorecombinationCheck_Callback(hObject, eventdata, handles)
% hObject    handle to NorecombinationCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of NorecombinationCheck



% --- Executes on button press in recombinationrhocheck.
function recombinationrhocheck_Callback(hObject, eventdata, handles)
% hObject    handle to recombinationrhocheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of recombinationrhocheck

if(get(handles.norecombinationcheck,'Value'))
    set(handles.norecombinationcheck,'Value',0);
end
set(handles.recombinationrho,'Enable','on');
set(handles.nsites,'Enable','on');



% --- Executes during object creation, after setting all properties.
function recombinationrho_CreateFcn(hObject, eventdata, handles)
% hObject    handle to recombinationrho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function recombinationrho_Callback(hObject, eventdata, handles)
% hObject    handle to recombinationrho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of recombinationrho as text
%        str2double(get(hObject,'String')) returns contents of recombinationrho as a double


% --- Executes on button press in norecombinationcheck.
function norecombinationcheck_Callback(hObject, eventdata, handles)
% hObject    handle to norecombinationcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of norecombinationcheck
if(get(handles.recombinationrhocheck,'Value'))
    set(handles.recombinationrhocheck,'Value',0);
end
set(handles.recombinationrho,'Enable','off');
set(handles.nsites,'Enable','off');


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5


% --- Executes during object creation, after setting all properties.
function nsites_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nsites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function nsites_Callback(hObject, eventdata, handles)
% hObject    handle to nsites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nsites as text
%        str2double(get(hObject,'String')) returns contents of nsites as a double


