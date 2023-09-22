function varargout = snp_genodropbox(varargin)
% SNP_GENODROPBOX M-file for snp_genodropbox.fig
%      SNP_GENODROPBOX, by itself, creates a new SNP_GENODROPBOX or raises the existing
%      singleton*.
%
%      H = SNP_GENODROPBOX returns the handle to a new SNP_GENODROPBOX or the handle to
%      the existing singleton*.
%
%      SNP_GENODROPBOX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SNP_GENODROPBOX.M with the given input arguments.
%
%      SNP_GENODROPBOX('Property','Value',...) creates a new SNP_GENODROPBOX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before snp_genodropbox_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to snp_genodropbox_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help snp_genodropbox

% Last Modified by GUIDE v2.5 15-Apr-2008 08:41:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name', mfilename, ...
    'gui_Singleton', gui_Singleton, ...
    'gui_OpeningFcn', @snp_genodropbox_OpeningFcn, ...
    'gui_OutputFcn', @snp_genodropbox_OutputFcn, ...
    'gui_LayoutFcn', [], ...
    'gui_Callback', []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before snp_genodropbox is made visible.
    function snp_genodropbox_OpeningFcn(hObject, eventdata, handles, varargin)
        % This function has no output args, see OutputFcn.
        % hObject    handle to figure
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)
        % varargin   command line arguments to snp_genodropbox (see VARARGIN)

        % Choose default command line output for snp_genodropbox
        handles.output = hObject;

        % Update handles structure
        guidata(hObject, handles);

        % UIWAIT makes snp_genodropbox wait for user response (see UIRESUME)
        % uiwait(handles.figure1);


        % --- Outputs from this function are returned to the command line.
            function varargout = snp_genodropbox_OutputFcn(hObject, eventdata, handles)
                % varargout  cell array for returning output args (see VARARGOUT);
                % hObject    handle to figure
                % eventdata  reserved - to be defined in a future version of MATLAB
                % handles    structure with handles and user data (see GUIDATA)

                % Get default command line output from handles structure
                varargout{1} = handles.output;


                    function edit1_Callback(hObject, eventdata, handles)
                        % hObject    handle to edit1 (see GCBO)
                        % eventdata  reserved - to be defined in a future version of MATLAB
                        % handles    structure with handles and user data (see GUIDATA)

                        % Hints: get(hObject,'String') returns contents of edit1 as text
                        %        str2double(get(hObject,'String')) returns contents of edit1 as a double


                        % --- Executes during object creation, after setting all properties.
                            function edit1_CreateFcn(hObject, eventdata, handles)
                                % hObject    handle to edit1 (see GCBO)
                                % eventdata  reserved - to be defined in a future version of MATLAB
                                % handles    empty - handles not created until after all CreateFcns called

                                % Hint: edit controls usually have a white background on Windows.
                                %       See ISPC and COMPUTER.
                                if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
                                    set(hObject, 'BackgroundColor', 'white');
                                end
