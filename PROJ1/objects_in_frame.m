function varargout = objects_in_frame(varargin)
% OBJECTS_IN_FRAME MATLAB code for objects_in_frame.fig
%      OBJECTS_IN_FRAME, by itself, creates a new OBJECTS_IN_FRAME or raises the existing
%      singleton*.
%
%      H = OBJECTS_IN_FRAME returns the handle to a new OBJECTS_IN_FRAME or the handle to
%      the existing singleton*.
%
%      OBJECTS_IN_FRAME('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OBJECTS_IN_FRAME.M with the given input arguments.
%
%      OBJECTS_IN_FRAME('Property','Value',...) creates a new OBJECTS_IN_FRAME or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before objects_in_frame_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to objects_in_frame_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help objects_in_frame

% Last Modified by GUIDE v2.5 12-Jun-2015 15:49:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @objects_in_frame_OpeningFcn, ...
                   'gui_OutputFcn',  @objects_in_frame_OutputFcn, ...
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


% --- Executes just before objects_in_frame is made visible.
function objects_in_frame_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to objects_in_frame (see VARARGIN)

% Choose default command line output for objects_in_frame
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes objects_in_frame wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = objects_in_frame_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output=hObject;
[a b]=uigetfile({'*.*'});
img=imread([b a]);
grayy=rgb2gray(img);
gr=graythresh(grayy);
handles.bw=im2bw(grayy,gr);
imshow(img,'Parent',handles.axes1);
guidata(hObject,handles);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output=hObject;
inverse_binary=not(handles.bw);
[handles.L handles.Num_object]=bwlabel(inverse_binary);
set(handles.text2,'string',handles.Num_object);

imshow(handles.L,'Parent',handles.axes2);


guidata(hObject,handles);