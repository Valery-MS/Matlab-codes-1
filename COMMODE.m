function varargout = COMMODE(varargin)
% COMMODE MATLAB code for COMMODE.fig
%   COMMODE, by itself, creates a new COMMODE or raises the existing
%   singleton*.

%   H = COMMODE returns the handle to a new COMMODE or the handle to
%   the existing singleton*.
%
%   COMMODE('CALLBACK',hObject,eventData,handles,...) calls the local
%   function named CALLBACK in COMMODE.M with the given input arguments.

%   COMMODE('Property','Value',...) creates a new COMMODE or raises the
%   existing singleton*.  Starting from the left, property value pairs are
%   applied to the GUI before COMMODE_OpeningFcn gets called.  An
%   unrecognized property name or invalid value makes property application
%   stop.  All inputs are passed to COMMODE_OpeningFcn via varargin.

%   *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%   instance to run (singleton)".

% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help COMMODE
% Last Modified by GUIDE v2.5 05-Sep-2016 19:56:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @COMMODE_OpeningFcn, ...
                   'gui_OutputFcn',  @COMMODE_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout,  [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else         gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before COMMODE is made visible.
function COMMODE_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to COMMODE (see VARARGIN)

handles.output = hObject;  % Choose default command line output for COMMODE
guidata(hObject, handles); % Update handles structure

% UIWAIT makes COMMODE wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = COMMODE_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% Get default command line output from handles structure
varargout{1} = handles.output;

                      % USER FUNCTIONS
% Info about Equation_Name(Air, Bessel, Shrod, DNK) is in the 
% handles.Equation_Name.UserData in the Menu Editor

function Problem_Callback(hObject, eventdata, handles)
% hObject    handle to Problem (see GCBO)
% Hints: get(hObject,'String') returns contents of Problem as text
%     str2double(get(hObject,'String')) returns contents of hs as a double
handles.DataBestC.Visible = 'off';

function Linear_Callback(hObject, eventdata, handles)
global yTA
handles.LiNo.String = 'linear';
handles.Methods.String = {'Magnus (lin)','ds10L (lin)','dp853','ode113',...
                          'ode45','odex',};

function Nonlinear_Callback(hObject, eventdata, handles)
global yTA
handles.LiNo.String = 'nonlinear';
handles.Methods.String = {'ds10','dp853','ode113','ode45','odex'};

function Airy_Callback(hObject, eventdata, handles)
global omegAiry
load(hObject.UserData{1})     % cel-array GenInpD, SpecInpD = {omegAiry=1}
FillGenInpD(handles,GenInpD); % fill GenInpD 
                     % fill SpecInpD
omegAiry = SpecInpD{1};  handles.omegAiry.String = num2str(omegAiry);   

handles.nSDE.Enable = 'off';
handles.RelT.Enable = 'off';
handles.AbsT.Enable = 'off';
handles.y0.Enable   = 'off';  % y0 = (1, 0)

handles.DataAiry.Visible   = 'on';  % turn on Airy data
handles.DataBess.Visible   = 'off'; % turn off other data
handles.DataShrod.Visible  = 'off';
handles.DataDNK.Visible    = 'off';
handles.DataArbitr.Visible = 'off';

handles.EqNam.String     = hObject.Label;
handles.EqNam.UserData   = hObject.UserData;
handles.Linear.UserData  = 1;

% --------------------------------------------------------------------
function Bess_Callback(hObject, eventdata, handles)
global nBess
load(hObject.UserData{1})     % cel-array GenInpD, SpecInpD = {nBess=1}
FillGenInpD(handles,GenInpD); % fill GenInpD   
                  % fill SpecInpD
nBess = SpecInpD{1};  handles.nBess.String = num2str(nBess);

handles.nSDE.Enable = 'off';
handles.RelT.Enable = 'off';
handles.AbsT.Enable = 'off';
handles.y0.Enable   = 'off';

handles.DataBess.Visible   = 'on';  % turn on Bess data
handles.DataAiry.Visible   = 'off'; % turn off other data
handles.DataShrod.Visible  = 'off';
handles.DataDNK.Visible    = 'off';
handles.DataArbitr.Visible = 'off';

handles.EqNam.String     = hObject.Label;
handles.EqNam.UserData   = hObject.UserData;
handles.Linear.UserData  = 2;

% --------------------------------------------------------------------
function Shrod_Callback(hObject, eventdata, handles)
global Pamp Pfreq
load(hObject.UserData{1})     % cel-array GenInpD, SpecInpD = {Pamp Pfreq}
FillGenInpD(handles,GenInpD); % fill GenInpD   
               % fill SpecInpD
Pamp  = SpecInpD{1};  handles.PampShrod.String  = num2str(Pamp);
Pfreq = SpecInpD{2};  handles.PfreqShrod.String = num2str(Pfreq);

handles.nSDE.Enable = 'off';
handles.RelT.Enable = 'off';
handles.AbsT.Enable = 'off';

handles.DataShrod.Visible  = 'on';  % turn on Shrod data
handles.DataAiry.Visible   = 'off'; % turn off other data
handles.DataBess.Visible   = 'off';
handles.DataDNK.Visible    = 'off';
handles.DataArbitr.Visible = 'off';

handles.EqNam.String   = hObject.Label;
handles.EqNam.UserData = hObject.UserData;
handles.LiNo.UserData  = 3;

% --------------------------------------------------------------------
function Eq4_Callback(hObject, eventdata, handles)
handles.EqNam.String    = hObject.Label;
handles.Linear.UserData = 4;

% --------------------------------------------------------------------
function DNK_Callback(hObject, eventdata, handles)
global SpecInpD a1 a IV NEq_ cT_  GenInpD yTA

load(hObject.UserData{1})     % cel-array GenInpD, SpecInpD = {Pamp Pfreq}
FillGenInpD(handles,GenInpD); % fill GenInpD   
                    % fill SpecInpD
                 
a = SpecInpD{1}; a1 = 1e4*a; handles.a1.String  = num2str(a1,'%g  ');
IV  = SpecInpD{2};     handles.IV.String  = ...
      sprintf('%g  %g  %g  %g  %.0e  %g  %g %g',IV);
NEq_= SpecInpD{3};   handles.NEq_.String = num2str(NEq_,'%g  ');
cT_ = SpecInpD{4};   handles.cT_.String = sprintf('%.3g  ',cT_);
                     handles.DelFig.String = '0';
                       
handles.nSDE.Enable = 'off';
handles.y0.Enable   = 'on';

handles.DataDNK.Visible    = 'on';  % turn on DNK data
handles.DataAiry.Visible   = 'off'; % turn off other data
handles.DataBess.Visible   = 'off';
handles.DataShrod.Visible  = 'off';
handles.DataArbitr.Visible = 'off';

handles.EqNam.String   = hObject.Label;
handles.EqNam.UserData = hObject.UserData;

% --------------------------------------------------------------------
function Arbitr_Callback(hObject, eventdata, handles)
global t0 tn h RelT AbsT

load(hObject.UserData{1})
handles.DataAiry.Visible   = 'off';
handles.DataBess.Visible   = 'off';
handles.DataShrod.Visible  = 'off';
handles.DataDNK.Visible    = 'off';
handles.DataArbitr.Visible = 'on';

handles.EqNam.String   = hObject.Label;
handles.EqNam.UserData = hObject.UserData;

% --------------------------------------------------------------------
function Tasks_Callback(hObject, eventdata, handles)
if hObject.Value == 3
   if strcmp(handles.LiNo.String,'nonlinear')
      errordlg('This task runs only for linear equations');
      return
   end
   handles.DataBestC.Visible = 'on';
   handles.Methods.Enable    = 'off';
else
   handles.DataBestC.Visible = 'off';
   handles.Methods.Enable    = 'on';   
   end

% --- Executes when Methods is resized.
%function Methods_SizeChangedFcn(hObject, eventdata, handles)

% --------------------------------------------------------------------
function hs_Callback(hObject, eventdata, handles)
if isempty(str2num(hObject.String)), errordlg('hs are not numbers');end

function Methods_Callback(hObject, eventdata, handles)
% Hints:
% c = cellstr(get(hObject,'String')) returns Methods contents as cell array
% c{get(hObject,'Value')} returns selected item from Methods

function Tools_Callback(hObject, eventdata, handles)

% --- Executes on button press in Run.
function Run_Callback(hObject, eventdata, handles)
switch handles.Tasks.Value
   case 1, Task_TACh(handles);
   case 2
      if strcmp(handles.LiNo.String,'linear'), Task_S_C_lin(handles);
      else                                   Task_S_C_nolin(handles);end
   case 3, Task_BestC(handles);
   case 4, Task_Solve(handles);  % Here should be the Task_Solve program 
   case 5, Task_nDSvsDP(handles);
end

% Check of input Matrix & Saving in UserData if it's without rough mistakes
function a = CheckMat(hObject,m,n)
a = str2num(hObject.String);
if     isempty(a),          errordlg('It is not a number');
elseif ~all(size(a)==[m n]),errordlg(sprintf('Size ~= %dx%d',m,n)); a=[]; 
end

% Check of input Scalar and Saving if it's without rough mistakes
function s = CheckScal(hObject)
s = str2num(hObject.String);
if       isempty(s), errordlg('It is not a number');
elseif ~isscalar(s), errordlg('It is not a scalar'); s = [];
end

function a1_Callback(hObject, eventdata, handles)
global a1 SpecInpD
a_ = CheckMat(hObject,1,4);
if ~isempty(a_), a1 = a_;  SpecInpD{1} = {1e-4*a1}; end

function IV_Callback(hObject, eventdata, handles)
% IV - Initial Value
% IVI = IV(1) - ? Initial Value Index
% IV(4) - devider of pi: Step of InitialValue hIV = pi/IV(4)
% (IV(2):IV(3))*d - range of IV
global IV SpecInpD
IV_ = CheckMat(hObject,1,8);
if ~isempty(IV_),  IV = IV_; SpecInpD{2} = IV; end

function NEq__Callback(hObject, eventdata, handles)
% NEq_ = 1,2,3 - № уравнения (1-Airy, 2- Bessel, 3 - Shrod)
% m   = 1,2,3 - параметр, определяющий порядок РС
% Nh  = 1...6 - № шага h из массива hs. Для шага с этим номером N
%            вычислялись лучшие к-ты BestC{m}{N}
% ks - подмножество из 1:К (К>100) - номера наборов Lmn по возраст погрешн
%   например, ks=[1 3 5] => счёт для 3х наборов (La,mu,nu): 1го,3го и 5го
%   из массива BestC{m}{Nh}, причём решение у-я NEq р.схемой с 1-м набором 
%   к-тов имело мин погрешн, с 5-м набором - 5-ю по величине погрешность
global NEq_ SpecInpD
NEq__ = CheckMat(hObject,1,4);
if ~isempty(NEq__), NEq_ = NEq__; SpecInpD{3} = NEq_; end

function cT__Callback(hObject, eventdata, handles)
global cT_ SpecInpD
cT__ = CheckMat(hObject,1,3);
if ~isempty(cT__), cT_ = cT__; SpecInpD{7} = cT_;end


function t0_Callback(hObject, eventdata, handles)
global t0 GenInpD
t0_ = CheckScal(hObject);
if ~isempty(t0_), t0 = t0_; GenInpD{1} = t0; end

function tn_Callback(hObject, eventdata, handles)
global tn GenInpD
tn_ = CheckScal(hObject);
if ~isempty(tn_), tn = tn_; GenInpD{2} = tn; end

function h1_h2_Callback(hObject, eventdata, handles)
global h1_h2 GenInpD
h = str2num(hObject.String);
s = size(h);    
if s(1) == 1 && 1 <= s(2) && s(2) <= 2,  h1_h2 = h; GenInpD{3} = h1_h2;
else  errordlg(sprintf('Size(h1_h2) = %dx%d',s));  end  

function RelT_Callback(hObject, eventdata, handles)
global RelT GenInpD
RelT_ = CheckScal(hObject);
if ~isempty(RelT_), RelT = RelT_; GenInpD{4} = RelT;end

function AbsT_Callback(hObject, eventdata, handles)
global AbsT GenInpD
AbsT_ = CheckScal(hObject);
if ~isempty(AbsT_), AbsT = AbsT_; GenInpD{5} = AbsT;end

function nSDE_Callback(hObject, eventdata, handles)
global nSDE GenInpD
nSDE_ = CheckScal(hObject);
if ~isempty(nSDE_), nSDE = nSDE_; GenInpD{6} = nSDE;end

function y0_Callback(hObject, eventdata, handles)
global nSDE y0 GenInpD
y0_ = CheckMat(hObject,1,nSDE);
if ~isempty(y0_), y0 = y0_(:); GenInpD{7} = y0;end

function EqNam_Callback(hObject, eventdata, handles)

function LiNo_Callback(hObject, eventdata, handles)


function omegAiry_Callback(hObject, eventdata, handles)
global omegAiry
omegAiry_ = CheckScal(hObject);
if ~isempty(omegAiry_), omegAiry = omegAiry_; end

function nBess_Callback(hObject, eventdata, handles)
global nBess
nBess_ = CheckScal(hObject);
if ~isempty(nBess_), nBess = nBess_; end

function PampShrod_Callback(hObject, eventdata, handles)

function PfreqShrod_Callback(hObject, eventdata, handles)

function DelFig_Callback(hObject, eventdata, handles)
global Figs
s = str2num(hObject.String);
n = numel(Figs);
if isempty(s),  errordlg('Wrong input');
elseif s == 0,  for i = 1:n, delete(Figs(i)); end
else
   for j = 1:n
      k = find(s==Figs(j));
      if isempty(k), delete(Figs(j)); 
      else
         s(k) = [];
         if isempty(s) && j<n
            for m = j+1:n, delete(Figs(m));end,end,end,end,end
    

                     % CreateFcn 
% Executes during object creation, after setting all properties
% hObject    handle to hs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.   

function Tasks_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function hs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function Methods_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function a1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white'); end

function y0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function t0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function tn_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function h1_h2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function RelT_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function AbsT_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function nSDE_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function EqNam_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function LiNo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function NEq__CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function omegAiry_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function nBess_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function PampShrod_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function PfreqShrod_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

function IV_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

% --- Executes during object creation, after setting all properties.
function DelFig_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');end

% --- Executes during object creation, after setting all properties.
function cT__CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white'); end
