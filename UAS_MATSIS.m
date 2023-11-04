function varargout = UAS_MATSIS(varargin)
% UAS_MATSIS M-file for UAS_MATSIS.fig
%      UAS_MATSIS, by itself, creates a new UAS_MATSIS or raises the existing
%      singleton*.
%
%      H = UAS_MATSIS returns the handle to a new UAS_MATSIS or the handle to
%      the existing singleton*.
%
%      UAS_MATSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UAS_MATSIS.M with the given input arguments.
%
%      UAS_MATSIS('Property','Value',...) creates a new UAS_MATSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before UAS_MATSIS_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UAS_MATSIS_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UAS_MATSIS

% Last Modified by GUIDE v2.5 26-May-2014 14:32:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @UAS_MATSIS_OpeningFcn, ...
                   'gui_OutputFcn',  @UAS_MATSIS_OutputFcn, ...
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


% --- Executes just before UAS_MATSIS is made visible.
function UAS_MATSIS_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to UAS_MATSIS (see VARARGIN)

% Choose default command line output for UAS_MATSIS
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes UAS_MATSIS wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = UAS_MATSIS_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function n0_Callback(hObject, eventdata, handles)
% hObject    handle to n0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n0 as text
%        str2double(get(hObject,'String')) returns contents of n0 as a double


% --- Executes during object creation, after setting all properties.
function n0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a0_Callback(hObject, eventdata, handles)
% hObject    handle to a0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a0 as text
%        str2double(get(hObject,'String')) returns contents of a0 as a double


% --- Executes during object creation, after setting all properties.
function a0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function b0_Callback(hObject, eventdata, handles)
% hObject    handle to b0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of b0 as text
%        str2double(get(hObject,'String')) returns contents of b0 as a double


% --- Executes during object creation, after setting all properties.
function b0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to b0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in process1.
function process1_Callback(hObject, eventdata, handles)
% hObject    handle to process1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
k1=0.5;k2=0.7;k3=0.5;k4=0.9;alpha=1.1;beta=2;
A=[(-1)*(k2+k3)*k4/alpha*beta*k2 (-1)*(k2+k3)/alpha 0;...
    (k2+k3)*k4/beta*k2 0 (-1)*(k2+k3)*k4/alpha*beta*k1;...
    -k4 beta*k2 0
    ];
B=[1;0;0];
C=[0 0 1];
%% NONLINEAR STEP
h=0.2;
t=0:h:30;
% n(1)=2;a(1)=1;b(1)=0.01;
handles.x=str2num(get(handles.n0,'String'));
handles.y=str2num(get(handles.a0,'String'));
handles.z=str2num(get(handles.b0,'String'));
n(1)=handles.x;a(1)=handles.y;b(1)=handles.z;
for i=1:length(t)-1
    u(i)=1;
    n1=h*(u(i)-k1*n(i)*a(i));
    a1=h*(k1*alpha*n(i)*a(i)-k2*a(i)*b(i)-k3*a(i));
    b1=h*(k2*beta*a(i)*b(i)-k4*n(i)*b(i));
    
    n2=h*(u(i)-k1*(n(i)+0.5*n1)*(a(i)+0.5*a1));
    a2=h*(k1*alpha*(n(i)+0.5*n1)*(a(i)+0.5*a1)-k2*(a(i)+0.5*a1)*(b(i)+0.5*b1)-k3*(a(i)+0.5*a1));
    b2=h*(k2*beta*(a(i)+0.5*a1)*(b(i)+0.5*b1)-k4*(n(i)+0.5*n1)*(b(i)+0.5*b1));
    
    n3=h*(u(i)-k1*(n(i)+0.5*n2)*(a(i)+0.5*a2));
    a3=h*(k1*alpha*(n(i)+0.5*n2)*(a(i)+0.5*a2)-k2*(a(i)+0.5*a2)*(b(i)+0.5*b2)-k3*(a(i)+0.5*a2));
    b3=h*(k2*beta*(a(i)+0.5*a2)*(b(i)+0.5*b2)-k4*(n(i)+0.5*n2)*(b(i)+0.5*b2));
    
    n4=h*(u(i)-k1*(n(i)+n3)*(a(i)+a3));
    a4=h*(k1*alpha*(n(i)+n3)*(a(i)+a3)-k2*(a(i)+a3)*(b(i)+b3)-k3*(a(i)+a3));
    b4=h*(k2*beta*(a(i)+a3)*(b(i)+b3)-k4*(n(i)+n3)*(b(i)+b3));
    
    n(i+1)=n(i)+(1/6)*(n1+2*n2+2*n3+n4);
    a(i+1)=a(i)+(1/6)*(a1+2*a2+2*a3+a4);
    b(i+1)=b(i)+(1/6)*(b1+2*b2+2*b3+b4);
end;
axes(handles.axes1);
plot(t,n,'b');
hold on;
plot(t,a,'g');
hold on;
plot(t,b,'r');
xlabel('Time (hours)');
ylabel('Concentration');
legend('Nutrient','Algae','Rotifer');
axes(handles.axes2);
plot(a,b,'m');
title('Rotifer Versus Algae');
xlabel('Algae Concentration');
ylabel('Rotifer Concentration');
%% LINEAR STEP
h_linear=0.2;
t_linear=0:h_linear:50;
A1=-(k2+k3)*k4/(alpha*beta*k2);
A2=-(k2+k3)/alpha;
A3=(k2+k3)*k4/(beta*k2);
A4=-(k2+k3)*k4/(alpha*beta*k1);
A5=-k4;
A6=beta*k2;
n_star=(k2+k3)/alpha*k1;
a_star=(k2+k3)*k4/alpha*beta*k1*k2;
b_star=1;
u_star=((k2+k3)^2)*k4/(alpha^2)*beta*(k2^2);
n_linear(1)=handles.x;a_linear(1)=handles.y;b_linear(1)=handles.z;
xe1(1)=n_linear(1)-n_star;
xe2(1)=a_linear(1)-a_star;
xe3(1)=b_linear(1)-b_star;
for i=1:length(t_linear)-1
    u(i)=1;
    ue(i)=u(i)-u_star;
    xe1_value=h_linear*(A1*xe1(i)+A2*xe2(i)+ue(i));
    xe2_value=h_linear*(A3*xe1(i)+A4*xe3(i));
    xe3_value=h_linear*(A5*xe1(i)+A6*xe2(i));
    
    xe1_value2=h_linear*(A1*(xe1(i)+0.5*xe1_value)+A2*(xe2(i)+0.5*xe2_value)+ue(i));
    xe2_value2=h_linear*(A3*(xe1(i)+0.5*xe1_value)+A4*(xe3(i)+0.5*xe3_value));
    xe3_value2=h_linear*(A5*(xe1(i)+0.5*xe1_value)+A6*(xe2(i)+0.5*xe2_value));
    
    xe1_value3=h_linear*(A1*(xe1(i)+0.5*xe1_value2)+A2*(xe2(i)+0.5*xe2_value2)+ue(i));
    xe2_value3=h_linear*(A3*(xe1(i)+0.5*xe1_value2)+A4*(xe3(i)+0.5*xe3_value2));
    xe3_value3=h_linear*(A5*(xe1(i)+0.5*xe1_value2)+A6*(xe2(i)+0.5*xe2_value2));
    
    xe1_value4=h_linear*(A1*(xe1(i)+xe1_value3)+A2*(xe2(i)+xe2_value3)+ue(i));
    xe2_value4=h_linear*(A3*(xe1(i)+xe1_value3)+A4*(xe3(i)+xe3_value3));
    xe3_value4=h_linear*(A5*(xe1(i)+xe1_value3)+A6*(xe2(i)+xe2_value3));
    
    xe1(i+1)=xe1(i)+(1/6)*(xe1_value+2*xe1_value2+2*xe1_value3+xe1_value4);
    xe2(i+1)=xe2(i)+(1/6)*(xe2_value+2*xe2_value2+2*xe2_value3+xe2_value4);
    xe3(i+1)=xe3(i)+(1/6)*(xe3_value+2*xe3_value2+2*xe3_value3+xe3_value4);
    
    n_linear(i+1)=xe1(i)+n_star;
    a_linear(i+1)=xe2(i)+a_star;
    b_linear(i+1)=xe3(i)+b_star;
end;
axes(handles.axes3);
plot(t_linear,n_linear,'b');
hold on;
plot(t_linear,a_linear,'g');
hold on;
plot(t_linear,b_linear,'r');
xlabel('Time (hours)');
ylabel('Concentration');
legend('Nutrient','Algae','Rotifer');
DATA1=[n_linear' a_linear' b_linear'];
set(handles.data1,'data',DATA1);
guidata(hObject, handles);



function f1_Callback(hObject, eventdata, handles)
% hObject    handle to f1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f1 as text
%        str2double(get(hObject,'String')) returns contents of f1 as a double


% --- Executes during object creation, after setting all properties.
function f1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f2_Callback(hObject, eventdata, handles)
% hObject    handle to f2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f2 as text
%        str2double(get(hObject,'String')) returns contents of f2 as a double


% --- Executes during object creation, after setting all properties.
function f2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f3_Callback(hObject, eventdata, handles)
% hObject    handle to f3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f3 as text
%        str2double(get(hObject,'String')) returns contents of f3 as a double


% --- Executes during object creation, after setting all properties.
function f3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in feedback.
function feedback_Callback(hObject, eventdata, handles)
% hObject    handle to feedback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
k1=0.5;k2=0.7;k3=0.5;k4=0.9;alpha=1.1;beta=2;
n_star=(k2+k3)/alpha*k1;
a_star=(k2+k3)*k4/alpha*beta*k1*k2;
b_star=1;
u_star=((k2+k3)^2)*k4/(alpha^2)*beta*(k2^2);
h_linear=0.2;
t_linear=0:h_linear:50;
A1=-(k2+k3)*k4/(alpha*beta*k2);
A2=-(k2+k3)/alpha;
A3=(k2+k3)*k4/(beta*k2);
A4=-(k2+k3)*k4/(alpha*beta*k1);
A5=-k4;
A6=beta*k2;
n_fb(1)=handles.x;a_fb(1)=handles.y;b_fb(1)=handles.z;
xe1_fb(1)=n_fb(1)-n_star;
xe2_fb(1)=a_fb(1)-a_star;
xe3_fb(1)=b_fb(1)-b_star;
% a_feedback=-5;
% b_feedback=-1;
% c_feedback=1;
a_feedback=str2num(get(handles.f1,'String'));
b_feedback=str2num(get(handles.f2,'String'));
c_feedback=str2num(get(handles.f3,'String'));
for i=1:length(t_linear)-1
   xe1_aksen=h_linear*((A1+a_feedback)*xe1_fb(i)+(A2+b_feedback)*xe2_fb(i)+c_feedback*xe3_fb(i));
   xe2_aksen=h_linear*(A3*xe1_fb(i)+A4*xe3_fb(i));
   xe3_aksen=h_linear*(A5*xe1_fb(i)+A6*xe2_fb(i));
   
   xe1_aksen2=h_linear*((A1+a_feedback)*(xe1_fb(i)+0.5*xe1_aksen)+(A2+b_feedback)*(xe2_fb(i)+0.5*xe2_aksen)+c_feedback*(xe3_fb(i)+0.5*xe3_aksen));
   xe2_aksen2=h_linear*(A3*(xe1_fb(i)+0.5*xe1_aksen)+A4*(xe3_fb(i)+0.5*xe3_aksen));
   xe3_aksen2=h_linear*(A5*(xe1_fb(i)+0.5*xe1_aksen)+A6*(xe2_fb(i)+0.5*xe2_aksen));   
   
   xe1_aksen3=h_linear*((A1+a_feedback)*(xe1_fb(i)+0.5*xe1_aksen2)+(A2+b_feedback)*(xe2_fb(i)+0.5*xe2_aksen2)+c_feedback*(xe3_fb(i)+0.5*xe3_aksen2));
   xe2_aksen3=h_linear*(A3*(xe1_fb(i)+0.5*xe1_aksen2)+A4*(xe3_fb(i)+0.5*xe3_aksen2));
   xe3_aksen3=h_linear*(A5*(xe1_fb(i)+0.5*xe1_aksen2)+A6*(xe2_fb(i)+0.5*xe2_aksen2));
   
   xe1_aksen4=h_linear*((A1+a_feedback)*(xe1_fb(i)+xe1_aksen3)+(A2+b_feedback)*(xe2_fb(i)+xe2_aksen3)+c_feedback*(xe3_fb(i)+xe3_aksen3));
   xe2_aksen4=h_linear*(A3*(xe1_fb(i)+xe1_aksen3)+A4*(xe3_fb(i)+xe3_aksen3));
   xe3_aksen4=h_linear*(A5*(xe1_fb(i)+xe1_aksen3)+A6*(xe2_fb(i)+xe2_aksen3));
   
   xe1_fb(i+1)=xe1_fb(i)+(1/6)*(xe1_aksen+2*xe1_aksen2+2*xe1_aksen3+xe1_aksen4);
   xe2_fb(i+1)=xe2_fb(i)+(1/6)*(xe2_aksen+2*xe2_aksen2+2*xe2_aksen3+xe2_aksen4);
   xe3_fb(i+1)=xe3_fb(i)+(1/6)*(xe3_aksen+2*xe3_aksen2+2*xe3_aksen3+xe3_aksen4);
   
   n_fb(i+1)=xe1_fb(i+1)+n_star;
   a_fb(i+1)=xe2_fb(i+1)+a_star;
   b_fb(i+1)=xe3_fb(i+1)+b_star;
end;
axes(handles.axes11);
plot(t_linear,n_fb,'b');
hold on;
plot(t_linear,a_fb,'g');
hold on;
plot(t_linear,b_fb,'r');
xlabel('Time (hours)');
ylabel('Concentration');
DATA2=[n_fb' a_fb' b_fb'];
set(handles.data2,'data',DATA2);
legend('Nutrient','Algae','Rotifer');


% --- Executes on button press in reset1.
function reset1_Callback(hObject, eventdata, handles)
% hObject    handle to reset1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
axes(handles.axes1);
cla reset;
axes(handles.axes2);
cla reset;
axes(handles.axes3);
cla reset;


% --- Executes on button press in reset2.
function reset2_Callback(hObject, eventdata, handles)
% hObject    handle to reset2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
axes(handles.axes11);
cla reset;


% --- Executes on button press in clear1.
function clear1_Callback(hObject, eventdata, handles)
% hObject    handle to clear1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.data1,'data',[]);


% --- Executes on button press in clear2.
function clear2_Callback(hObject, eventdata, handles)
% hObject    handle to clear2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.data2,'data',[]);
