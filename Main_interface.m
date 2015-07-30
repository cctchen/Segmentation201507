function varargout = Main_interface(varargin)
% MAIN_INTERFACE MATLAB code for Main_interface.fig
%      MAIN_INTERFACE, by itself, creates a new MAIN_INTERFACE or raises the existing
%      singleton*.
%
%      H = MAIN_INTERFACE returns the handle to a new MAIN_INTERFACE or the handle to
%      the existing singleton*.
%
%      MAIN_INTERFACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN_INTERFACE.M with the given input arguments.
%
%      MAIN_INTERFACE('Property','Value',...) creates a new MAIN_INTERFACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Main_interface_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Main_interface_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Main_interface

% Last Modified by GUIDE v2.5 30-Jul-2015 23:51:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Main_interface_OpeningFcn, ...
                   'gui_OutputFcn',  @Main_interface_OutputFcn, ...
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


% --- Executes just before Main_interface is made visible.
function Main_interface_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Main_interface (see VARARGIN)

% Choose default command line output for Main_interface
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Main_interface wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Main_interface_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
Dicom_folder = uigetdir();
[tmp,list]=system(['dir ' Dicom_folder '\*CTXA.dcm /S/B']);
Dicom_list = strsplit(list);
[temp,filenumber] = size(Dicom_list);
for ii=1:filenumber
    [pathstr,dname,ext] = fileparts(Dicom_list{ii});
    Dicom_name{ii}=dname;
end
set(handles.listbox1,'string',Dicom_name);
handles.Dicom_list=Dicom_list;
guidata(hObject,handles);
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
index_selected = get(handles.listbox1,'Value');
axes(handles.axes1);
dcm_pic=dicomread(handles.Dicom_list{index_selected});
imshow(dcm_pic);
p_Dilate=get(handles.slider1,'value');
handles.p_Dilate=p_Dilate;
handles.index_selected=index_selected;
guidata(hObject,handles);
get_outline(hObject, handles, p_Dilate);

% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in slider1.
function slider1_Callback(hObject, eventdata, handles)
p_Dilate=get(handles.slider1,'value');
handles.p_Dilate=p_Dilate;
guidata(hObject,handles);
set(handles.text5,'string',p_Dilate);
get_outline(hObject, handles, p_Dilate);
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of slider1


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function get_outline(hObject, handles, p_Dilate)
[I]=dicomread(handles.Dicom_list{handles.index_selected});
PSF = fspecial('gaussian',5,5);
luc1 = deconvlucy(I,PSF,5);
BW2 = edge(luc1,'canny');

%Dilate the Image
se90 = strel('line', p_Dilate, 90);
se0 = strel('line', p_Dilate, 0);
BWsdil = imdilate(BW2, [se90 se0]);

%Fill Interior Gaps
BWdfill = imfill(BWsdil, 'holes');

%Smoothen the Object
seD = strel('diamond',1);
BWfinal = imerode(BWdfill,seD);
BWfinal = imerode(BWfinal,seD);

%Get outline
BWoutline = bwperim(BWfinal);
Segout = I;
Segout(BWoutline) = 255;
axes(handles.axes2);
handles.BWoutline=BWoutline;
handles.I=I;
guidata(hObject,handles);
imshow(BWoutline); hold on


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
axes(handles.axes2);
% read two curves
set(handles.text1,'string','Select low neck arc.');
Low_rect = getrect;
Low_rect = round(Low_rect); %transfer into integer
Lminx=Low_rect(1);
Lmaxx=Low_rect(1) + Low_rect(3);
Lminy=Low_rect(2);
Lmaxy=Low_rect(2) + Low_rect(4);
[Ly,Lx]=find(handles.BWoutline([Lminy:Lmaxy],[Lminx:Lmaxx]));%find the location where the value>0
Lx=Lx+Lminx-1; %transfer into the local coordinate
Ly=Ly+Lminy-1;


set(handles.text1,'string','Select high neck arc.');
High_rect = getrect;
High_rect = round(High_rect);
Hminx=High_rect(1);
Hmaxx=High_rect(1) + High_rect(3);
Hminy=High_rect(2);
Hmaxy=High_rect(2) + High_rect(4);
[Hy,Hx]=find(handles.BWoutline([Hminy:Hmaxy],[Hminx:Hmaxx]));
Hx=Hx+Hminx-1;
Hy=Hy+Hminy-1;


%Find two closest points

Px_L=0;
Py_L=0;
Px_H=0;
Py_H=0;
size_L=size(Lx);
size_H=size(Hx);
d_min=Inf;

for i_L=1:size_L
    for j_H=1:size_H
        tmpd = (Lx(i_L) - Hx(j_H))^2 + (Ly(i_L) - Hy(j_H))^2;%distance between two points
        if(tmpd<d_min)
            d_min=tmpd;
            Px_L=Lx(i_L);
            Px_H=Hx(j_H);
            Py_L=Ly(i_L);
            Py_H=Hy(j_H);
        end       
    end
end

plot([Px_L, Px_H], [Py_L, Py_H], 'r*-'); hold on

seta = atan(abs(Py_H-Py_L)/abs(Px_H-Px_L));
seta = (pi/2) - seta; %calculat the angle of the perpendicular line
M_x=round((Px_H + Px_L)/2); % the center of the nearest femur neck
M_y=round((Py_H + Py_L)/2);

%get perpendicular line
Length_F=70;
End_x1 = M_x - round(cos(seta)*Length_F);%sin(seta)=(delta y)/70,so delta y=70*sin(seta),similarly, delta x=70*cos(seta)
End_y1 = M_y - round(sin(seta)*Length_F);
End_x2 = M_x + round(cos(seta)*Length_F);
End_y2 = M_y + round(sin(seta)*Length_F);

plot([End_x1, End_x2], [End_y1, End_y2], 'r-');


%Find axis of the femur shaft
set(handles.text1,'string','Select Shaft area');
End_rect = getrect;
set(handles.text1,'string','');
End_rect = round(End_rect);
End_minx=End_rect(1);
End_maxx=End_rect(1) + End_rect(3);
End_miny=End_rect(2);
End_maxy=End_rect(2) + End_rect(4);
%[Endy,Endx]=find(BWoutline([End_miny:End_maxy],[End_minx:End_maxx]));
%Endx=Endx+End_minx-1;
%Endy=Endy+End_miny-1;
E_n=1;
for i_E=End_miny:End_maxy
    O=nnz(handles.BWoutline(i_E,:));
    if (O==2)
        [TP_y,TP_x]=find(handles.BWoutline(i_E,:));
        Mid_x(E_n)=(TP_x(1)+TP_x(2))/2;
        Mid_y(E_n)=i_E;
        E_n=E_n+1;
    end
end
E_n=E_n-1;
Mid_x(1:E_n)=round(min(Mid_x));
plot([Mid_x(1),Mid_x(E_n)],[Mid_y(1)-70,Mid_y(E_n)],  'g-');

%find crosspoint between shaft axis and neck axis
%using two points deciding one line, creat two line functions, then calculate their public solution which is the crosspoint(using matris method) 
PA=[Mid_y(E_n)-Mid_y(1), Mid_x(1)-Mid_x(E_n); End_y2 - End_y1, End_x1 - End_x2];
PB=[((Mid_y(E_n)-Mid_y(1))*Mid_x(1) - (Mid_x(E_n)-Mid_x(1))*Mid_y(1)); ((End_y2 - End_y1)*End_x1 -(End_x2 - End_x1)*End_y1)];
[Cp]=PA\PB;% left multiply
Cp=round(abs(Cp));
plot(Cp(1),Cp(2),'*');

%V1=(Mid_x(1) - Cp(1) + (Cp(2) - (Mid_y(1)-70))*i)/abs(Mid_x(1) - Cp(1) + (Cp(2) - (Mid_y(1)-70))*i);
%V2=(End_x2 - Cp(1) + ( Cp(2) - End_y2 )*i)/abs(End_x2 - Cp(1) + (Cp(2) - End_y2 )*i);
%Vplus1=50*(V1+V2);
%Intetrochanter_x1=real(Vplus1)+Cp(1);
%Intetrochanter_y1=Cp(2)-imag(Vplus1);
%line([Cp(1),Intetrochanter_x1],[Cp(2),Intetrochanter_y1]);

V1=(Mid_x(E_n) - Cp(1) + ((Mid_y(E_n))-Cp(2))*i)/abs(Mid_x(E_n) - Cp(1) + ((Mid_y(E_n))-Cp(2))*i);
V2=(End_x1 - Cp(1) + (End_y1-Cp(2))*i)/abs(End_x1 - Cp(1) + (End_y1-Cp(2))*i);

[tmp,size_outl]=size(handles.BWoutline);
Intetrochanter_x1=Cp(1);
Intetrochanter_y1=Cp(2);
r=1;
while((Intetrochanter_x1 < size_outl) && (Intetrochanter_x1 > 0) && (handles.BWoutline(Intetrochanter_y1,Intetrochanter_x1) ~= 1))
    Vplus1=r*(V1+V2);
    Intetrochanter_x1=round(Cp(1)+real(Vplus1));
    Intetrochanter_y1=round(Cp(2)+imag(Vplus1));
    r=r+0.5;
end
%line([Cp(1),Intetrochanter_x1],[Cp(2),Intetrochanter_y1]);
r=-1;
Intetrochanter_x2=Cp(1);
Intetrochanter_y2=Cp(2);
while((Intetrochanter_x2 < size_outl) && (Intetrochanter_x2 > 0) && handles.BWoutline(Intetrochanter_y2,Intetrochanter_x2) ~= 1)
    Vplus2=r*(V1+V2);
    Intetrochanter_x2=round(Cp(1)+real(Vplus2));
    Intetrochanter_y2=round(Cp(2)+imag(Vplus2));
    r=r-0.5;
end

plot([Intetrochanter_x1,Intetrochanter_x2],[Intetrochanter_y1,Intetrochanter_y2], '*-');

[shafty,shaftx]=find(handles.BWoutline(Cp(2)+round(1.5*sqrt(d_min)),:));

shafty=shafty+Cp(2)+round(1.5*sqrt(d_min));

%line([Cp(1)-20,Cp(1)+20],[Cp(2)+1.5*sqrt(d_min),Cp(2)+1.5*sqrt(d_min)]);
plot(shaftx,shafty,'*');

axes(handles.axes2);
imshow(handles.I+im2uint16(handles.BWoutline));hold on
plot([Px_L, Px_H], [Py_L, Py_H], 'r*-'); hold on
plot([End_x1, End_x2], [End_y1, End_y2], 'r-');
plot([Cp(1),Mid_x(E_n)],[Cp(2)-20,Mid_y(E_n)],  'g-');
plot(Cp(1),Cp(2),'*');
plot([Intetrochanter_x1,Intetrochanter_x2],[Intetrochanter_y1,Intetrochanter_y2], '*-');
%line([Cp(1)-20,Cp(1)+20],[Cp(2)+1.5*sqrt(d_min),Cp(2)+1.5*sqrt(d_min)]);
plot(shaftx,shafty,'r*-');
axes(handles.axes3);
improfile(handles.I, [Px_L, Px_H], [Py_L, Py_H]);
axes(handles.axes4);
improfile(handles.I, [Intetrochanter_x1,Intetrochanter_x2],[Intetrochanter_y1,Intetrochanter_y2]);
axes(handles.axes5);
improfile(handles.I, shaftx,shafty);
set(handles.edit1,'string',strcat('NFN_L_x=',int2str(Px_L),', NFN_L_y=',int2str(Py_L),'; NFN_H_x=',int2str(Px_H),', NFN_H_y=',int2str(Py_H)));
set(handles.edit2,'string',strcat('IT_L_x=',int2str(Intetrochanter_x1),', IT_L_y=',int2str(Intetrochanter_y1),'; IT_H_x=',int2str(Intetrochanter_x2),', IT_H_y=',int2str(Intetrochanter_y2)));
set(handles.edit3,'string',strcat('FS_L_x=',int2str(shaftx(1)),', FS_L_y=',int2str(shafty(1)),'; FS_R_x=',int2str(shaftx(2)),', FS_R_y=',int2str(shafty(2))));

% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
