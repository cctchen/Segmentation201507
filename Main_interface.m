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

% Last Modified by GUIDE v2.5 03-Aug-2015 16:44:49

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
Dicom_folder = uigetdir();% after select a fold, the path will be assigned to Dicom_folder
[tmp,list]=system(['dir ' Dicom_folder '\*CTXA.dcm /S/B']);%get the path
Dicom_list = strsplit(list);%split the string
[temp,filenumber] = size(Dicom_list);
for ii=1:filenumber
    [pathstr,dname,ext] = fileparts(Dicom_list{ii});%seperate 'Dicom_list{ii}'into path, name and format.
    Dicom_name{ii}=dname;
end
set(handles.listbox1,'string',Dicom_name);%set the name of string to 'listbox1',which can display the name list in the box1


handles.Dicom_list=Dicom_list;%seems like put 'Dicom_list'into 'handls'.handles like a capacity ,can include lots things,and it serves for all function.
handles.Dicom_name=Dicom_name;
guidata(hObject,handles);%save those things into this two capacity.

axes(handles.axes1);%display the pic in ases1
index_selected=1;
dcm_pic=dicomread(Dicom_list{1});
imshow(dcm_pic);
p_Dilate=get(handles.slider1,'value');%when creat 'slider1',it has been given a default value which is 4(dilate parameter).
handles.p_Dilate=p_Dilate;
handles.index_selected=index_selected;
guidata(hObject,handles);
handles.Image_name=Dicom_name(index_selected);
guidata(hObject,handles);%save those things into this two capacity.
get_outline(hObject, handles, p_Dilate);

% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
index_selected = get(handles.listbox1,'Value');%  get the number in the order of the cases when click the case name in the box1
axes(handles.axes1);%display the pic in ases1
dcm_pic=dicomread(handles.Dicom_list{index_selected});
imshow(dcm_pic);
p_Dilate=get(handles.slider1,'value');%when creat 'slider1',it has been given a default value which is 4(dilate parameter).
handles.p_Dilate=p_Dilate;
handles.index_selected=index_selected;
guidata(hObject,handles);
handles.Image_name=handles.Dicom_name(index_selected);
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
p_Dilate=get(handles.slider1,'value');%when clike slider1, will get the number of dilate parameter
handles.p_Dilate=p_Dilate;
guidata(hObject,handles);
set(handles.text5,'string',p_Dilate);%display the dilate parameter
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
[tmp,size_outl]=size(handles.BWoutline);
handles.size_outl=size_outl;
handles.I=I;
guidata(hObject,handles);
imshow(BWoutline); hold on


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
FUC=1.073;
INTERCEPT=-1024.3;
SLOP=1.2882;
Pixel_scale=0.97656e-3;
mu=0.3;

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

%calculate EI_FN
V1=(Px_L - M_x + (Py_L - M_y)*i)/abs(Px_L - M_x + (Py_L - M_y)*i);
V2=(Px_H - M_x + (Py_H - M_y)*i)/abs(Px_H - M_x + (Py_H - M_y)*i);

ex_L=10;
V_hight=5/0.97656;

V_AL_1=ex_L*V1+V_hight*V1*(-i);
V_AL_2=ex_L*V1+V_hight*V1*i;

V_AH_1=ex_L*V2+V_hight*V2*(-i);
V_AH_2=ex_L*V2+V_hight*V2*i;

POLY_x=[real(V_AL_1)+Px_L, real(V_AL_2)+Px_L, real(V_AH_1)+Px_H, real(V_AH_2)+Px_H];
POLY_y=[imag(V_AL_1)+Py_L, imag(V_AL_2)+Py_L, imag(V_AH_1)+Py_H, imag(V_AH_2)+Py_H];
POLY_x_show=[real(V_AL_1)+Px_L, real(V_AL_2)+Px_L, real(V_AH_1)+Px_H, real(V_AH_2)+Px_H, real(V_AL_1)+Px_L];
POLY_y_show=[imag(V_AL_1)+Py_L, imag(V_AL_2)+Py_L, imag(V_AH_1)+Py_H, imag(V_AH_2)+Py_H, imag(V_AL_1)+Py_L];
plot(POLY_x_show, POLY_y_show, 'r-');
xa = 1 : size(handles.I,2);
ya = 1 : size(handles.I,1);
[xq,yq] = meshgrid(xa,ya);
[in]=inpolygon(xq,yq,POLY_x,POLY_y);
I_rec = and(in, handles.I);
[Rec_y,Rec_x]=find(I_rec);
Rec_size=size(Rec_x);
Rec_pixel_sum=0;
for(ii=1:Rec_size(1))
    Rec_pixel_sum=double(Rec_pixel_sum) + double(handles.I(Rec_y(ii),Rec_x(ii)));
end

Ave_Pixel= Rec_pixel_sum/(Rec_size(1)*Pixel_scale^2);

set(handles.text19,'string',int2str(Rec_pixel_sum));
set(handles.text20,'string',int2str(Ave_Pixel));

r=1;
nfn_x1=M_x;
nfn_y1=M_y;
ei_px=nfn_x1;
ei_py=nfn_y1;
rou=(double(FUC)*(double(handles.I(ei_py,ei_px))-double(INTERCEPT)))/double(SLOP);%mg/cm3
E=double(29.8)*(double((1e-3))*double(rou))^1.56*1e9;
EA_FN=double(double(E)*double(Pixel_scale^2));
EI_FN = 0;
while((nfn_x1 < handles.size_outl) && (nfn_x1 > 0) && (handles.BWoutline(nfn_y1,nfn_x1) ~= 1))
    Vplus1=r*V1;
    nfn_x1=round(M_x + real(Vplus1));
    nfn_y1=round(M_y + imag(Vplus1)); % along the vector to find the edge
    if((nfn_x1 ~=ei_px) || (nfn_y1 ~= ei_py ))
        ei_px=nfn_x1;
        ei_py=nfn_y1;
        if (handles.I(ei_py,ei_px) > 0)
            rou=(double(FUC)*(double(handles.I(ei_py,ei_px))-double(INTERCEPT)))/double(SLOP);%mg/cm3
            E=double(29.8)*(double((1e-3))*double(rou))^1.56*1e9;
            EA_FN=EA_FN + double(double(E)*double(Pixel_scale^2));
            EI_FN = EI_FN + E*((r*Pixel_scale)^2)*double(Pixel_scale^2);   
        end
    end
    if ((handles.BWoutline(nfn_y1-1,nfn_x1) == 1) && (handles.BWoutline(nfn_y1,nfn_x1+1) == 1))
        break;
    end
    r=r+0.5;
end

r=1;
nfn_x2=M_x;
nfn_y2=M_y;
ei_px=nfn_x2;
ei_py=nfn_y2;
while((nfn_x2 < handles.size_outl) && (nfn_x2 > 0) && handles.BWoutline(nfn_y2,nfn_x2) ~= 1)%??????????????????????????
    Vplus2=r*V2;
    nfn_x2=round(M_x + real(Vplus2));
    nfn_y2=round(M_y + imag(Vplus2));
    if((nfn_x2 ~=ei_px) || (nfn_y2 ~= ei_py ))% if they equal, that means they are in the same pixel value, so do not need to calculate again.
        ei_px=nfn_x2;
        ei_py=nfn_y2;
        if(handles.I(ei_py,ei_px) > 0)
            rou=(double(FUC)*(double(handles.I(ei_py,ei_px))-double(INTERCEPT)))/double(SLOP);%mg/cm3
            E=(double(29.8)*(double((1e-3))*double(rou))^1.56)*1e9;
            EA_FN=EA_FN + double(double(E)*double(Pixel_scale^2));
            EI_FN = EI_FN + E*((r*Pixel_scale)^2)*double(Pixel_scale^2);
        end
    end
    if ((handles.BWoutline(nfn_y2,nfn_x2-1) == 1) && (handles.BWoutline(nfn_y2+1,nfn_x2) == 1))
        break;
    end
    r=r+0.5;
end
GA_FN=EA_FN/(2*(1+mu));
set(handles.text15,'string',num2str(EI_FN));

%get perpendicular line

V_neck=(Px_H - M_x + (Py_H - M_y)*i)/abs((Px_H - M_x + (Py_H - M_y)*i));
V_axis=-V_neck*i;

perpen_x1 = round((Px_L + Px_H)/2);
perpen_y1 = round((Py_L + Py_H)/2);
perpen_x2 = round((Px_L + Px_H)/2);
perpen_y2 = round((Py_L + Py_H)/2);

r=1;
while((perpen_x1 > 1) && (perpen_x1 < handles.size_outl) && (handles.BWoutline(perpen_y1,perpen_x1) ~= 1))
    Vd=r*V_axis;
    perpen_x1=round((Px_L + Px_H)/2 + real(Vd));
    perpen_y1=round((Py_L + Py_H)/2 + imag(Vd)); % along the vector to find the edge
    if ((handles.BWoutline(perpen_y1+1,perpen_x1) == 1) && (handles.BWoutline(perpen_y1,perpen_x1+1) == 1))
        break;
    end
    r=r+0.5;
end

r=-1;
while((perpen_x2 > 1) && (perpen_x2 < handles.size_outl) && (handles.BWoutline(perpen_y2,perpen_x2) ~= 1))
    Vd=r*V_axis;
    perpen_x2=round((Px_L + Px_H)/2+real(Vd));
    perpen_y2=round((Py_L + Py_H)/2+imag(Vd)); % along the vector to find the edge
    if ((handles.BWoutline(perpen_y2+1,perpen_x2) == 1) && (handles.BWoutline(perpen_y2,perpen_x2-1) == 1))
        break;
    end
    r=r-0.5;
end

%Length_F=70;
%perpen_x1 = M_x - round(cos(seta)*Length_F);%sin(seta)=(delta y)/70,so delta y=70*sin(seta),similarly, delta x=70*cos(seta)
%perpen_y1 = M_y - round(sin(seta)*Length_F);
%perpen_x2 = M_x + round(cos(seta)*Length_F);
%perpen_y2 = M_y + round(sin(seta)*Length_F);

plot([perpen_x1, perpen_x2], [perpen_y1, perpen_y2], 'r-');

%Find axis of the femur shaft
set(handles.text1,'string','Select Shaft area');
shaft_rect = getrect;
set(handles.text1,'string','');
shaft_rect = round(shaft_rect);
shaft_minx=shaft_rect(1);
shaft_maxx=shaft_rect(1) + shaft_rect(3);
shaft_miny=shaft_rect(2);
shaft_maxy=shaft_rect(2) + shaft_rect(4);
%[Endy,Endx]=find(BWoutline([shaft_miny:shaft_maxy],[shaft_minx:shaft_maxx]));
%Endx=Endx+shaft_minx-1;
%Endy=Endy+shaft_miny-1;
E_n=1;
for i_E=shaft_miny:shaft_maxy
    O=nnz(handles.BWoutline(i_E,:));%find the how many none zero pixel in each line
    if (O==2)
        [TP_y,TP_x]=find(handles.BWoutline(i_E,:)); %if the none zero number=2, means they are on the outline edge, and get their coordinate.
        Mid_x(E_n)=(TP_x(1)+TP_x(2))/2;
        Mid_y(E_n)=i_E;
        E_n=E_n+1;
    end
end
E_n=E_n-1;
Mid_x(1:E_n)=round(min(Mid_x));
%plot([Mid_x(1),Mid_x(E_n)],[Mid_y(1)-70,Mid_y(E_n)],  'g-');

%find crosspoint between shaft axis and neck axis
%using two points deciding one line, creat two line functions, then calculate their public solution which is the crosspoint(using matris method) 
PA=[Mid_y(E_n)-Mid_y(1), Mid_x(1)-Mid_x(E_n); perpen_y2 - perpen_y1, perpen_x1 - perpen_x2];
PB=[((Mid_y(E_n)-Mid_y(1))*Mid_x(1) - (Mid_x(E_n)-Mid_x(1))*Mid_y(1)); ((perpen_y2 - perpen_y1)*perpen_x1 -(perpen_x2 - perpen_x1)*perpen_y1)];
[Cp]=PA\PB;% left multiply
Cp=round(abs(Cp));
plot(Cp(1),Cp(2),'*');

%V1=(Mid_x(1) - Cp(1) + (Cp(2) - (Mid_y(1)-70))*i)/abs(Mid_x(1) - Cp(1) + (Cp(2) - (Mid_y(1)-70))*i);
%V2=(End_x2 - Cp(1) + ( Cp(2) - End_y2 )*i)/abs(End_x2 - Cp(1) + (Cp(2) - End_y2 )*i);
%Vplus1=50*(V1+V2);
%Intetrochanter_x1=real(Vplus1)+Cp(1);
%Intetrochanter_y1=Cp(2)-imag(Vplus1);
%line([Cp(1),Intetrochanter_x1],[Cp(2),Intetrochanter_y1]);

V1=(Mid_x(E_n) - Cp(1) + ((Mid_y(E_n))-Cp(2))*i)/abs(Mid_x(E_n) - Cp(1) + ((Mid_y(E_n))-Cp(2))*i);%get the unit vector of the two line
V2=(perpen_x1 - Cp(1) + (perpen_y1-Cp(2))*i)/abs(perpen_x1 - Cp(1) + (perpen_y1-Cp(2))*i);

r=1;
Intetrochanter_x1=Cp(1);
Intetrochanter_y1=Cp(2);
ei_px=Intetrochanter_x1;
ei_py=Intetrochanter_y1;
rou=(double(FUC)*(double(handles.I(ei_py,ei_px))-double(INTERCEPT)))/double(SLOP);%mg/cm3
E=double(29.8)*(double((1e-3))*double(rou))^1.56*1e9;
EA_IT=double(double(E)*double(Pixel_scale^2));
EI_IT = 0;
while((Intetrochanter_x1 < handles.size_outl) && (Intetrochanter_x1 > 0) && (handles.BWoutline(Intetrochanter_y1,Intetrochanter_x1) ~= 1))
    Vplus1=r*(V1+V2);
    Intetrochanter_x1=round(Cp(1)+real(Vplus1));
    Intetrochanter_y1=round(Cp(2)+imag(Vplus1)); % along the vector to find the edge
    if((Intetrochanter_x1 ~=ei_px) || (Intetrochanter_y1 ~= ei_py ))
        ei_px=Intetrochanter_x1;
        ei_py=Intetrochanter_y1;
        if(handles.I(ei_py,ei_px) > 0)
            rou=(double(FUC)*(double(handles.I(ei_py,ei_px))-double(INTERCEPT)))/double(SLOP);%mg/cm3
            E=double(29.8)*(double((1e-3))*double(rou))^1.56*1e9;
            EA_IT=EA_IT + double(double(E)*double(Pixel_scale^2));
            EI_IT = EI_IT + E*((r*Pixel_scale)^2)*double(Pixel_scale^2); 
        end
    end
    if ((handles.BWoutline(Intetrochanter_y1-1,Intetrochanter_x1) == 1) && (handles.BWoutline(Intetrochanter_y1,Intetrochanter_x1+1) == 1))
        break;
    end
    r=r+0.5;
end
%line([Cp(1),Intetrochanter_x1],[Cp(2),Intetrochanter_y1]);

r=-1;
Intetrochanter_x2=Cp(1);
Intetrochanter_y2=Cp(2);
ei_px=Intetrochanter_x2;
ei_py=Intetrochanter_y2;
while((Intetrochanter_x2 < handles.size_outl) && (Intetrochanter_x2 > 0) && handles.BWoutline(Intetrochanter_y2,Intetrochanter_x2) ~= 1)
    Vplus2=r*(V1+V2);
    Intetrochanter_x2=round(Cp(1)+real(Vplus2));
    Intetrochanter_y2=round(Cp(2)+imag(Vplus2));
    if((Intetrochanter_x2 ~=ei_px) || (Intetrochanter_y2 ~= ei_py ))
        ei_px=Intetrochanter_x2;
        ei_py=Intetrochanter_y2;
        if(handles.I(ei_py,ei_px) > 0)
            rou=(double(FUC)*(double(handles.I(ei_py,ei_px))-double(INTERCEPT)))/double(SLOP);%mg/cm3
            E=double(29.8)*(double((1e-3))*double(rou))^1.56*1e9;
            EA_IT=EA_IT + double(double(E)*double(Pixel_scale^2));
            EI_IT = EI_IT + E*((r*Pixel_scale)^2)*double(Pixel_scale^2);
        end
    end
    if ((handles.BWoutline(Intetrochanter_y2,Intetrochanter_x2-1) == 1) && (handles.BWoutline(Intetrochanter_y2+1,Intetrochanter_x2) == 1))
        break;
    end
    r=r-0.5;
end
GA_IT=EA_IT/(2*(1+mu));
set(handles.text16,'string',num2str(EI_IT));%???????????
plot([Intetrochanter_x1,Intetrochanter_x2],[Intetrochanter_y1,Intetrochanter_y2], '*-');

[shafty,shaftx]=find(handles.BWoutline(Cp(2)+round(1.5*sqrt(d_min)),:));%find the two points of the shaft line

shafty=shafty+Cp(2)+round(1.5*sqrt(d_min))-1;

%calculate EI_FS
V1=1;
V2=-1;

r=1;
fs_x1=Cp(1);
fs_y1=shafty(1);
ei_px=fs_x1;
ei_py=fs_y1;
rou=(double(FUC)*(double(handles.I(ei_py,ei_px))-double(INTERCEPT)))/double(SLOP);%mg/cm3
E=double(29.8)*(double((1e-3))*double(rou))^1.56*1e9;
EA_FS=double(double(E)*double(Pixel_scale^2));
EI_FS = 0;
while((fs_x1 < handles.size_outl) && (fs_x1 > 0) && (handles.BWoutline(fs_y1,fs_x1) ~= 1))
    Vplus1=r*V1;
    fs_x1=round(Cp(1) + real(Vplus1));%??????????????
    fs_y1=round(shafty(1)); % along the vector to find the edge
    if((fs_x1 ~=ei_px) || (fs_y1 ~= ei_py ))
        ei_px=fs_x1;
        ei_py=fs_y1;
        if(handles.I(ei_py,ei_px) > 0)
            rou=(double(FUC)*(double(handles.I(ei_py,ei_px))-double(INTERCEPT)))/double(SLOP);%mg/cm3
            E=double(29.8)*(double((1e-3))*double(rou))^1.56*1e9;
            EA_FS=EA_FS + double(double(E)*double(Pixel_scale^2));
            EI_FS = EI_FS + E*((r*Pixel_scale)^2)*double(Pixel_scale^2);
        end
    end
    r=r+0.5;
end

r=1;
fs_x2=Cp(1);
fs_y2=shafty(1);
ei_px=fs_x2;
ei_py=fs_y2;
while((fs_x2 < handles.size_outl) && (fs_x2 > 0) && handles.BWoutline(fs_y2,fs_x2) ~= 1)
    Vplus2=r*V2;
    fs_x2=round(Cp(1) + real(Vplus2));
    fs_y2=round(shafty(1));
    if((fs_x2 ~=ei_px) || (fs_y2 ~= ei_py )) 
        ei_px=fs_x2;
        ei_py=fs_y2;
        if(handles.I(ei_py,ei_px) > 0)
            rou=(double(FUC)*(double(handles.I(ei_py,ei_px))-double(INTERCEPT)))/double(SLOP);%mg/cm3
            E=double(29.8)*(double((1e-3))*double(rou))^1.56*1e9;
            EA_FS=EA_FS + double(double(E)*double(Pixel_scale^2));
            EI_FS = EI_FS + E*((r*Pixel_scale)^2)*double(Pixel_scale^2);
        end
    end
    r=r+0.5;
end
GA_FS=EA_FS/(2*(1+mu));
set(handles.text17,'string',num2str(EI_FS));%??????

%line([Cp(1)-20,Cp(1)+20],[Cp(2)+1.5*sqrt(d_min),Cp(2)+1.5*sqrt(d_min)]);
axes(handles.axes2);
imshow(handles.I+im2uint16(handles.BWoutline));hold on
plot([Px_L, Px_H], [Py_L, Py_H], 'r*-'); hold on
plot([perpen_x1, perpen_x2], [perpen_y1, perpen_y2], 'r-');
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

%find head coordinator
axes(handles.axes2);
set(handles.text1,'string','Select head area');
head_rec=getrect;
head_rec = round(head_rec);
head_minx=head_rec(1);
head_maxx=head_rec(1) + head_rec(3);
head_miny=head_rec(2);
head_maxy=head_rec(2) + head_rec(4);
to_head_x = round((Px_L + Px_H)/2);
to_head_y = round((Py_L + Py_H)/2);
r=1;
min_p=inf;
head_x=0;
head_y=0;
while((to_head_x > 1) && (to_head_x < handles.size_outl) && (handles.BWoutline(to_head_y,to_head_x) ~= 1))
    Vd=r*V_axis;
    to_head_x=round((Px_L + Px_H)/2 + real(Vd));
    to_head_y=round((Py_L + Py_H)/2 + imag(Vd)); % along the vector to find the edge
    if ((handles.BWoutline(to_head_y+1,to_head_x) == 1) && (handles.BWoutline(to_head_y,to_head_x+1) == 1))
        break;
    end
    if((to_head_x >= head_minx) && (to_head_x <= head_maxx) && (to_head_y >= head_miny) && (to_head_y <= head_maxy) && (handles.I(to_head_y,to_head_x)<min_p))
        min_p=handles.I(to_head_y,to_head_x);
        head_x=to_head_x;
        head_y=to_head_y;
    end
    r=r+0.5;
end
set(handles.text1,'string','');
plot(head_x,head_y,'*'); 

set(handles.edit1,'string',strcat('NFN_L_x=',int2str(Px_L),', NFN_L_y=',int2str(Py_L),'; NFN_H_x=',int2str(Px_H),', NFN_H_y=',int2str(Py_H)));
set(handles.edit2,'string',strcat('IT_L_x=',int2str(Intetrochanter_x1),', IT_L_y=',int2str(Intetrochanter_y1),'; IT_H_x=',int2str(Intetrochanter_x2),', IT_H_y=',int2str(Intetrochanter_y2)));
set(handles.edit3,'string',strcat('FS_L_x=',int2str(shaftx(1)),', FS_L_y=',int2str(shafty(1)),'; FS_R_x=',int2str(shaftx(2)),', FS_R_y=',int2str(shafty(2))));
set(handles.edit4,'string',strcat('HEAD_x=',int2str(head_x),', HEAD_y=',int2str(head_y),'; MIT_x=',int2str(Cp(1)),', MIT_y=',int2str(Cp(2))));

handles.Image_name=char(handles.Dicom_name(handles.index_selected));
handles.EA_FN=num2str(EA_FN);
handles.EA_IT=num2str(EA_IT);
handles.EA_FS=num2str(EA_FS);
handles.EI_FN=num2str(EI_FN);
handles.EI_IT=num2str(EI_IT);
handles.EI_FS=num2str(EI_FS);
handles.GA_FN=num2str(GA_FN);
handles.GA_IT=num2str(GA_IT);
handles.GA_FS=num2str(GA_FS);
handles.FN1_x=num2str(Px_L);
handles.FN1_y=num2str(Py_L);
handles.FN2_x=num2str(Px_H);
handles.FN2_y=num2str(Py_H);
handles.IT1_x=num2str(Intetrochanter_x1);
handles.IT1_y=num2str(Intetrochanter_y1);
handles.IT2_x=num2str(Intetrochanter_x2);
handles.IT2_y=num2str(Intetrochanter_y2);
handles.FS1_x=num2str(shaftx(1));
handles.FS1_y=num2str(shafty(1));
handles.FS2_x=num2str(shaftx(2));
handles.FS2_y=num2str(shafty(2));
handles.HEAD_x=num2str(head_x);
handles.HEAD_y=num2str(head_y);
handles.CP_x=num2str(Cp(1));
handles.CP_y=num2str(Cp(2));
guidata(hObject,handles);

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



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)


% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
conn = database('seg_results','admin','123456');
handles.conn=conn;
guidata(hObject,handles);

% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
close(handles.conn);
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)


colnames = {'Image_Name'  'EI_FN' 'EI_IT' 'EI_FS' 'EA_FN' 'EA_IT' 'EA_FS' 'GA_FN' 'GA_IT' 'GA_FS' 'FN1_x' 'FN1_y' 'FN2_x' 'FN2_y' 'IT1_x' 'IT1_y' 'IT2_x' 'IT2_y' 'FS1_x' 'FS1_y' 'FS2_x' 'FS2_y' 'HEAD_x' 'HEAD_y' 'CP_x' 'CP_y'};
edata={handles.Image_name, handles.EA_FN, handles.EA_IT, handles.EA_FS, handles.EI_FN, handles.EI_IT, handles.EI_FS, handles.GA_FN, handles.GA_IT, handles.GA_FS, handles.FN1_x, handles.FN1_y, ...
    handles.FN2_x, handles.FN2_y, handles.IT1_x, handles.IT1_y, handles.IT2_x, handles.IT2_y, handles.FS1_x, handles.FS1_y, handles.FS2_x, handles.FS2_y, handles.HEAD_x, handles.HEAD_y, handles.CP_x, handles.CP_y};
curs=exec(handles.conn, strcat('select Image_name from results where Image_name= ''', handles.Image_name, ''''));
curs=fetch(curs);
isem=strcmp(curs.Data, 'No Data');
if (1 == isem)
    insert(handles.conn, 'results', colnames, edata);
else
    whereclause = strcat('where Image_Name = ''', handles.Image_name, '''');
    update(handles.conn,'results',colnames,edata,whereclause);    
end
sqlquery = 'commit';
exec(handles.conn,sqlquery);
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
