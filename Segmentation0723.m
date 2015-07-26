clear;
%http://www.mathworks.com/help/images/examples/detecting-a-cell-using-image-segmentation.html
clear;
[I]=dicomread('D:\CODE\yhjmatlab\matlab_yhj\data\BS037\BS037-CTXA.dcm');
PSF = fspecial('gaussian',5,5);
luc1 = deconvlucy(I,PSF,5);
BW2 = edge(luc1,'canny');

%Dilate the Image
se90 = strel('line', 4, 90);
se0 = strel('line', 4, 0);
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
imshow(BWoutline); hold on
%figure;
%imshow(I+im2uint16(BWoutline));

% read two curves
Low_rect = getrect;
Low_rect = round(Low_rect);
Lminx=Low_rect(1);
Lmaxx=Low_rect(1) + Low_rect(3);
Lminy=Low_rect(2);
Lmaxy=Low_rect(2) + Low_rect(4);
[Ly,Lx]=find(BWoutline([Lminy:Lmaxy],[Lminx:Lmaxx]));
Lx=Lx+Lminx-1;
Ly=Ly+Lminy-1;

High_rect = getrect;
High_rect = round(High_rect);
Hminx=High_rect(1);
Hmaxx=High_rect(1) + High_rect(3);
Hminy=High_rect(2);
Hmaxy=High_rect(2) + High_rect(4);
[Hy,Hx]=find(BWoutline([Hminy:Hmaxy],[Hminx:Hmaxx]));
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
        tmpd = (Lx(i_L) - Hx(j_H))^2 + (Ly(i_L) - Hy(j_H))^2;
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
seta = (pi/2) - seta;
M_x=round((Px_H + Px_L)/2);
M_y=round((Py_H + Py_L)/2);

%get perpendicular line
Length_F=70;
End_x1 = M_x - round(cos(seta)*Length_F);
End_y1 = M_y - round(sin(seta)*Length_F);
End_x2 = M_x + round(cos(seta)*Length_F);
End_y2 = M_y + round(sin(seta)*Length_F);

plot([End_x1, End_x2], [End_y1, End_y2], 'r*-');

%Find midle line at the end
End_rect = getrect;
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
    O=nnz(BWoutline(i_E,:));
    if (O==2)
        [TP_y,TP_x]=find(BWoutline(i_E,:));
        Mid_x(E_n)=(TP_x(1)+TP_x(2))/2;
        Mid_y(E_n)=i_E;
        E_n=E_n+1;
    end
end
E_n=E_n-1;
Mid_x(1,1:E_n-1)=round(min(Mid_x));
plot([Mid_x(1),Mid_x(E_n)],[Mid_y(1)-70,Mid_y(E_n)],  'g-');

%find crosspoint
PA=[Mid_y(E_n)-Mid_y(1), Mid_x(1)-Mid_x(E_n); End_y2 - End_y1, End_x1 - End_x2];
PB=[((Mid_y(E_n)-Mid_y(1))*Mid_x(1) - (Mid_x(E_n)-Mid_x(1))*Mid_y(1)); ((End_y2 - End_y1)*End_x1 -(End_x2 - End_x1)*End_y1)];
[Cp]=PA\PB;
Cp=round(abs(Cp));
plot(Cp(1),Cp(2),'*');

V1=(Mid_x(1) - Cp(1) + (Cp(2) - (Mid_y(1)-70))*i)/abs(Mid_x(1) - Cp(1) + (Cp(2) - (Mid_y(1)-70))*i);
V2=(End_x2 - Cp(1) + ( Cp(2) - End_y2 )*i)/abs(End_x2 - Cp(1) + (Cp(2) - End_y2 )*i);
Vplus1=50*(V1+V2);
Intetrochanter_x1=real(Vplus1)+Cp(1);
Intetrochanter_y1=Cp(2)-imag(Vplus1);
line([Cp(1),Intetrochanter_x1],[Cp(2),Intetrochanter_y1]);

V3=(Mid_x(E_n) - Cp(1) + (Cp(2) - (Mid_y(E_n)))*i)/abs(Mid_x(E_n) - Cp(1) + (Cp(2) - (Mid_y(E_n)))*i);
V4=(End_x1 - Cp(1) + ( Cp(2) - End_y1 )*i)/abs(End_x1 - Cp(1) + (Cp(2) - End_y1 )*i);
Vplus2=50*(V3+V4);
Intetrochanter_x2=real(Vplus2)+Cp(1);
Intetrochanter_y2=Cp(2)-imag(Vplus2);
line([Cp(1),Intetrochanter_x2],[Cp(2),Intetrochanter_y2]);


line([Cp(1)-20,Cp(1)+20],[Cp(2)+1.5*sqrt(d_min),Cp(2)+1.5*sqrt(d_min)]);













