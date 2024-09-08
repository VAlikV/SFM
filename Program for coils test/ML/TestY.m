clc
close all 
clear all

% ------------------------------------------------------------------- Контур

R = 0.05*1000; % Радиус

data1 = readtable('Y1.txt');
data = table2array(data1)./1;
Y1 = zeros(length(data(:,1)),3);
Y1(:,1) = data(:,2);
Y1(:,2) = data(:,1);
Y1 = transZ(Y1,R);
Y1 = flex(Y1,R);

data1 = readtable('Y2.txt');
data = table2array(data1)./1;
Y2 = zeros(length(data(:,1)),3);
Y2(:,1) = data(:,2);
Y2(:,2) = data(:,1);
Y2 = transZ(Y2,R);
Y2 = flex(Y2,R);

data1 = readtable('Y3.txt');
data = table2array(data1)./1;
Y3 = zeros(length(data(:,1)),3);
Y3(:,1) = data(:,2);
Y3(:,2) = data(:,1);
Y3 = transZ(Y3,R);
Y3 = flex(Y3,R);

data1 = readtable('Y4.txt');
data = table2array(data1)./1;
Y4 = zeros(length(data(:,1)),3);
Y4(:,1) = data(:,2);
Y4(:,2) = data(:,1);
Y4 = transZ(Y4,R);
Y4 = flex(Y4,R);

% w = 1.5;
% figure('Name','Контур1','NumberTitle','off'); 
% plot3(Y1(:,1),Y1(:,2),Y1(:,3),'r-','LineWidth',w);
% hold on
% plot3(Y2(:,1),Y2(:,2),Y2(:,3),'b-','LineWidth',w);
% plot3(Y3(:,1),Y3(:,2),Y3(:,3),'r-','LineWidth',w);
% plot3(Y4(:,1),Y4(:,2),Y4(:,3),'b-','LineWidth',w);
% xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]')
% axis equal
% grid on 

nnn = 5;

Y1(:,1) = smooth(Y1(:,1),nnn);
Y1(:,2) = smooth(Y1(:,2),nnn);
Y1(:,3) = smooth(Y1(:,3),nnn);

Y2(:,1) = smooth(Y2(:,1),nnn);
Y2(:,2) = smooth(Y2(:,2),nnn);
Y2(:,3) = smooth(Y2(:,3),nnn);

Y3(:,1) = smooth(Y3(:,1),nnn);
Y3(:,2) = smooth(Y3(:,2),nnn);
Y3(:,3) = smooth(Y3(:,3),nnn);

Y4(:,1) = smooth(Y4(:,1),nnn);
Y4(:,2) = smooth(Y4(:,2),nnn);
Y4(:,3) = smooth(Y4(:,3),nnn);

% ------------------------------------------------------------------- ROI

nx = 15;    % Кол-во точек вдоль Х
ny = 15;    % Кол-во точек вдоль Y
nz = 15;    % Кол-во точек вдоль Z
K = nx*ny*nz;  % Полное число точек ROI

lx = 30/1000;  % Длина области вдоль Х
ly = 30/1000;  % Длина области вдоль Y
lz = 30/1000;  % Длина области вдоль Z

CenterROI = [0 0 0]; % Положение центра ROI

ROI = CreateCubeROI(nx, ny, nz, lx, ly, lz, CenterROI);

level = 0.0;  % Положение плоскости
plane = "XY";  % Выбор плоскости XY, XZ, YZ

if plane == "XY"
    disp(linspace(-lz/2, lz/2, nz)+CenterROI(3));
elseif plane == "XZ"
    disp(linspace(-ly/2, ly/2, ny)+CenterROI(2));
elseif plane == "YZ"
    disp(linspace(-lx/2, lx/2, nx)+CenterROI(1));
end

I = -10;

% ------------------------------------------------------------------- Расчет

B = zeros(length(ROI(:,1)),3);
B = B + BSL(Y1, I, ROI);
B = B + BSL(Y2, -I, ROI);
B = B + BSL(Y3, I, ROI);
B = B + BSL(Y4, -I, ROI);

PrintFieldCube(lx, ly, lz, nx, ny, nz, CenterROI, B(:,3)*10000, level, plane, lx, ly, lz, CenterROI, 'Полученное поле', 'Полученное поле, Гс', "contourf") % Отрисовка поля

% ------------------------------------------------------------------- Отрисовка

w = 1.5;
figure('Name','Контур','NumberTitle','off'); 
plot3(Y1(:,1),Y1(:,2),Y1(:,3),'r-','LineWidth',w);
hold on
plot3(Y2(:,1),Y2(:,2),Y2(:,3),'b-','LineWidth',w);
plot3(Y3(:,1),Y3(:,2),Y3(:,3),'r-','LineWidth',w);
plot3(Y4(:,1),Y4(:,2),Y4(:,3),'b-','LineWidth',w);
xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]')
axis equal
grid on 


