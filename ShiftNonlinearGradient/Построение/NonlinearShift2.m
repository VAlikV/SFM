clear all
close all
clc

tic

%--------------------------------------------------------------------------

R = 0.07; % Радиус
R1 = 0.08;
R2 = 0.09;
nn = 2; % Дополнительные точки

%uiimport('1.txt');
data1 = readtable('YZ_1.txt');
data = table2array(data1)./1000;
YZ1 = zeros(length(data(:,1)),3);
YZ1(:,1) = data(:,2);
YZ1(:,2) = data(:,1);
YZ1 = transZ(YZ1,R);
YZ1 = flex(YZ1,R);
YZ1 = InterpCoil(YZ1, nn);

data1 = readtable('YZ_2.txt');
data = table2array(data1)./1000;
YZ2 = zeros(length(data(:,1)),3);
YZ2(:,1) = data(:,2);
YZ2(:,2) = data(:,1);
YZ2 = transZ(YZ2,R);
YZ2 = flex(YZ2,R);
YZ2 = InterpCoil(YZ2, nn);

data1 = readtable('YZ_3.txt');
data = table2array(data1)./1000;
YZ3 = zeros(length(data(:,1)),3);
YZ3(:,1) = data(:,2);
YZ3(:,2) = data(:,1);
YZ3 = transZ(YZ3,R);
YZ3 = flex(YZ3,R);
YZ3 = InterpCoil(YZ3, nn);

data1 = readtable('YZ_4.txt');
data = table2array(data1)./1000;
YZ4 = zeros(length(data(:,1)),3);
YZ4(:,1) = data(:,2);
YZ4(:,2) = data(:,1);
YZ4 = transZ(YZ4,R);
YZ4 = flex(YZ4,R);
YZ4 = InterpCoil(YZ4, nn);

data1 = readtable('YZ_5.txt');
data = table2array(data1)./1000;
YZ5 = zeros(length(data(:,1)),3);
YZ5(:,1) = data(:,2);
YZ5(:,2) = data(:,1);
YZ5 = transZ(YZ5,R);
YZ5 = flex(YZ5,R);
YZ5 = InterpCoil(YZ5, nn);

data1 = readtable('YZ_6.txt');
data = table2array(data1)./1000;
YZ6 = zeros(length(data(:,1)),3);
YZ6(:,1) = data(:,2);
YZ6(:,2) = data(:,1);
YZ6 = transZ(YZ6,R);
YZ6 = flex(YZ6,R);
YZ6 = InterpCoil(YZ6, nn);

% YZ1 = rotX(YZ1, -0.5*pi/180);
% YZ2 = rotX(YZ2, -0.5*pi/180);
% YZ3 = rotX(YZ3, -0.5*pi/180);
% YZ4 = rotX(YZ4, -0.5*pi/180);
% YZ5 = rotX(YZ5, -0.5*pi/180);
% YZ6 = rotX(YZ6, -0.5*pi/180);

% ------------------------------------------------------- Z

data1 = readtable('Z_1.txt');
data = table2array(data1)./1000;
Z1 = zeros(length(data(:,1)),3);
Z1(:,1) = data(:,2);
Z1(:,2) = data(:,1);
Z1 = transZ(Z1,R1);
Z1 = flex(Z1,R1);
Z1 = InterpCoil(Z1, nn);

data1 = readtable('Z_2.txt');
data = table2array(data1)./1000;
Z2 = zeros(length(data(:,1)),3);
Z2(:,1) = data(:,2);
Z2(:,2) = data(:,1);
Z2 = transZ(Z2,R1);
Z2 = flex(Z2,R1);
Z2 = InterpCoil(Z2, nn);

data1 = readtable('Z_3.txt');
data = table2array(data1)./1000;
Z3 = zeros(length(data(:,1)),3);
Z3(:,1) = data(:,2);
Z3(:,2) = data(:,1);
Z3 = transZ(Z3,R1);
Z3 = flex(Z3,R1);
Z3 = InterpCoil(Z3, nn);

data1 = readtable('Z_4.txt');
data = table2array(data1)./1000;
Z4 = zeros(length(data(:,1)),3);
Z4(:,1) = data(:,2);
Z4(:,2) = data(:,1);
Z4 = transZ(Z4,R1);
Z4 = flex(Z4,R1);
Z4 = InterpCoil(Z4, nn);

Z1 = rotX(Z1, -0.2*pi/180);
Z2 = rotX(Z2, -0.2*pi/180);
Z3 = rotX(Z3, -0.2*pi/180);
Z4 = rotX(Z4, -0.2*pi/180);

% ------------------------------------------------------- Y

data1 = readtable('Y_1.txt');
data = table2array(data1)./1000;
Y1 = zeros(length(data(:,1)),3);
Y1(:,1) = data(:,2);
Y1(:,2) = data(:,1);
Y1 = transZ(Y1,R2);
Y1 = flex(Y1,R2);
Y1 = InterpCoil(Y1, nn);

data1 = readtable('Y_2.txt');
data = table2array(data1)./1000;
Y2 = zeros(length(data(:,1)),3);
Y2(:,1) = data(:,2);
Y2(:,2) = data(:,1);
Y2 = transZ(Y2,R2);
Y2 = flex(Y2,R2);
Y2 = InterpCoil(Y2, nn);

data1 = readtable('Y_3.txt');
data = table2array(data1)./1000;
Y3 = zeros(length(data(:,1)),3);
Y3(:,1) = data(:,2);
Y3(:,2) = data(:,1);
Y3 = transZ(Y3,R2);
Y3 = flex(Y3,R2);
Y3 = InterpCoil(Y3, nn);

data1 = readtable('Y_4.txt');
data = table2array(data1)./1000;
Y4 = zeros(length(data(:,1)),3);
Y4(:,1) = data(:,2);
Y4(:,2) = data(:,1);
Y4 = transZ(Y4,R2);
Y4 = flex(Y4,R2);
Y4 = InterpCoil(Y4, nn);

Y1 = rotX(Y1, -0.2*pi/180);
Y2 = rotX(Y2, -0.2*pi/180);
Y3 = rotX(Y3, -0.2*pi/180);
Y4 = rotX(Y4, -0.2*pi/180);

% ------------------------------------------------------- 

Iyz = 10; % Величина тока в катушках
Iz = 5;
Iy = 5;

lx = 0.05;
ly = 0.05;
lz = 0.05;

nx = 25;
ny = 25;
nz = 25;

centerROI = [0 0 0];

plane = "YZ";

switch plane
    case "XY"
        Levels = linspace(lz/2, -lz/2, nz) + centerROI(3);
    case "XZ"
        Levels = linspace(ly/2, -ly/2, ny) + centerROI(2);
    case "YZ"
        Levels = linspace(lx/2, -lx/2, nx) + centerROI(1);
end

disp("Доступные уровни:");
disp(Levels);

level = 0.0;

ROI = CreateCubeROI(nx, ny, nz, lx, ly, lz, centerROI);
 
%--------------------------------------------------------------------------

% Создание катушек

figure(3)
plot3(YZ1(:,1),YZ1(:,2),YZ1(:,3),'b-')  % 'b-'
%plot3(Z1(:,1),Z1(:,2),Z1(:,3),'b-')    % 'b-'
%plot3(Y1(:,1),Y1(:,2),Y1(:,3),'b-')`   % 'b-'
hold on

plot3(YZ2(:,1),YZ2(:,2),YZ2(:,3),'b-')  % 'r-'
plot3(YZ3(:,1),YZ3(:,2),YZ3(:,3),'b-')  % 'b-'
plot3(YZ4(:,1),YZ4(:,2),YZ4(:,3),'b-')  % 'r-'
plot3(YZ5(:,1),YZ5(:,2),YZ5(:,3),'b-')  % 'b-'
plot3(YZ6(:,1),YZ6(:,2),YZ6(:,3),'b-')  % 'r-'

plot3(Z1(:,1),Z1(:,2),Z1(:,3),'r-')  % 'b-'
plot3(Z2(:,1),Z2(:,2),Z2(:,3),'r-')  % 'b-'
plot3(Z3(:,1),Z3(:,2),Z3(:,3),'r-')  % 'r-'
plot3(Z4(:,1),Z4(:,2),Z4(:,3),'r-')  % 'r-'

plot3(Y1(:,1),Y1(:,2),Y1(:,3),'g-')  % 'b-'
plot3(Y2(:,1),Y2(:,2),Y2(:,3),'g-')  % 'r-'
plot3(Y3(:,1),Y3(:,2),Y3(:,3),'g-')  % 'r-'
plot3(Y4(:,1),Y4(:,2),Y4(:,3),'g-')  % 'b-'
xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]')
axis equal
grid on

%--------------------------------------------------------------------------

Byz = zeros(length(ROI(:,1)),3);
Bz = zeros(length(ROI(:,1)),3);
By = zeros(length(ROI(:,1)),3);

Byz = Byz + BSL(YZ1, Iyz, ROI);
Byz = Byz + BSL(YZ2, Iyz, ROI);
Byz = Byz + BSL(YZ3, Iyz, ROI);
Byz = Byz + BSL(YZ4, Iyz, ROI);
Byz = Byz + BSL(YZ5, Iyz, ROI);
Byz = Byz + BSL(YZ6, Iyz, ROI);

Bz = Bz + BSL(Z1, Iz, ROI);
Bz = Bz + BSL(Z2, Iz, ROI);
Bz = Bz + BSL(Z3, Iz, ROI);
Bz = Bz + BSL(Z4, Iz, ROI);

By = By + BSL(Y1, Iy, ROI);
By = By + BSL(Y2, Iy, ROI);
By = By + BSL(Y3, Iy, ROI);
By = By + BSL(Y4, Iy, ROI);

B = Byz + By + Bz;
%--------------------------------------------------------------------------

PrintFieldCube(lx, ly, lz, nx, ny, nz, centerROI, B(:,3)*1000, level, plane, lx, ly, lz, centerROI, "Полученное поле", "Полученное поле, мТ");

toc
