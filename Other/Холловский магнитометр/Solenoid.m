clc
close all 
clear all

% ------------------------------------------------------------------- Контур

R = 150/1000; % Радиус
L = 1500/1000;
n_points = 100;
n_vit = 200;

[Coil, S] = Spiral(n_points, n_vit, R, L);

% ------------------------------------------------------------------- ROI

nx = 15;    % Кол-во точек вдоль Х
ny = 15;    % Кол-во точек вдоль Y
nz = 15;    % Кол-во точек вдоль Z
K = nx*ny*nz;  % Полное число точек ROI

lx = 50/1000;  % Длина области вдоль Х
ly = 50/1000;  % Длина области вдоль Y
lz = 50/1000;  % Длина области вдоль Z

CenterROI = [0 0 1500/2000]; % Положение центра ROI

ROI = CreateCubeROI(nx, ny, nz, lx, ly, lz, CenterROI);

level = 0.75;  % Положение плоскости
plane = "XY";  % Выбор плоскости XY, XZ, YZ

if plane == "XY"
    disp(linspace(-lz/2, lz/2, nz)+CenterROI(3));
elseif plane == "XZ"
    disp(linspace(-ly/2, ly/2, ny)+CenterROI(2));
elseif plane == "YZ"
    disp(linspace(-lx/2, lx/2, nx)+CenterROI(1));
end

I = 50;

% ------------------------------------------------------------------- Расчет

B = BSL(Coil, I, ROI);

PrintFieldCube(lx, ly, lz, nx, ny, nz, CenterROI, B(:,3)*10000, level, plane, lx, ly, lz, CenterROI, 'Полученное поле', 'Полученное поле, Гс', "surf") % Отрисовка поля

% ------------------------------------------------------------------- Отрисовка

figure('Name','Контур','NumberTitle','off'); 
plot3(Coil(:,1),Coil(:,2),Coil(:,3),'r-');
hold on
plot3(ROI(:,1),ROI(:,2),ROI(:,3),'k.');
xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]')
axis equal
grid on 


