clc
clear all
close all
% Метод функции потока для цилиндрических катушек
%% ------------------------------------------------------------------------- % Задание параметров

tic

Calc = 1;
OutData = 0;

Web = struct('R', 0.0665, 'L', 0.15, 'nfi', 15, 'nL', 11, 'axial', 0, 'Node', [], 'S', [], 'I', []);
            
CubeROI = struct('nx', 8, 'ny', 8, 'nz', 8, 'lx', 0.05, 'ly', 0.05, 'lz', 0.05, 'center', [0 0 0]);

SphereROI = struct('nR', 5, 'nksi', 14, 'ntetta', 7, 'Rsp', 0.025, 'center', [0 0 0]);

FOV = struct('nx', 15, 'ny', 15, 'nz', 15, 'lx', 0.05, 'ly', 0.05, 'lz', 0.05, 'center', [0 0 0], 'level', 0, 'plane', 'XY', 'FOV', [], 'B_target', [], 'B_received', [], 'B_difference', []);

Grad = struct('type', 'MY', 'value', 0.05, 'boff', 0, 'ROI', [], 'B_target', [], 'B_received', [], 'B_difference', []);

Optimisation = struct('ro', 0.017, 't', 2, 'alpfa', 1500, 'beta', 0.0005);

Contour = struct('zero', 0, 'Nc', 5, 'offb', 0.005, 'offr', 0.005, 'Red_2D', [], 'Blue_2D', [], 'Red_3D', [], 'Blue_3D', []);

Results = struct('Fi', [], 'Ll', [], 'II', [], 'px', [], 'py', []);

ROI_form = 0;   % 0) Куб; 1) Сфера; 2) Произвольно

switch FOV.plane
    case "XY"
        Levels = linspace(FOV.lz/2, -FOV.lz/2, FOV.nz) + FOV.center(3);
    case "XZ"
        Levels = linspace(FOV.ly/2, -FOV.ly/2, FOV.ny) + FOV.center(2);
    case "YZ"
        Levels = linspace(FOV.lx/2, -FOV.lx/2, FOV.nx) + FOV.center(1);
end

disp("Доступные уровни:");
disp(Levels);

FOV.level = Levels(8);
disp("Уровень:");
disp(FOV.level);

%% --------------------------------------------------------------------------- Создание сетки

Web.Node = CreateNode(Web.R, Web.L, Web.nfi, Web.nL, Web.axial);
Web.S = CreateS(Web.nfi*Web.nL, Web.nfi);

FOV.FOV = CreateCubeROI(FOV.nx, FOV.ny, FOV.nz, FOV.lx, FOV.ly, FOV.lz, FOV.center);

switch ROI_form
    case 0
        Grad.ROI = CreateCubeROI(CubeROI.nx, CubeROI.ny, CubeROI.nz, CubeROI.lx, CubeROI.ly, CubeROI.lz, CubeROI.center);
        PointWebROI(Web.Node*1000, Web.S, Grad.ROI*1000, FOV.lx*1000, FOV.ly*1000, FOV.lz*1000, FOV.center*1000)
    case 1
        Grad.ROI = CreateSphereROI(SphereROI.nR, SphereROI.nksi, SphereROI.ntetta, SphereROI.Rsp, SphereROI.center);
        PointWebROI(Web.Node*1000, Web.S, Grad.ROI*1000, FOV.lx*1000, FOV.ly*1000, FOV.lz*1000, FOV.center*1000)
    case 2
        Grad.ROI = [];
        PointWebROI(Web.Node*1000, Web.S, Grad.ROI*1000, FOV.lx*1000, FOV.ly*1000, FOV.lz*1000, FOV.center*1000)  
end

disp(strcat("Узлов: ",num2str(Web.nfi*Web.nL)));
disp(strcat("Элементов: ",num2str(length(Web.S))));

%% --------------------------------------------------------------------------- Создание целевого поля

if ROI_form ~= 2
    switch Grad.type
        case "X"
            Grad.B_target = Grad.value*Grad.ROI(:,1) + Grad.boff;
            FOV.B_target = Grad.value*FOV.FOV(:,1) + Grad.boff;
        case "Y"
            Grad.B_target = Grad.value*Grad.ROI(:,2) + Grad.boff;
            FOV.B_target = Grad.value*FOV.FOV(:,2) + Grad.boff;
        case "Z"
            Grad.B_target = Grad.value*Grad.ROI(:,3) + Grad.boff;
            FOV.B_target = Grad.value*FOV.FOV(:,3) + Grad.boff;
        case "XY"
            Grad.B_target = Grad.value*Grad.ROI(:,1).*Grad.ROI(:,2) + Grad.boff;
            FOV.B_target = Grad.value*FOV.FOV(:,1).*FOV.FOV(:,2) + Grad.boff;
        case "XZ"
            Grad.B_target = Grad.value*Grad.ROI(:,1).*Grad.ROI(:,3) + Grad.boff;
            FOV.B_target = Grad.value*FOV.FOV(:,1).*FOV.FOV(:,3) + Grad.boff;
        case "YZ"
            Grad.B_target = Grad.value*Grad.ROI(:,2).*Grad.ROI(:,3) + Grad.boff;
            FOV.B_target = Grad.value*FOV.FOV(:,2).*FOV.FOV(:,3) + Grad.boff;
        case "X2-Y2"
            Grad.B_target = Grad.value*(Grad.ROI(:,1).*Grad.ROI(:,1) - Grad.ROI(:,2).*Grad.ROI(:,2)) + Grad.boff;
            FOV.B_target = Grad.value*(FOV.FOV(:,1).*FOV.FOV(:,1) - FOV.FOV(:,2).*FOV.FOV(:,2)) + Grad.boff;
        case "X2"
            Grad.B_target = Grad.value*abs(Grad.ROI(:,1)).*Grad.ROI(:,1) + Grad.boff;
            FOV.B_target = Grad.value*abs(FOV.FOV(:,1)).*FOV.FOV(:,1) + Grad.boff;
        case "MY"
            Grad.B_target = Grad.value*(Grad.ROI(:,1)) + Grad.boff;
            FOV.B_target = Grad.value*(FOV.FOV(:,1)) + Grad.boff;
    end
else 
    Grad.B_target = [];
    FOV.B_target = [];
end

switch ROI_form
    case 0 
        PrintFieldCube(FOV.lx, FOV.ly, FOV.lz, FOV.nx, FOV.ny, FOV.nz, FOV.center, FOV.B_target*1000, FOV.level, FOV.plane, CubeROI.lx, CubeROI.ly, CubeROI.lz, CubeROI.center, 'Целевое поле', 'Целевое поле, мТ')
    case 1
        PrintFieldSphere(FOV.lx, FOV.ly, FOV.lz, FOV.nx, FOV.ny, FOV.nz, FOV.center, FOV.B_target*1000, FOV.level, FOV.plane, SphereROI.Rsp, SphereROI.center, 'Целевое поле', 'Целевое поле, мТ')
    case 2 
        disp("Отображение целевого поля не доступно")
end

%% --------------------------------------------------------------------------- Расчет

if ~OutData
    Web.I = ABI(Web.L, Web.R, CubeROI.lx, CubeROI.ly, CubeROI.lz, Web.nfi*Web.nL, Web.Node, Web.S, SphereROI.Rsp, Grad.ROI, Optimisation.ro, Optimisation.t, Optimisation.alpfa, Optimisation.beta, Grad.B_target, ROI_form);
else
    I = [];
end

%% --------------------------------------------------------------------------- Расчет полученного поля

switch ROI_form
    case 0
        [Grad.B_received, Grad.B_difference] = ReceivedField(Web.I, CubeROI.nx*CubeROI.ny*CubeROI.nz, Web.nfi*Web.nL, Web.S, Web.Node, Grad.ROI, Grad.B_target);
    case 1
        [Grad.B_received, Grad.B_difference] = ReceivedField(Web.I, SphereROI.nR*SphereROI.nksi*SphereROI.ntetta, Web.nfi*Web.nL, Web.S, Web.Node, Grad.ROI, Grad.B_target);
    case 2
        [Grad.B_received, Grad.B_difference] = ReceivedField(Web.I, length(Grad.ROI(:,1)), Web.nfi*Web.nL, Web.S, Web.Node, Grad.ROI, Grad.B_target);
end

[FOV.B_received, FOV.B_difference] = ReceivedField(Web.I, FOV.nx*FOV.ny*FOV.nz, Web.nfi*Web.nL, Web.S, Web.Node, FOV.FOV, FOV.B_target);

Krit = max(abs(Grad.B_difference)./max(Grad.B_target, [], 'all'), [], 'all')*100;
FOV.B_difference = 100*(FOV.B_difference)./max(FOV.B_target, [], 'all');

disp(" ")
disp(strcat("Максамальное отклонение: ",num2str(Krit)," %"))

%% --------------------------------------------------------------------------- Результаты

[Results.px, Results.py, Results.Fi, Results.Ll, Results.II, Contour.Red_3D, Contour.Blue_3D, Contour.Red_2D, Contour.Blue_2D] = PrintResults(Web.R, Web.L, Web.nfi, Web.nL, Web.axial, Web.I, Contour.Nc, Contour.offb, Contour.offr, Contour.zero);
            
%% --------------------------------------------------------------------------- Полученное поле

switch ROI_form
    case 0 
        PrintFieldCube(FOV.lx, FOV.ly, FOV.lz, FOV.nx, FOV.ny, FOV.nz, FOV.center, FOV.B_received*1000, FOV.level, FOV.plane, CubeROI.lx, CubeROI.ly, CubeROI.lz, CubeROI.center, 'Полученное поле', 'Полученное поле, мТ')
    case 1
        PrintFieldSphere(FOV.lx, FOV.ly, FOV.lz, FOV.nx, FOV.ny, FOV.nz, FOV.center, FOV.B_received*1000, FOV.level, FOV.plane, SphereROI.Rsp, SphereROI.center, 'Полученное поле', 'Полученное поле, мТ')
    case 2 
        PrintFieldSphere(FOV.lx, FOV.ly, FOV.lz, FOV.nx, FOV.ny, FOV.nz, FOV.center, FOV.B_received*1000, FOV.level, FOV.plane, 0, SphereROI.center, 'Полученное поле', 'Полученное поле, мТ')
end 

switch ROI_form
    case 0 
        PrintFieldCube(FOV.lx, FOV.ly, FOV.lz, FOV.nx, FOV.ny, FOV.nz, FOV.center, FOV.B_difference, FOV.level, FOV.plane, CubeROI.lx, CubeROI.ly, CubeROI.lz, CubeROI.center, 'Разность', 'Разность, %')
    case 1
        PrintFieldSphere(FOV.lx, FOV.ly, FOV.lz, FOV.nx, FOV.ny, FOV.nz, FOV.center, FOV.B_difference, FOV.level, FOV.plane, SphereROI.Rsp, SphereROI.center, 'Разность', 'Разность, %')
    case 2 
        PrintFieldSphere(FOV.lx, FOV.ly, FOV.lz, FOV.nx, FOV.ny, FOV.nz, FOV.center, FOV.B_difference, FOV.level, FOV.plane, 0, SphereROI.center, 'Разность', 'Разность, %')
end 

Time = toc;
disp(strcat("Затраченное время: ", num2str(Time), " с"))


