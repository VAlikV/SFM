clc
clear all
close all
% Метод функции потока для цилиндрических катушек
%% ------------------------------------------------------------------------- % Задание параметров

tic

% Размеры
R = 0.0665; % Радиус катушки
nfi = 25;  % Количество точек окружности
%nfi = 14;

L = 0.15; % Половина длины катушки (от 0 в одну сторону)
nL = 21;  % Количество точек по длинне
%nL = 6;

axial = 0;

% ROI
nx = 7;    % Кол-во точек вдоль Х
ny = 7;    % Кол-во точек вдоль Y
nz = 7;    % Кол-во точек вдоль Z

lx = 0.066;  % Длина области вдоль Х
ly = 0.066;  % Длина области вдоль Y
lz = 0.066;  % Длина области вдоль Z

CenterROI = [0 0 0]; % Положение центра ROI

level = 0.0;  % Положение плоскости
plane = "XY";  % Выбор плоскости 0)XY, 1)XZ, 2)YZ

gamma_x = 0;
gamma_y = 0;
gamma_z = 1;

%% ------------------------------------------------------------------------- % Создание массива узлов

N = nfi*nL; % Количество узлов
K = nx*ny*nz;  % Полное число точек ROI

Node = CreateNode(R, L, nfi, nL, axial); % Создание массива узлов

%% ------------------------------------------------------------------------- % Создание массива ROI

ROI = CreateCubeROI(nx, ny, nz, lx, ly, lz, CenterROI);

% [file,path3] = uigetfile('*.txt');
% path3 = strcat(path3,file);
% fileID = fopen(path3,'r');
% 
% if fileID ~= -1
%     data = table2array(readtable(path3));
%     if length(data(1,:)) == 3
%         ROI = data;
%     end          
% end

%% ------------------------------------------------------------------------- % Создание магнитных моментов

Magnets = zeros(N,3,2);
Magnets(:,:,1) = Node;
Magnets(:,2,2) = Node(:,2)/R;
Magnets(:,3,2) = Node(:,3)/R;

%% ------------------------------------------------------------------------- % Задание целевого поля

%Bdes = ROI(:,3).*ROI(:,3) - (ROI(:,1).*ROI(:,1)+ROI(:,2).*ROI(:,2))/2; % Вычисление целевого поля

Bx = zeros(length(ROI(:,3)),1);
By = zeros(length(ROI(:,3)),1);
Bz = ROI(:,1);

[X, Y, Z, BBx, BBy, BBz] = Transform3DMagnets(lx, ly, lz, nx, ny, nz, CenterROI, Bx, By, Bz);

[X1, Y1, Z1] = Display(X, Y, Z, BBx, nx, ny, nz, level, plane, "Целевое поле Bx");
[X1, Y1, Z1] = Display(X, Y, Z, BBy, nx, ny, nz, level, plane, "Целевое поле By");
[X1, Y1, Z1] = Display(X, Y, Z, BBz, nx, ny, nz, level, plane, "Целевое поле Bz");

%data1 = readtable('B_target.txt');
%Bdes = 0.153346313636364 - table2array(data1)./1000;

%M = mean(Bdes,"all");
%TargetField(lx, ly, lz, nx, ny, nz, CenterROI, zeros(K,1), zeros(K,1), Bdes, level, plane);

%% ------------------------------------------------------------------------- % Расчет поля

I = OptM(N, Magnets, ROI, Bx, By, Bz, gamma_x, gamma_y, gamma_z);
Magnets(:,1,2) = Magnets(:,1,2).*I;
Magnets(:,2,2) = Magnets(:,2,2).*I;
Magnets(:,3,2) = Magnets(:,3,2).*I;

%% ------------------------------------------------------------------------- % Отрисовка сетки

PointMagnetsROI(Node, Magnets, ROI)

%% ------------------------------------------------------------------------- % Полученное поле

B = ReceivedMagnetsField(Magnets, ROI, length(ROI(:,1)));

[X, Y, Z, Hx, Hy, Hz] = Transform3DMagnets(lx, ly, lz, nx, ny, nz, CenterROI, B(:,1), B(:,2), B(:,3));

[X1, Y1, Z1] = Display(X, Y, Z, Hx, nx, ny, nz, level, plane, "Полученное поле Bx");
[X1, Y1, Z1] = Display(X, Y, Z, Hy, nx, ny, nz, level, plane, "Полученное поле By");
[X1, Y1, Z1] = Display(X, Y, Z, Hz, nx, ny, nz, level, plane, "Полученное поле Bz");

DES = Bz - B(:,3);

krit = max(abs(DES)./max(Bz, [], 'all'),[],'all')*100;
disp(strcat("Максимальное отклонение: ", num2str(krit), ' %'));

toc

%% ------------------------------------------------------------------------- % Функции

function PointMagnetsROI(Node, Magnets, ROI) % Отрисовка точек, магнитных моментов и ROI

    figure('Name','Сетка','NumberTitle','off'); 
    movegui([0 560]);
    plot3(Node(:,1),Node(:,2),Node(:,3),'b.'); % Отрисовка узлов
    hold on
    xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]')
    axis equal
    grid on 
    %set(gca,'Visible', 'off')
    quiver3(Magnets(:,1,1), Magnets(:,2,1),Magnets(:,3,1),Magnets(:,1,2),Magnets(:,2,2),Magnets(:,3,2), 1,'r')
    plot3(ROI(:,1),ROI(:,2),ROI(:,3),'k.'); % Отрисовка ROI
    
end

function I = OptM(N, Magnets, ROI, Bx, By, Bz, gamma_x, gamma_y, gamma_z) % Создание массива узлов

    A = zeros(N,N);
    parfor j = 1:N
        for i = 1:N
            bz1 = CalcMagnetsbz(i, Magnets, ROI);
            bz2 = CalcMagnetsbz(j, Magnets, ROI);
            
            by1 = CalcMagnetsby(i, Magnets, ROI);
            by2 = CalcMagnetsby(j, Magnets, ROI);
            
            bx1 = CalcMagnetsbx(i, Magnets, ROI);
            bx2 = CalcMagnetsbx(j, Magnets, ROI);
            
            ArrZ = bz1.*bz2;
            ArrY = by1.*by2;
            ArrX = bx1.*bx2;
            
            A(i,j) = gamma_x*zum(ArrX) + gamma_y*zum(ArrY) + gamma_z*zum(ArrZ);
        end
    end

    clear bz1 bz2

    % Задание матрицы B

    B = zeros(N,1);

    for i = 1:N
        bz1 = CalcMagnetsbz(i, Magnets, ROI);
        by1 = CalcMagnetsby(i, Magnets, ROI);
        bx1 = CalcMagnetsbx(i, Magnets, ROI);
        
        ArrZ = Bz(:,1).*bz1;
        ArrY = By(:,1).*by1;
        ArrX = Bx(:,1).*bx1;
        B(i,1) = gamma_x*zum(ArrX) + gamma_y*zum(ArrY) + gamma_z*zum(ArrZ);
    end 

    I = pinv(A)*B;
end

function B = ReceivedMagnetsField(Magnets, ROI, K) % Расчет и отрисовка полученного поля

    B = zeros(K,3);
    Temp = B;

    for i=1:length(Magnets(:,1,1))
        Rx = ROI(:,1) - Magnets(i,1,1);
        Ry = ROI(:,2) - Magnets(i,2,1);
        Rz = ROI(:,3) - Magnets(i,3,1);
        Temp(:,1) = ((3*Rx.*ScalarP(Magnets(i,:,2),[Rx Ry Rz]) - Magnets(i,1,2).*RadiusV([Rx Ry Rz]).^2)./RadiusV([Rx Ry Rz]).^5)*1e-7;
        Temp(:,2) = ((3*Ry.*ScalarP(Magnets(i,:,2),[Rx Ry Rz]) - Magnets(i,2,2).*RadiusV([Rx Ry Rz]).^2)./RadiusV([Rx Ry Rz]).^5)*1e-7;
        Temp(:,3) = ((3*Rz.*ScalarP(Magnets(i,:,2),[Rx Ry Rz]) - Magnets(i,3,2).*RadiusV([Rx Ry Rz]).^2)./RadiusV([Rx Ry Rz]).^5)*1e-7;
        B = B + Temp;
    end

end

function [X1, Y1, Z1] = Display(X, Y, Z, H, nx, ny, nz, level, plane, Title) % Отрисовка поля

    %k = level;
    figure('Name',Title,'NumberTitle','off')
    movegui([1130 560]);
    if (plane == "XY") % ---------------------------------------------------------- XY
        k = find(Z(1,1,:)==level);
        X1 = X(:,:,k);
        Y1 = Y(:,:,k);
        Z1 = Z(:,:,k);
        contourf(X(:,:,k), Y(:,:,k), H(:,:,k), 20)
        xlabel ('x [m]'), ylabel ('y [m]'), title(strcat(Title,', z = ', num2str(level))) 
    elseif (plane == "XZ") % ------------------------------------------------------ XZ 
        k = find(Y(1,:,1)==level);
        XX = reshape(X(:,k,:),[nx,nz,1]);
        ZZ = reshape(Z(:,k,:),[nx,nz,1]);
        HH = reshape(H(:,k,:),[nx,nz,1]);
        X1 = X(:,k,:);
        Y1 = Y(:,k,:);
        Z1 = Z(:,k,:);
        contourf(XX, ZZ, HH, 20)
        xlabel ('x [m]'), ylabel ('z [m]'), title(strcat(Title,', y = ', num2str(level)))
    elseif (plane == "YZ") % ------------------------------------------------------ YZ
        k = find(X(:,1,1)==level);
        YY = reshape(Y(k,:,:),[ny,nz,1]);
        ZZ = reshape(Z(k,:,:),[ny,nz,1]);
        HH = reshape(H(k,:,:),[ny,nz,1]);
        X1 = X(k,:,:);
        Y1 = Y(k,:,:);
        Z1 = Z(k,:,:);
        contourf(YY, ZZ, HH, 20)
        xlabel ('y [m]'), ylabel ('z [m]'), title(strcat(Title,', x = ', num2str(level)))
    else
        disp('___ОШИБКА ВЫБОРА ПЛОСКОСТИ___');
    end
    colorbar
    box on
    hold on
    grid on

end

% -------------------------------------------------------------------------

function bz = CalcMagnetsbz(n, Magnets, ROI) % Параметр bz 
    bz = (1e-7)*(3*((Magnets(n,1,2).*(ROI(:,1)-Magnets(n,1,1)) + Magnets(n,2,2).*(ROI(:,2)-Magnets(n,2,1)) + Magnets(n,3,2).*(ROI(:,3)-Magnets(n,3,1))).*(ROI(:,3)-Magnets(n,3,1)))./(DistanceROI(Magnets(n,:,1), ROI).^5) - (Magnets(n,3,2)./(DistanceROI(Magnets(n,:,1), ROI).^3)));
end

function by = CalcMagnetsby(n, Magnets, ROI) % Параметр bz 
    by = (1e-7)*(3*((Magnets(n,1,2).*(ROI(:,1)-Magnets(n,1,1)) + Magnets(n,2,2).*(ROI(:,2)-Magnets(n,2,1)) + Magnets(n,3,2).*(ROI(:,3)-Magnets(n,3,1))).*(ROI(:,2)-Magnets(n,2,1)))./(DistanceROI(Magnets(n,:,1), ROI).^5) - (Magnets(n,2,2)./(DistanceROI(Magnets(n,:,1), ROI).^3)));
end

function bx = CalcMagnetsbx(n, Magnets, ROI) % Параметр bz 
    bx = (1e-7)*(3*((Magnets(n,1,2).*(ROI(:,1)-Magnets(n,1,1)) + Magnets(n,2,2).*(ROI(:,2)-Magnets(n,2,1)) + Magnets(n,3,2).*(ROI(:,3)-Magnets(n,3,1))).*(ROI(:,1)-Magnets(n,1,1)))./(DistanceROI(Magnets(n,:,1), ROI).^5) - (Magnets(n,1,2)./(DistanceROI(Magnets(n,:,1), ROI).^3)));
end

function di = DistanceROI(Node, ROI) % Расстояние между двумя узлами
    di = sqrt((Node(1)-ROI(:,1)).^2+(Node(2)-ROI(:,2)).^2+(Node(3)-ROI(:,3)).^2);
end

function [X, Y, Z, Bx, By, Bz] = Transform3DMagnets(lx, ly, lz, nx, ny, nz, Centre, Bx1, By1, Bz1) % Преобразование массива вида Node или ROI в трехмерный массив для построения
    [X, Y, Z] = ndgrid(linspace(-lx/2,lx/2,nx)+Centre(1),linspace(-ly/2,ly/2,ny)+Centre(2),linspace(-lz/2,lz/2,nz)+Centre(3));
    Bx = reshape(Bx1,[nx,ny,nz]);   
    By = reshape(By1,[nx,ny,nz]);  
    Bz = reshape(Bz1,[nx,ny,nz]);  
end

function p = ScalarP(Vector1, Vector2) % Скалярное произведение
    p = Vector1(:,1).*Vector2(:,1) + Vector1(:,2).*Vector2(:,2) + Vector1(:,3).*Vector2(:,3);
end

function R = RadiusV(A) % Радиус вектор точки
    R = sqrt(A(:,1).^2 + A(:,2).^2 + A(:,3).^2);
end



