clear all
close all
clc
BSmag = BSmag_init(); % Initialize BSmag analysis

%--------------------------------------------------------------------------
R = 0.07; % Радиус
R1 = 0.08;
R2 = 0.09;

%uiimport('1.txt');
data1 = readtable('YZ_1.txt');
data = table2array(data1)./1000;
YZ1 = zeros(length(data(:,1)),3);
YZ1(:,1) = data(:,2);
YZ1(:,2) = data(:,1);
YZ1 = transZ(YZ1,R);
YZ1 = flex(YZ1,R);

data1 = readtable('YZ_2.txt');
data = table2array(data1)./1000;
YZ2 = zeros(length(data(:,1)),3);
YZ2(:,1) = data(:,2);
YZ2(:,2) = data(:,1);
YZ2 = transZ(YZ2,R);
YZ2 = flex(YZ2,R);
%Blue = rotX(Blue, -pi/2);

data1 = readtable('YZ_3.txt');
data = table2array(data1)./1000;
YZ3 = zeros(length(data(:,1)),3);
YZ3(:,1) = data(:,2);
YZ3(:,2) = data(:,1);
YZ3 = transZ(YZ3,R);
YZ3 = flex(YZ3,R);
%Green = invY(Green);

data1 = readtable('YZ_4.txt');
data = table2array(data1)./1000;
YZ4 = zeros(length(data(:,1)),3);
YZ4(:,1) = data(:,2);
YZ4(:,2) = data(:,1);
YZ4 = transZ(YZ4,R);
YZ4 = flex(YZ4,R);
% Viol = rotX(Viol, pi/2);

data1 = readtable('YZ_5.txt');
data = table2array(data1)./1000;
YZ5 = zeros(length(data(:,1)),3);
YZ5(:,1) = data(:,2);
YZ5(:,2) = data(:,1);
YZ5 = transZ(YZ5,R);
YZ5 = flex(YZ5,R);

data1 = readtable('YZ_6.txt');
data = table2array(data1)./1000;
YZ6 = zeros(length(data(:,1)),3);
YZ6(:,1) = data(:,2);
YZ6(:,2) = data(:,1);
YZ6 = transZ(YZ6,R);
YZ6 = flex(YZ6,R);

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

data1 = readtable('Z_2.txt');
data = table2array(data1)./1000;
Z2 = zeros(length(data(:,1)),3);
Z2(:,1) = data(:,2);
Z2(:,2) = data(:,1);
Z2 = transZ(Z2,R1);
Z2 = flex(Z2,R1);

data1 = readtable('Z_3.txt');
data = table2array(data1)./1000;
Z3 = zeros(length(data(:,1)),3);
Z3(:,1) = data(:,2);
Z3(:,2) = data(:,1);
Z3 = transZ(Z3,R1);
Z3 = flex(Z3,R1);

data1 = readtable('Z_4.txt');
data = table2array(data1)./1000;
Z4 = zeros(length(data(:,1)),3);
Z4(:,1) = data(:,2);
Z4(:,2) = data(:,1);
Z4 = transZ(Z4,R1);
Z4 = flex(Z4,R1);

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

data1 = readtable('Y_2.txt');
data = table2array(data1)./1000;
Y2 = zeros(length(data(:,1)),3);
Y2(:,1) = data(:,2);
Y2(:,2) = data(:,1);
Y2 = transZ(Y2,R2);
Y2 = flex(Y2,R2);

data1 = readtable('Y_3.txt');
data = table2array(data1)./1000;
Y3 = zeros(length(data(:,1)),3);
Y3(:,1) = data(:,2);
Y3(:,2) = data(:,1);
Y3 = transZ(Y3,R2);
Y3 = flex(Y3,R2);

data1 = readtable('Y_4.txt');
data = table2array(data1)./1000;
Y4 = zeros(length(data(:,1)),3);
Y4(:,1) = data(:,2);
Y4(:,2) = data(:,1);
Y4 = transZ(Y4,R2);
Y4 = flex(Y4,R2);

Y1 = rotX(Y1, -0.2*pi/180);
Y2 = rotX(Y2, -0.2*pi/180);
Y3 = rotX(Y3, -0.2*pi/180);
Y4 = rotX(Y4, -0.2*pi/180);

% ------------------------------------------------------- 

app = 0; % Аппроксимация: 1 - вкл; 0 - выкл
Grad = 0; % 0 - X; 1 - Y; 2 - Z; 3 - XY; 4 - XZ; 5 - YZ

plane = 2; % Выбор плоскости 0 - XY; 1 - XZ; 2 - YZ

level_plane = 0.0; % Выбор начального уровня (Для слайдшоу автоматически)

movie = 0; % 1 - по срезам; 0 - поле на уровне level_plane
nm = 11; % Кол-во срезов

Iyz = 10; % Величина тока в катушках
Iz = -5;
Iy = -5;

b = 0; % Угол повората в радианах

% Размеры области расчета
fovX = 0.05;
fovY = 0.05;
fovZ = 0.05;

% Кол-во точек 
fovP = 21;

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

dGamma1 = 1e-3; % filament max discretization step [mm]

[BSmag] = BSmag_add_filament(BSmag,YZ1,Iyz,dGamma1);
[BSmag] = BSmag_add_filament(BSmag,YZ2,Iyz,dGamma1);
[BSmag] = BSmag_add_filament(BSmag,YZ3,Iyz,dGamma1);
[BSmag] = BSmag_add_filament(BSmag,YZ4,Iyz,dGamma1);
[BSmag] = BSmag_add_filament(BSmag,YZ5,Iyz,dGamma1);
[BSmag] = BSmag_add_filament(BSmag,YZ6,Iyz,dGamma1);

[BSmag] = BSmag_add_filament(BSmag,Z1,Iz,dGamma1);
[BSmag] = BSmag_add_filament(BSmag,Z2,Iz,dGamma1);
[BSmag] = BSmag_add_filament(BSmag,Z3,Iz,dGamma1);
[BSmag] = BSmag_add_filament(BSmag,Z4,Iz,dGamma1);

[BSmag] = BSmag_add_filament(BSmag,Y1,Iy,dGamma1);
[BSmag] = BSmag_add_filament(BSmag,Y2,Iy,dGamma1);
[BSmag] = BSmag_add_filament(BSmag,Y3,Iy,dGamma1);
[BSmag] = BSmag_add_filament(BSmag,Y4,Iy,dGamma1);

%--------------------------------------------------------------------------

% Точки (где хотим считать поле)
x_M = linspace(-fovX/2, fovX/2, fovP); % x [m]
y_M = linspace(-fovY/2, fovY/2, fovP); % y [m]
z_M = linspace(-fovZ/2, fovZ/2, fovP); % z [m]

% Расчет-------------------------------------------------------------------

if plane == 0 % XY --------------------------------------------------------
    if movie
        level_plane = -fovZ/2;
        [X_M,Y_M] = ndgrid(x_M,y_M);
        Z_M = zeros(fovP,fovP)+level_plane;
        dd = fovZ/(nm-1);
        for m = 1:nm
            [BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,X_M,Y_M,Z_M); % Био-савар

            figure(1)
            quiver3(X,Y,Z,BX,BY,BZ,'b')

            figure(2) % Отображение
            box on
            grid on
            contourf(X, Y, BZ, 20) 
            xlabel ('x [m]'), ylabel ('y [m]'), title (strcat('Bz [T], z = ', num2str(level_plane + dd*(m-1)))) 
            colorbar
            axis equal
            %saveas(gcf,strcat(num2str(m),'.png')) % Сохранение
            Z_M = Z_M + dd;
        end
    else
        [X_M,Y_M] = ndgrid(x_M,y_M);
        Z_M = zeros(fovP,fovP)+level_plane;
        
        [BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,X_M,Y_M,Z_M); % Био-савар

        figure(1)
        quiver3(X,Y,Z,BX,BY,BZ,'b')

        figure(2) % Отображение
        box on
        grid on
        contourf(X, Y, BZ,20) 
        xlabel ('x [m]'), ylabel ('y [m]'), title (strcat('Bz [T], z = ', num2str(level_plane))) 
        colorbar
        axis equal
    end
    
elseif plane == 1 % XZ ----------------------------------------------------
    
    if movie
        level_plane = -fovZ/2;
        [X_M,Z_M] = ndgrid(x_M,z_M);
        Y_M = zeros(fovP,fovP)+level_plane;
        dd = fovY/(nm-1);
        for m = 1:nm
            [BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,X_M,Y_M,Z_M); % Био-савар

            figure(1)
            quiver3(X,Y,Z,BX,BY,BZ,'b')

            figure(2) % Отображение
            box on
            grid on
            contourf(X, Z, BZ,20) 
            xlabel ('x [m]'), ylabel ('z [m]'), title (strcat('Bz [T], y = ', num2str(level_plane + dd*(m-1)))) 
            colorbar
            axis equal
            %saveas(gcf,strcat(num2str(m),'.png')) % Сохранение
            Y_M = Y_M + dd;
        end
    else
        [X_M,Z_M] = ndgrid(x_M,z_M);
        Y_M = zeros(fovP,fovP)+level_plane;
        
        [BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,X_M,Y_M,Z_M); % Био-савар

        figure(1)
        quiver3(X,Y,Z,BX,BY,BZ,'b')

        figure(2) % Отображение
        box on
        grid on
        contourf(X, Z, BZ,20) 
        xlabel ('x [m]'), ylabel ('z [m]'), title (strcat('Bz [T], y = ', num2str(level_plane))) 
        colorbar 
        axis equal
    end
    
elseif plane == 2 % YZ ----------------------------------------------------
    
    if movie
        level_plane = -fovZ/2;
        [Y_M,Z_M] = ndgrid(y_M,z_M);
        X_M = zeros(fovP,fovP)+level_plane;
        dd = fovX/(nm-1);
        for m = 1:nm
            [BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,X_M,Y_M,Z_M); % Био-савар

            figure(1)
            quiver3(X,Y,Z,BX,BY,BZ,'b')

            figure(2) % Отображение
            box on
            grid on
            contourf(Y, Z, BZ, 20) 
            xlabel ('y [m]'), ylabel ('z [m]'), title (strcat('Bz [T], x = ', num2str(level_plane + dd*(m-1)))) 
            colorbar
            axis equal
            %saveas(gcf,strcat(num2str(m),'.png')) % Сохранение
            X_M = X_M + dd;
        end
    else
        [Y_M,Z_M] = ndgrid(y_M,z_M);
        X_M = zeros(fovP,fovP)+level_plane;
        
        [BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,X_M,Y_M,Z_M); % Био-савар

        figure(1)
        quiver3(X,Y,Z,BX,BY,BZ,'b')

        figure(2) % Отображение
        box on
        grid on
        contourf(Y, Z, BZ, 20) 
        %surf(Y, Z, BZ)
        xlabel ('y [m]'), ylabel ('z [m]'), title (strcat('Bz [T], x = ', num2str(level_plane))) 
        colorbar 
        axis equal
    end
    
else 
    ERROR = 404
end
%--------------------------------------------------------------------------

% "Аппроксимация"
if app
    switch Grad
        case 0 % X --------------------------------------------------------------------------------------
            
            G = [-0.1:0.0001:0.1];
            B0 = [-0.0005:0.00001:0.0005];

            Bt = G(1)*X+B0(1);
            E = 100*(Bt-BZ)./max(max(Bt));
            a = max(max(abs(E)));
            G0 = G(1);
            B00 = B0(1);

            for i = B0
                for j = G
                    Bt = j*X+i;
                    E = 100*(Bt-BZ)./max(max(Bt));
                    if max(max(abs(E))) < a
                        G0 = j;
                        B00 = i;
                        a = max(max(abs(E)));
                    end
                end
            end

            disp(strcat('Градиент Gx:', ' ', num2str(G0), ' Т/м'))
            disp(strcat('Фон:', ' ', num2str(B00)))
            disp(strcat('Отклонение:', ' ',num2str(a),' %'))
            
            Bt = G0*X + B00;
            E = 100*(Bt-BZ)./max(max(Bt));
            
            figure(4)
            box on
            grid on
            if plane == 0
                surf(X, Y, E) 
                xlabel ('x [m]'), ylabel ('y [m]'), title ('Error')
            elseif plane == 1 
                surf(X, Z, E) 
                xlabel ('x [m]'), ylabel ('z [m]'), title ('Error')
            elseif plane == 2 
                surf(Y, Z, E) 
                xlabel ('y [m]'), ylabel ('z [m]'), title ('Error')
            else
                ERROR = 404
            end
            colorbar 
            
        case 1 % Y --------------------------------------------------------------------------------------
            
            G = [-0.1:0.0001:0.1];
            B0 = [-0.0005:0.00001:0.0005];

            Bt = G(1)*Y+B0(1);
            E = 100*(Bt-BZ)./max(max(Bt));
            a = max(max(abs(E)));
            G0 = G(1);
            B00 = B0(1);

            for i = B0
                for j = G
                    Bt = j*Y+i;
                    E = 100*(Bt-BZ)./max(max(Bt));
                    if max(max(abs(E))) < a
                        G0 = j;
                        B00 = i;
                        a = max(max(abs(E)));
                    end
                end
            end

            disp(strcat('Градиент Gy:', ' ', num2str(G0), ' Т/м'))
            disp(strcat('Фон:', ' ', num2str(B00)))
            disp(strcat('Отклонение:', ' ',num2str(a),' %'))
            
            Bt = G0*Y + B00;
            E = 100*(Bt-BZ)./max(max(Bt));
            
            figure(4)
            box on
            grid on
            if plane == 0
                surf(X, Y, E) 
                xlabel ('x [m]'), ylabel ('y [m]'), title ('Error')
            elseif plane == 1 
                surf(X, Z, E) 
                xlabel ('x [m]'), ylabel ('z [m]'), title ('Error')
            elseif plane == 2 
                surf(Y, Z, E) 
                xlabel ('y [m]'), ylabel ('z [m]'), title ('Error')
            else
                ERROR = 404
            end
            colorbar
            
        case 2 % Z  --------------------------------------------------------------------------------------
            
            G = [-0.1:0.0001:0.1];
            B0 = [-0.0005:0.00001:0.0005];

            Bt = G(1)*Z+B0(1);
            E = 100*(Bt-BZ)./max(max(Bt));
            a = max(max(abs(E)));
            G0 = G(1);
            B00 = B0(1);

            for i = B0
                for j = G
                    Bt = j*Z+i;
                    E = 100*(Bt-BZ)./max(max(Bt));
                    if max(max(abs(E))) < a
                        G0 = j;
                        B00 = i;
                        a = max(max(abs(E)));
                    end
                end
            end

            disp(strcat('Градиент Gz:', ' ', num2str(G0), ' Т/м'))
            disp(strcat('Фон:', ' ', num2str(B00)))
            disp(strcat('Отклонение:', ' ',num2str(a),' %'))
            
            Bt = G0*Z + B00;
            E = 100*(Bt-BZ)./max(max(Bt));
            
            figure(4)
            box on
            grid on
            if plane == 0
                surf(X, Y, E) 
                xlabel ('x [m]'), ylabel ('y [m]'), title ('Error')
            elseif plane == 1 
                surf(X, Z, E) 
                xlabel ('x [m]'), ylabel ('z [m]'), title ('Error')
            elseif plane == 2 
                surf(Y, Z, E) 
                xlabel ('y [m]'), ylabel ('z [m]'), title ('Error')
            else
                ERROR = 404
            end
            colorbar
            
        case 3 % XY  --------------------------------------------------------------------------------------
            
            G = [-0.2:0.0001:0.2];
            B0 = [-0.0005:0.00001:0.0005];

            Bt = G(1).*X.*Y+B0(1);
            E = 100*(Bt-BZ)./max(max(Bt));
            a = max(max(abs(E)));
            G0 = G(1);
            B00 = B0(1);

            for i = B0
                for j = G
                    Bt = j.*X.*Y+i;
                    E = 100*(Bt-BZ)./max(max(Bt));
                    if max(max(abs(E))) < a
                        G0 = j;
                        B00 = i;
                        a = max(max(abs(E)));
                    end
                end
            end

            disp(strcat('Градиент Gxy:', ' ', num2str(G0), ' Т/м2'))
            disp(strcat('Фон:', ' ', num2str(B00)))
            disp(strcat('Отклонение:', ' ',num2str(a),' %'))
            
            Bt = G0.*X.*Y + B00;
            E = 100*(Bt-BZ)./max(max(Bt));
            
            figure(4)
            box on
            grid on
            if plane == 0
                surf(X, Y, E) 
                xlabel ('x [m]'), ylabel ('y [m]'), title ('Error')
            elseif plane == 1 
                surf(X, Z, E) 
                xlabel ('x [m]'), ylabel ('z [m]'), title ('Error')
            elseif plane == 2 
                surf(Y, Z, E) 
                xlabel ('y [m]'), ylabel ('z [m]'), title ('Error')
            else
                ERROR = 404
            end
            colorbar
            
        case 4 % XZ  --------------------------------------------------------------------------------------
            
            G = [-0.2:0.0001:0.2];
            B0 = [-0.0005:0.00001:0.0005];

            Bt = G(1).*X.*Z+B0(1);
            E = 100*(Bt-BZ)./max(max(Bt));
            a = max(max(abs(E)));
            G0 = G(1);
            B00 = B0(1);

            for i = B0
                for j = G
                    Bt = j.*X.*Z+i;
                    E = 100*(Bt-BZ)./max(max(Bt));
                    if max(max(abs(E))) < a
                        G0 = j;
                        B00 = i;
                        a = max(max(abs(E)));
                    end
                end
            end

            disp(strcat('Градиент Gxz:', ' ', num2str(G0), ' Т/м2'))
            disp(strcat('Фон:', ' ', num2str(B00)))
            disp(strcat('Отклонение:', ' ',num2str(a),' %'))
            
            Bt = G0.*X.*Z + B00;
            E = 100*(Bt-BZ)./max(max(Bt));
            
            figure(4)
            box on
            grid on
            if plane == 0
                surf(X, Y, E) 
                xlabel ('x [m]'), ylabel ('y [m]'), title ('Error')
            elseif plane == 1 
                surf(X, Z, E) 
                xlabel ('x [m]'), ylabel ('z [m]'), title ('Error')
            elseif plane == 2 
                surf(Y, Z, E) 
                xlabel ('y [m]'), ylabel ('z [m]'), title ('Error')
            else
                ERROR = 404
            end
            colorbar
            
        case 5 % YZ  --------------------------------------------------------------------------------------
            
            G = [-0.2:0.0001:0.2];
            B0 = [-0.0005:0.00001:0.0005];

            Bt = G(1).*Y.*Z+B0(1);
            E = 100*(Bt-BZ)./max(max(Bt));
            a = max(max(abs(E)));
            G0 = G(1);
            B00 = B0(1);

            for i = B0
                for j = G
                    Bt = j.*Y.*Z+i;
                    E = 100*(Bt-BZ)./max(max(Bt));
                    if max(max(abs(E))) < a
                        G0 = j;
                        B00 = i;
                        a = max(max(abs(E)));
                    end
                end
            end

            disp(strcat('Градиент Gyz:', ' ', num2str(G0), ' Т/м2'))
            disp(strcat('Фон:', ' ', num2str(B00)))
            disp(strcat('Отклонение:', ' ',num2str(a),' %'))
            
            Bt = G0.*Y.*Z + B00;
            E = 100*(Bt-BZ)./max(max(Bt));
            
            figure(4)
            box on
            grid on
            if plane == 0
                surf(X, Y, E) 
                xlabel ('x [m]'), ylabel ('y [m]'), title ('Error')
            elseif plane == 1 
                surf(X, Z, E) 
                xlabel ('x [m]'), ylabel ('z [m]'), title ('Error')
            elseif plane == 2 
                surf(Y, Z, E) 
                xlabel ('y [m]'), ylabel ('z [m]'), title ('Error')
            else
                ERROR = 404
            end
            colorbar
            
        otherwise
            disp('Градиент не указан')
    end
else 
    disp('Аппроксимация отключена')
end

%--------------------------------------------------------------------------

function [Coil, S] = wend(n_points, n_vit, R, alfa, len, h, dy, dz, sm_points)

S = 0; % Длина катушки 

alfaI = alfa/2;
alfaF = -alfa/2;

Coil(1,1) = len/2;
Coil(1,2) = R*cos(alfaF); 
Coil(1,3) = R*sin(alfaF);

last_point = 0;
for j = 1:n_vit
    
    R1 = R;
    Y = 0;
    Z = 0;
    
    % Изменение угла на виток 
    if j > 1
        alfaF = alfaF + 2*h(j)/R;
    end
    dalfa = (alfaI-alfaF)/(n_points);
    
    % Отрисовка первой дуги
    for i = 1:(n_points+1)
    
        Coil(i+last_point,1) = len/2-h(j)*(j-1);
        Coil(i+last_point,2) = R*cos(alfaI - dalfa*(i-1));
        Coil(i+last_point,3) = R*sin(alfaI - dalfa*(i-1));
    
    end

    last_point = last_point + i;
    dx = (len-2*(j-1)*h(j))/(n_points+1);
    
    % Отрисовка первой прямой
    for i = 1:(n_points+1)
    
        Coil(i+last_point,1) = len/2-h(j)*(j-1)-dx*i;
        Coil(i+last_point,2) = Coil(i+last_point-1,2)+dz;
        Coil(i+last_point,3) = Coil(i+last_point-1,3)-dy;
        Y = Y + dy;
        Z = Z + dz;
    end

    last_point = last_point + i;
    
    R = R1 + Y;
    
    if j > 1
        alfaI = alfaI - 2*h(j)/R;
    end
    
    alfaI1 = alfaI;
    alfaF1 = alfaF;
    alfaI = alfaI + (atan(tan(alfaI)+Y/(R*cos(alfaI)))-alfaI);
    alfaF = alfaF - (atan(tan(alfaF)+Y/(R*cos(alfaF)))-alfaF);
    
    dalfa = (alfaI-alfaF)/(n_points);    
    
    % Отрисовка второй дуги
    for i = 1:(n_points)
    
        Coil(i+last_point,1) = Coil(i+last_point-1,1);
        Coil(i+last_point,2) = R*cos(alfaF + dalfa*(i))+Z;
        Coil(i+last_point,3) = R*sin(alfaF + dalfa*(i));
    
    end

    last_point = last_point + i;
    
    if j>1
        k = 2*j;
    else 
        k = j;
    end
    
    dx = (len-k*h(j))/(n_points+1);
    
    % Отрисовка второй прямой
    for i = 1:(n_points)
    
        Coil(i+last_point,1) = -(len/2-h(j)*(j-1))+dx*i;
        Coil(i+last_point,2) = Coil(i+last_point-1,2)-dz;
        Coil(i+last_point,3) = Coil(i+last_point-1,3)-dy;
    
    end
    
    last_point = last_point + i;
    R = R1;
    alfaI = alfaI1;
    alfaF = alfaF1;
end

% Сглаживание
Coil(:,1) = smooth(Coil(:,1),sm_points);
Coil(:,2) = smooth(Coil(:,2),sm_points);
Coil(:,3) = smooth(Coil(:,3),sm_points);


for i = 2:length(Coil(:,1)) % Определение длины
    ds = sqrt((Coil(i-1,1)-Coil(i,1))^2 + (Coil(i-1,2)-Coil(i,2))^2 + (Coil(i-1,2)-Coil(i,2))^2);
    S = S + ds;
end

end

function Coil = rotX(Coil, a)
    M = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)];
    for i = 1:length(Coil(:,1))
        Coil(i,:) = Coil(i,:)*M;
    end
end

function Coil = rotY(Coil, a)
    M = [cos(a) 0 sin(a); 0 1 0; -sin(a) 0 cos(a)];
    for i = 1:length(Coil(:,1))
        Coil(i,:) = Coil(i,:)*M;
    end
end

function Coil = rotZ(Coil, a)
    M = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
    for i = 1:length(Coil(:,1))
        Coil(i,:) = Coil(i,:)*M;
    end
end

function Coil = invX(Coil)
    Coil(:,1) = -Coil(:,1);
end

function Coil = invY(Coil)
    Coil(:,2) = -Coil(:,2);
end

function Coil = invZ(Coil)
    Coil(:,3) = -Coil(:,3);
end

function Coil = transX(Coil, X)
    Coil(:,1) = Coil(:,1) + X;
end

function Coil = transY(Coil, Y)
    Coil(:,2) = Coil(:,2) + Y;
end

function Coil = transZ(Coil, Z)
    Coil(:,3) = Coil(:,3) + Z;
end

function [Coil, S] = wendPlanar(n_points, n_vit, width, len, h1, h2, h3, h4, sm_points, pa, pb, pc, pd)

    S = 0;
    
    Coil(1,1) = len/2;
    Coil(1,2) = width/2; 
    Coil(1,3) = 0;

    last_point = 2;
    for j = 1:n_vit

        %---------------------------- Первый торец

        if (j == 1)
            dy = (width)/n_points;
        else
            dy = (width-h1(j)*(j-1)-h3(j)*(j-2))/n_points;
        end

        for i = last_point:n_points+last_point-1
            Coil(i,1) = Coil(i-1,1);
            Coil(i,2) = Coil(i-1,2)-dy;
            Coil(i,3) = Coil(i-1,3);
            Coil(i,1) = Coil(i-1,1)-pa*(i-(50+last_point))*dy/101;
        end

        last_point = i+1;

        %---------------------------- Второй торец

        if (j == 1)
            dx = (len)/n_points;
        else
            dx = (len-(h2(j)+h4(j))*(j-1))/n_points;
        end

        for i = last_point:n_points+last_point-1
            Coil(i,1) = Coil(i-1,1)-dx;
            Coil(i,2) = Coil(i-1,2);
            Coil(i,3) = Coil(i-1,3);
            Coil(i,2) = Coil(i-1,2)-pb*(i-(50+last_point))*dx/101;
        end

        last_point = i+1;

        %---------------------------- Третий торец

        if (j == 1)
            dy = (width)/n_points;
        else
            dy = (width-(h1(j)+h3(j))*(j-1))/n_points;
        end

        for i = last_point:n_points+last_point-1
            Coil(i,1) = Coil(i-1,1);
            Coil(i,2) = Coil(i-1,2)+dy;
            Coil(i,3) = Coil(i-1,3);
            Coil(i,1) = Coil(i-1,1)+pc*(i-(50+last_point))*dy/101;
        end

        last_point = i+1;

        %---------------------------- Четвертый торец

        if (j == 1)
            dx = (len-h4(j))/n_points;
        else
            dx = (len-(h2(j)*(j-1)+h4(j)*j))/n_points;
        end

        for i = last_point:n_points+last_point-1
            Coil(i,1) = Coil(i-1,1)+dx;
            Coil(i,2) = Coil(i-1,2);
            Coil(i,3) = Coil(i-1,3);
            Coil(i,2) = Coil(i-1,2)+pd*(i-(50+last_point))*dx/101;
        end

        last_point = i+1;
    end
    
    Coil(:,1) = smooth(Coil(:,1),sm_points);
    Coil(:,2) = smooth(Coil(:,2),sm_points);
    Coil(:,3) = smooth(Coil(:,3),sm_points);
    
    for i = 2:length(Coil(:,1)) % Определение длины
        ds = sqrt((Coil(i-1,1)-Coil(i,1))^2 + (Coil(i-1,2)-Coil(i,2))^2 + (Coil(i-1,2)-Coil(i,2))^2);
        S = S + ds;
    end
end

function [Coil, S] = wendSpiral(n_points, n_vit, R, h)

S = 0; % Длина катушки 

Coil(1,1) = R;
Coil(1,2) = 0; 
Coil(1,3) = 0;
    
R1 = R;
Y = 0;
Z = 0;
alfa = 0;
dalfa = 2*pi/n_points;

for i = 2:(n_vit*n_points)
    
    Coil(i,1) = R*(alfa*h(1)/(2*pi*R)+1)*cos(alfa);
    Coil(i,2) = R*(alfa*h(1)/(2*pi*R)+1)*sin(alfa); 
    Coil(i,3) = 0;
    alfa = alfa+dalfa;
    
end

for i = 2:length(Coil(:,1)) % Определение длины
    ds = sqrt((Coil(i-1,1)-Coil(i,1))^2 + (Coil(i-1,2)-Coil(i,2))^2 + (Coil(i-1,2)-Coil(i,2))^2);
    S = S + ds;
end

end

function Coil = flex(Coil,R)
    Len = length(Coil(:,1));
    for i = 1:(Len)
        y = Coil(i,2);
        %deltaZ = y*sin(asin(y/(2*R)));
        %deltaY = y*(1-cos(asin(y/(2*R))));
        Coil(i,1)=Coil(i,1);
        %Coil(i,2)=Coil(i,2)-deltaY;%*sign(Coil(i,2));
        %Coil(i,3)=Coil(i,3)-deltaZ;%*sign(Coil(i,2));
        Coil(i,2)=R*sin(y/R);
        Coil(i,3)=R*(cos(y/R));
    end
    Coil = real(Coil);
end
