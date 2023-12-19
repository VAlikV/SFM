clc
clear all
close all
% Метод функции потока для цилиндрических катушек
%% ------------------------------------------------------------------------- % Задание параметров

tic

% Размеры
R = 0.05; % Радиус катушки
nfi = 20;  % Количество точек окружности

L = 0.05; % Половина длины катушки (от 0 в одну сторону)
nL = 11;  % Количество точек по длинне

axial = 0;

%% ------------------------------------------------------------------------- % Создание массива узлов

ROI = CreateNode(R, L, nfi, nL, axial); % Создание массива узлов

figure(1)
plot3(ROI(:,1),ROI(:,2),ROI(:,3),'k.');
xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]')
axis equal
grid on 

%% ------------------------------------------------------------------------- % ROI

data1 = readtable('B_target.txt');
Bt = 0.153 - table2array(data1)./1000;

% [file,path3] = uigetfile('*.txt');
% path3 = strcat(path3,file);
% fileID = fopen(path3,'r');
% 
% if fileID ~= -1
%     data = table2array(readtable(path3));
%     if length(data(1,:)) == 3
%         A = data;
%     end          
% end
            
%% ------------------------------------------------------------------------- % Функции

function Node = CreateNode(R, L, nfi, nL, axial) % Создание массива узлов
    Node = zeros(nfi,3); % Сетка
    Temp = zeros(nfi,3); 

    f = linspace(-pi+2*pi/nfi,pi,nfi); % Массив значений угла
    l = linspace(-L,L,nL); % Массив значений длины

    Node(:,1) = l(1);      % X
    Node(:,2) = R*sin(f);  % Y
    Node(:,3) = R*cos(f);  % Z

    for i = 2:nL
        Temp(:,1) = l(i); % X
        Temp(:,2) = R*sin(f); % Y        
        Temp(:,3) = R*cos(f); % Z
        Node = [Node; Temp];
    end

    if axial
        Node = rotY(Node, pi/2); % Поворот катушки
    end
end






