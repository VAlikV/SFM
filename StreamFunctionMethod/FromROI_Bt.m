clc 
clear all
close all

CubeROI = struct('nx', 8, 'ny', 8, 'nz', 8, 'lx', 0.05, 'ly', 0.05, 'lz', 0.01, 'center', [0 0 0]);

ROI = CreateCubeROI(CubeROI.nx, CubeROI.ny, CubeROI.nz, CubeROI.lx, CubeROI.ly, CubeROI.lz, CubeROI.center);

G = 1;
Boff = 0;

B_target = G*(ROI(:,2).*ROI(:,2)) + Boff;

save("ROI.mat","ROI","-ascii")

save("B_t.txt","B_target","-ascii")

figure('Name', 'Сетка' ,'NumberTitle','off');
plot3(ROI(:,1),ROI(:,2),ROI(:,3),'k.'); % Отрисовка ROI
hold( "on")
xlabel('x, мм'), ylabel('y, мм'), zlabel('z, мм')
axis("equal")
grid( "on")

%pause(3)

%S = load("Out_I__Opt_150_66.5_15_11_1_15_5.txt", "-ascii");
