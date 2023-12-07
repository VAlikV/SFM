function PrintFieldSphere(lx, ly, lz, nx, ny, nz, centerFOV, Bdes, level, plane, Rsp, centerROI, Title, name) % Отрисовка целевого поля

[X,Y,Z,B] = Transform3D(lx, ly, lz, nx, ny, nz, centerFOV, Bdes); % Для отрисовки поля

figure('Name', Title ,'NumberTitle','off');
if (plane == "XY") % ---------------------------------------------------------- XY
    k = find(abs(Z(1,1,:)-level) < 0.00001);

    r = sqrt(Rsp^2 - (level-centerROI(3))^2);
    theta=linspace(0,2*pi,200); 
    x=r*cos(theta) + centerROI(1); 
    y=r*sin(theta) + centerROI(2);

    contourf(X(:,:,k)*1000, Y(:,:,k)*1000, B(:,:,k), 20)
    hold("on")
    plot(x*1000,y*1000,'--k');
    xlabel('x, мм'), ylabel('y, мм'), title(strcat(name,', z = ', num2str(level*1000), 'мм')) 
elseif (plane == "XZ") % ------------------------------------------------------ XZ 
    k = find(abs(Y(1,:,1)-level) < 0.00001);
    
    r=sqrt(Rsp^2 - (level-centerROI(2))^2);
    theta=linspace(0,2*pi,200); 
    x=r*cos(theta) + centerROI(1); 
    y=r*sin(theta) + centerROI(3);
    
    XX = reshape(X(:,k,:),[nx,nz,1]);
    ZZ = reshape(Z(:,k,:),[nx,nz,1]);
    BB = reshape(B(:,k,:),[nx,nz,1]);
    
    contourf(XX*1000, ZZ*1000, BB, 20)
    hold("on")
    plot(x*1000,y*1000,'--k');
    xlabel('x, мм'), ylabel('z, мм'), title(strcat(name,', y = ', num2str(level*1000), 'мм'))
elseif (plane == "YZ") % ------------------------------------------------------ YZ
    k = find(abs(X(:,1,1)-level) < 0.00001);
    
    r=sqrt(Rsp^2 - (level-centerROI(1))^2);
    theta=linspace(0,2*pi,200); 
    x=r*cos(theta) + centerROI(2); 
    y=r*sin(theta) + centerROI(3);
    
    YY = reshape(Y(k,:,:),[ny,nz,1]);
    ZZ = reshape(Z(k,:,:),[ny,nz,1]);
    BB = reshape(B(k,:,:),[ny,nz,1]);
    
    contourf(YY*1000, ZZ*1000, BB, 20)
    hold("on")
    plot(x*1000,y*1000,'--k');
    xlabel('y, мм'), ylabel('z, мм'), title(strcat(name,', x = ', num2str(level*1000), 'мм'))
else
    disp('___ОШИБКА ВЫБОРА ПЛОСКОСТИ___');
end
colorbar
axis('equal')
box("on")
grid("on") 
hold("off")
