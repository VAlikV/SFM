function PrintField(lx, ly, lz, nx, ny, nz, centerROI, Bz, level, plane) % Отрисовка полученного поля 
    
    [X1,Y1,Z1,Bz] = Transform3D(lx, ly, lz, nx, ny, nz, centerROI, Bz);

    figure('Name','Полученное поле','NumberTitle','off'); 
    %movegui([570 30]);
    if (plane == 0) % ---------------------------------------------------------- XY
        k = find(abs(Z1(1,1,:)-level) < 0.00001);
        contourf(X1(:,:,k), Y1(:,:,k), Bz(:,:,k), 20)
        xlabel ('x [m]'), ylabel ('y [m]'), title(strcat('Bz [T], z = ', num2str(level))) 
    elseif (plane == 1) % ------------------------------------------------------ XZ
        k = find(abs(Y1(1,:,1)-level) < 0.00001);
        XX1 = reshape(X1(:,k,:),[nx,nz,1]);
        ZZ1 = reshape(Z1(:,k,:),[nx,nz,1]);
        BB1 = reshape(Bz(:,k,:),[nx,nz,1]);
        contourf(XX1, ZZ1, BB1, 20)
        xlabel ('x [m]'), ylabel ('z [m]'), title(strcat('Bz [T], y = ', num2str(level)))
    elseif (plane == 2) % ------------------------------------------------------ YZ
        k = find(abs(X1(:,1,1)-level) < 0.00001);
        YY1 = reshape(Y1(k,:,:),[ny,nz,1]);
        ZZ1 = reshape(Z1(k,:,:),[ny,nz,1]);
        BB1 = reshape(Bz(k,:,:),[ny,nz,1]);
        contourf(YY1, ZZ1, BB1, 20)
        xlabel ('y [m]'), ylabel ('z [m]'), title(strcat('Bz [T], x = ', num2str(level)))
    else
        disp('___ОШИБКА ВЫБОРА ПЛОСКОСТИ___');
    end
    colorbar
    axis equal
    box on
    grid on

end