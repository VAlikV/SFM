function [Lx, Ly, BB] = Transform2D(lx, ly, nx, ny, B) % Преобразование массива вида Node или ROI в двумерный массив для построения
    [Lx, Ly] = meshgrid(linspace(-lx/2,lx/2,nx), linspace(-ly/2,ly/2,ny));
    BB = reshape(B,[nx,ny]);
end