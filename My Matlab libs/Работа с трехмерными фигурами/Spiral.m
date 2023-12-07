function [Coil, S] = Spiral(n_points, n_vit, R, L)
% n_points - кол-во точек в витке
% n_vit - кол-во витков
S = 0;                  % Длина катушки 
N = n_points*n_vit;     % Кол-во точек в катушке

Coil = zeros(N,3);

fi = linspace(0, 2*pi*(n_vit-1), N);
len = linspace(0, L, N);

Coil(:,1) = R*cos(fi);
Coil(:,2) = R*sin(fi);
Coil(:,3) = len;

for i = 2:length(Coil(:,1)) % Определение длины
    ds = sqrt((Coil(i-1,1)-Coil(i,1))^2 + (Coil(i-1,2)-Coil(i,2))^2 + (Coil(i-1,2)-Coil(i,2))^2);
    S = S + ds;
end

end
