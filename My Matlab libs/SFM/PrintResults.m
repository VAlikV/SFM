function [px, py, Fi, Ll, II, Red1, Blue1, Red2, Blue2] = PrintResults( R, L, nfi, nL, axial, I, Nc, eb, er, zeror) % Инетерполяция потоковой функции и контур с током

[Fi, Ll, II] = TransformICilinder(L, nfi, nL, I);
[px,py] = gradient(II);
p = max(sqrt(px.^2+py.^2), [],'all');
px = px/(p);
py = py/(p);
            
[Fi_int, L_int, I_int, Red, Blue, Red1, Blue1, Red2, Blue2] = InterpCircuit(R, L, nfi, nL, I, Nc, eb, er, zeror);
W = 2*pi*R*1000;
LL = L*1000;

% "Потоковая функция"
figure('Name', "Потоковая функция" ,'NumberTitle','off');
contourf(Fi,Ll*1000,II);
box("on")
xlabel('fi'), ylabel('L, мм')
colorbar
%title( "Потоковая функция")

% "Потоковая функция + Grad"
figure('Name', "Потоковая функция + Grad" ,'NumberTitle','off');
contourf(Fi,Ll,II);
box("on")
hold("on")
quiver(Fi,Ll,px,py,0.2,'*r', 'MarkerSize', 2)
xlabel('fi'), ylabel('L, м')
colorbar
%title("Потоковая функция + Grad")
hold("off")

% "Намотка 2D"
figure('Name', "Намотка 2D" ,'NumberTitle','off');
if ~isempty(Red2)
    plot(Red2(:,2)*1000, Red2(:,1)*1000, '.r','MarkerSize',1)
    hold("on")
end
if ~isempty(Blue2)
    plot(Blue2(:,2)*1000, Blue2(:,1)*1000, '.b','MarkerSize',1)
end

if axial
    xlabel('r, мм'), ylabel('z, мм')
else 
    xlabel('r, мм'), ylabel('x, мм')
end

rectangle('Position', [-W/2 -LL W 2*LL])
%title("Намотка 2D")
axis("equal")
grid("on")
hold("off")

% "Намотка 3D"
figure('Name', "Намотка 3D" ,'NumberTitle','off');
if axial  
    if ~isempty(Red1)
        Red1 = rotY(Red1, pi/2);
    end
    if ~isempty(Blue1)
        Blue1 = rotY(Blue1, pi/2);
    end
end
if ~isempty(Red1)
    plot3( Red1(:,1)*1000, Red1(:,2)*1000, Red1(:,3)*1000, '.r','MarkerSize',1)
    hold("on")
end
if ~isempty(Blue1)
    plot3(Blue1(:,1)*1000, Blue1(:,2)*1000, Blue1(:,3)*1000, '.b','MarkerSize',1)
end
box("on")
xlabel('x, мм'), ylabel('y, мм'), zlabel('z, мм')
%title("Намотка 3D")
axis("equal")
grid("on")
hold("off")

end
     
