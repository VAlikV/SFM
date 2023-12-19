% Две катушки (без момента сил)
clc
clear all
close all

%% ------------------------------------------------------------------------- % Задание параметров

tic
Calc = 1;
Interp = 1;
OutData = 0;

R1 = 0.055; % Радиус 1ой катушки
R2 = 0.085; % Радиус 2ой катушки
%nfi = 25; % Количество точек окружности
nfi = 25;
%nfi = 15;

L1 = 0.15; % Половина длины катушки (от 0 в одну сторону)
L2 = 0.15; % Половина длины катушки (от 0 в одну сторону)
%nL = 15;  % Количество точек по длинне
nL = 21;
%nL = 11;

nx = 7;    % Кол-во точек вдоль Х
ny = 7;    % Кол-во точек вдоль Y
nz = 7;    % Кол-во точек вдоль Z
lx = 0.066;  % Длина области вдоль Х
ly = 0.066;  % Длина области вдоль Y
lz = 0.066;  % Длина области вдоль Z

level = 0.0;  % Положение плоскости
plane = 0;  % Выбор плоскости 0)XY, 1)XZ, 2)YZ

%B0 = 10^(0); % Внешнее поле
alpfa = 50;  % Вклад индукции
betta = 0.001; % Вклад мощности

ro = 0.017; % Удельное сопротивление (Ом*мм^2/м)
t = 2;      % Тольщина слоя (мм) 

Nc1 = 20;        % Количество витков 1

er1 = 0.005;    % Отклонение от уровня 1
eb1 = 0.005;    % Отклонение от уровня 1 

Nc2 = 20;        % Количество витков 2

er2 = 0.005;    % Отклонение от уровня 2
eb2 = 0.005;    % Отклонение от уровня 2

%% ------------------------------------------------------------------------- % Создание массива узлов
% Первая катушка
N1 = nfi*nL; % Количество узлов

Node1 = zeros(nfi,3); % Сетка
Temp = zeros(nfi,3); 

f = linspace(-pi+2*pi/nfi,pi,nfi); % Массив значений угла
l = linspace(-L1,L1,nL); % Массив значений длины

Node1(:,1) = l(1);      % X
Node1(:,2) = R1*sin(f);  % Y
Node1(:,3) = R1*cos(f);  % Z

for i = 2:nL
    Temp(:,1) = l(i); % X
    Temp(:,2) = R1*sin(f); % Y        
    Temp(:,3) = R1*cos(f); % Z
    Node1 = [Node1; Temp];
end
clear f l Temp

% Вторая катушка

N2 = nfi*nL; % Количество узлов

Node2 = zeros(nfi,3); % Сетка
Temp = zeros(nfi,3); 

f = linspace(-pi+2*pi/nfi,pi,nfi); % Массив значений угла
l = linspace(-L2,L2,nL); % Массив значений длины

Node2(:,1) = l(1);      % X
Node2(:,2) = R2*sin(f);  % Y
Node2(:,3) = R2*cos(f);  % Z

for i = 2:nL
    Temp(:,1) = l(i); % X
    Temp(:,2) = R2*sin(f); % Y        
    Temp(:,3) = R2*cos(f); % Z
    Node2 = [Node2; Temp];
end
clear f l Temp

N = N1+N2;
Node = [Node1; Node2];

%% ------------------------------------------------------------------------- Создание массива треугольных элементов

% Первая катушка

S1 = zeros(N1,1);             % Массив элементов (хранится информация о номерах узлов, составляющих элемент)

i = 1;
while i <= N1-nfi            % Заолнение массива элементов     
    if mod(i,nfi) ~= 0 
        S1(2*i-1,1) = i;     % По два треугольника в одном квадрате
        S1(2*i-1,2) = i+nfi;
        S1(2*i-1,3) = i+nfi+1;
        S1(2*i,1) = i;
        S1(2*i,2) = i+nfi+1;
        S1(2*i,3) = i+1;
    else                    % Последние элементы в кольце
        S1(2*i-1,1) = i;     % По два треугольника в одном квадрате
        S1(2*i-1,2) = i+nfi;
        S1(2*i-1,3) = i+1;
        S1(2*i,1) = i;
        S1(2*i,2) = i+1;
        S1(2*i,3) = i-nfi+1;
    end
    i = i + 1;
end

T1 = size(S1,1); % Длинна массива элементов

% Вторая катушка

S2 = zeros(N1,1);             % Массив элементов (хранится информация о номерах узлов, составляющих элемент)

i = 1;
while i <= N2-nfi            % Заолнение массива элементов     
    if mod(i,nfi) ~= 0 
        S2(2*i-1,1) = i;     % По два треугольника в одном квадрате
        S2(2*i-1,2) = i+nfi;
        S2(2*i-1,3) = i+nfi+1;
        S2(2*i,1) = i;
        S2(2*i,2) = i+nfi+1;
        S2(2*i,3) = i+1;
    else                    % Последние элементы в кольце
        S2(2*i-1,1) = i;     % По два треугольника в одном квадрате
        S2(2*i-1,2) = i+nfi;
        S2(2*i-1,3) = i+1;
        S2(2*i,1) = i;
        S2(2*i,2) = i+1;
        S2(2*i,3) = i-nfi+1;
    end
    i = i + 1;
end

T2 = size(S2,1); % Длинна массива элементов
S2 = S2 + N2;

T = T1 + T2;
S = [S1; S2];

%% ------------------------------------------------------------------------- % Создание массива ROI (FOV)
 
ROI = [];

K = nx*ny*nz;  % Полное число точек ROI
dy = ly/(ny-1);
dz = lz/(nz-1);
for j = 1:nz
    z = -lz/2+(j-1)*dz;
    for i = 1:ny
        y = -ly/2+(i-1)*dy;
        Temp = [transpose(linspace(-lx/2,lx/2,nx)) y*ones(nx,1) z*ones(nx,1)];
        ROI = [ROI; Temp];
    end
end

clear x y z dx dy dz Temp

%% ------------------------------------------------------------------------- % Отрисовка сетки

n = 25;
Si = N2T(n,S); % Проверка функции N2T

%Ar = Area(S(St,:), Node);
%ri = N2Tri(n,S, nfi, N);

figure('Name','Сетка','NumberTitle','off'); 
movegui([0 560]);
plot3(Node(:,1),Node(:,2),Node(:,3),'b.'); % Отрисовка узлов
hold on
xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]')
axis equal
grid on 
plot3(Node(S(Si,:),1),Node(S(Si,:),2),Node(S(Si,:),3),'g*'); % Отрисовка элемента (N2T) 
% for i=1:N
%     text(Node(i,1),Node(i,2),Node(i,3),num2str(i),'Color','black','FontSize',8) % Нумерация узлов
% end
for i = 1:T
    plot3(Node(S(i,:),1),Node(S(i,:),2),Node(S(i,:),3),'r-'); % Отрисовка сетки
end
plot3(ROI(:,1),ROI(:,2),ROI(:,3),'k.'); % Отрисовка ROI

% for St = 1:length(Si)
%     Jn = CalcJn(n, S(Si(St),:), Node);
%     quiver3(Jn(2,1),Jn(2,2),Jn(2,3),Jn(1,1),Jn(1,2),Jn(1,3),'Color','b','LineWidth',1,'ShowArrowHead', 'on')
% end

%% ------------------------------------------------------------------------- % Задание целевого поля

Bdes = ROI(:,1) + ROI(:,2);%abs(ROI(:,1)).*ROI(:,1); % ==================================================================================================== Вычисление целевого поля

[X,Y,Z,B] = Transform3D(lx, ly, lz, nx, ny, nz, Bdes); % Для отрисовки поля

figure('Name','Целевое поле','NumberTitle','off'); 
movegui([570 560]);
if (plane == 0) % ---------------------------------------------------------- XY
    k = find(Z(1,1,:)==level);
    surf(X(:,:,k), Y(:,:,k), B(:,:,k))
    xlabel ('x [m]'), ylabel ('y [m]'), title(strcat('Bz [T], z = ', num2str(level))) 
elseif (plane == 1) % ------------------------------------------------------ XZ 
    k = find(Y(1,:,1)==level);
    XX = reshape(X(:,k,:),[nx,nz,1]);
    ZZ = reshape(Z(:,k,:),[nx,nz,1]);
    BB = reshape(B(:,k,:),[nx,nz,1]);
    surf(XX, ZZ, BB)
    xlabel ('x [m]'), ylabel ('z [m]'), title(strcat('Bz [T], y = ', num2str(level)))
elseif (plane == 2) % ------------------------------------------------------ YZ
    k = find(X(:,1,1)==level);
    YY = reshape(Y(k,:,:),[ny,nz,1]);
    ZZ = reshape(Z(k,:,:),[ny,nz,1]);
    BB = reshape(B(k,:,:),[ny,nz,1]);
    surf(YY, ZZ, BB)
    xlabel ('y [m]'), ylabel ('z [m]'), title(strcat('Bz [T], x = ', num2str(level)))
else
    disp('___ОШИБКА ВЫБОРА ПЛОСКОСТИ___');
end
colorbar
box on
hold on
grid on 


%% ------------------------------------------------------------------------- % Задание матрицы А

A = zeros(N,N);

if Calc
   parfor j = 1:N
        for i = 1:N
            bz1 = Calcbz(i, S, Node, ROI);
            bz2 = Calcbz(j, S, Node, ROI);
            Lmn = CalcLmn(i, j, S, Node);
            pmn = CalcPmn(i, j, ro, t, S, Node);
            Arr = bz1.*bz2;
            A(j,i) = zum(Arr) + alpfa*Lmn + betta*pmn;
        end
    end 
end

clear Temp bz1 bz2

%% ------------------------------------------------------------------------- % Задание матрицы B

B = zeros(N,1);
    
if Calc
    for i = 1:N
        bz1 = Calcbz(i, S, Node, ROI);
        Arr = Bdes(:,1).*bz1;
        B(i,1) = zum(Arr);
    end 
end 

clear Temp i k j

%% ------------------------------------------------------------------------- % Задание матрицы I

I = pinv(A)*B;

if OutData && ~Calc
    
    % Тесты
    I = [-0.105886074097665;-0.109207658610500;-0.114081049827258;-0.119654620024240;-0.124942373371644;-0.129007261357578;-0.131171087781932;-0.131092280104974;-0.128770808079363;-0.124562376279448;-0.119183223527272;-0.113606924000769;-0.108833289457430;-0.105683508477224;-0.104652102604198;-0.0956596916052332;-0.101697162913220;-0.110648097840511;-0.120906434280707;-0.130635854396229;-0.138082047021806;-0.142041930823492;-0.141943672275944;-0.137777252095138;-0.130126663785451;-0.120267059262040;-0.110016333252128;-0.101230876617420;-0.0954277989834693;-0.0934836197218663;-0.0782417124993309;-0.0887194733188487;-0.104707942510022;-0.123048203260148;-0.140454353560107;-0.153497130313997;-0.160309529691129;-0.160236144598969;-0.153194690730624;-0.139787387682477;-0.122111357852467;-0.103761629831250;-0.0881626982629934;-0.0780727451961741;-0.0746924509239331;-0.0545799176527920;-0.0702103836079819;-0.0961569582944817;-0.125299561231733;-0.153657752223458;-0.173723754656374;-0.183675051437902;-0.183925577885437;-0.174218936053524;-0.153719907329693;-0.124956266963022;-0.0955159902423346;-0.0705983091675963;-0.0554359143823082;-0.0504912120306182;-0.0775244002415974;-0.0840249766483679;-0.0996743717502087;-0.114607562907232;-0.132163093580628;-0.150950596558607;-0.160716549568221;-0.159812500500165;-0.155700105814864;-0.144012024307752;-0.126607582594910;-0.112978488156446;-0.0926160012324548;-0.0777089765509541;-0.0748317194292793;-0.130170453861140;-0.130838388971388;-0.139991606718739;-0.108620978452117;-0.0910124354286091;-0.109400567573046;-0.0962562492394419;-0.0962562492394420;-0.109400567573046;-0.0910124354286090;-0.108620978452117;-0.139991606718738;-0.130838388971388;-0.130170453861140;-0.140360848691529;-0.0777089765509545;-0.0926160012324549;-0.112978488156447;-0.126607582594911;-0.144012024307752;-0.155700105814864;-0.159812500500165;-0.160716549568222;-0.150950596558607;-0.132163093580628;-0.114607562907233;-0.0996743717502091;-0.0840249766483680;-0.0775244002415975;-0.0748317194292796;-0.0554359143823083;-0.0705983091675966;-0.0955159902423348;-0.124956266963022;-0.153719907329694;-0.174218936053524;-0.183925577885437;-0.183675051437902;-0.173723754656375;-0.153657752223459;-0.125299561231734;-0.0961569582944820;-0.0702103836079821;-0.0545799176527920;-0.0504912120306183;-0.0780727451961743;-0.0881626982629937;-0.103761629831250;-0.122111357852467;-0.139787387682477;-0.153194690730624;-0.160236144598969;-0.160309529691130;-0.153497130313998;-0.140454353560107;-0.123048203260148;-0.104707942510023;-0.0887194733188491;-0.0782417124993312;-0.0746924509239334;-0.0954277989834696;-0.101230876617420;-0.110016333252129;-0.120267059262040;-0.130126663785452;-0.137777252095138;-0.141943672275945;-0.142041930823492;-0.138082047021806;-0.130635854396229;-0.120906434280707;-0.110648097840512;-0.101697162913221;-0.0956596916052336;-0.0934836197218667;-0.105683508477225;-0.108833289457431;-0.113606924000770;-0.119183223527273;-0.124562376279448;-0.128770808079363;-0.131092280104974;-0.131171087781932;-0.129007261357579;-0.124942373371645;-0.119654620024241;-0.114081049827259;-0.109207658610501;-0.105886074097665;-0.104652102604198;0.0893987186601965;0.0862456271645018;0.0816467304556641;0.0764312052568382;0.0715151757257921;0.0677571149795416;0.0657858342227284;0.0659000124798106;0.0680835476661810;0.0720063228513800;0.0769987312427955;0.0821739368682659;0.0866273644319025;0.0895868287013571;0.0905649171772857;0.0934672915521701;0.0892345004588302;0.0830182614653844;0.0759596871136453;0.0692957652363916;0.0641969813928404;0.0615137991583356;0.0616302072002097;0.0645289411852652;0.0697988881215495;0.0765418626776325;0.0835481441763445;0.0895946081523346;0.0936256307523881;0.0949892987538333;0.0983599227701986;0.0928425292914974;0.0846673900790169;0.0754297568175956;0.0666798620773685;0.0599862399957675;0.0564879285084952;0.0565962039377324;0.0602847841399762;0.0671428404920876;0.0759743802769008;0.0851453581324477;0.0931043171097982;0.0984082142117727;0.100239238773916;0.101483827009139;0.0951228972549156;0.0857165686933622;0.0753738017724916;0.0653738517184774;0.0575667904372371;0.0535458133625621;0.0536106722091800;0.0576732114407371;0.0654981998079629;0.0754979191216963;0.0857939272342936;0.0950354403503187;0.101248020064123;0.103419562668464;0.0974911626041478;0.0916700819271100;0.0839380329643230;0.0766398481043163;0.0687460545817883;0.0615982983902783;0.0580130468900620;0.0579802009224210;0.0613739673380740;0.0682588812678865;0.0758546444084173;0.0832237392140224;0.0913144322825820;0.0970869373761375;0.0989492844110434;0.0922490798213434;0.0871649712738974;0.0809460139698923;0.0771981645370776;0.0722195955565074;0.0662063206439295;0.0635014957187669;0.0635014957187672;0.0662063206439302;0.0722195955565081;0.0771981645370781;0.0809460139698926;0.0871649712738974;0.0922490798213432;0.0933195131130446;0.0970869373761375;0.0913144322825818;0.0832237392140221;0.0758546444084170;0.0682588812678860;0.0613739673380735;0.0579802009224209;0.0580130468900625;0.0615982983902791;0.0687460545817891;0.0766398481043168;0.0839380329643233;0.0916700819271102;0.0974911626041479;0.0989492844110433;0.101248020064123;0.0950354403503185;0.0857939272342933;0.0754979191216960;0.0654981998079626;0.0576732114407369;0.0536106722091799;0.0535458133625625;0.0575667904372377;0.0653738517184780;0.0753738017724921;0.0857165686933625;0.0951228972549157;0.101483827009139;0.103419562668464;0.0984082142117725;0.0931043171097979;0.0851453581324474;0.0759743802769005;0.0671428404920873;0.0602847841399761;0.0565962039377324;0.0564879285084955;0.0599862399957681;0.0666798620773690;0.0754297568175961;0.0846673900790173;0.0928425292914977;0.0983599227701987;0.100239238773916;0.0936256307523880;0.0895946081523343;0.0835481441763442;0.0765418626776323;0.0697988881215494;0.0645289411852652;0.0616302072002099;0.0615137991583360;0.0641969813928409;0.0692957652363922;0.0759596871136458;0.0830182614653849;0.0892345004588305;0.0934672915521704;0.0949892987538334;0.0895868287013570;0.0866273644319023;0.0821739368682656;0.0769987312427953;0.0720063228513799;0.0680835476661811;0.0659000124798108;0.0657858342227287;0.0677571149795421;0.0715151757257927;0.0764312052568388;0.0816467304556645;0.0862456271645022;0.0893987186601967;0.0905649171772857];

end

I1 = I(1:N1);
I2 = I(N1+1:N);

[Fi1, Ll1, II1] = Transform2D(L1, nfi, nL, I1);
[Fi2, Ll2, II2] = Transform2D(L2, nfi, nL, I2);

%% ------------------------------------------------------------------------- % Отрисовка потоковой функции первой катушки

rgb = ones(T1,3);
if Calc || OutData
    figure('Name','Потоковая функция первой катушки','NumberTitle','off'); 
    movegui([1145 560]);
    contourf(Fi1,Ll1,II1);
    %surf(fi,XX,II);
    box on
    hold on
    xlabel('fi'), ylabel('L'), title('Stream Function')
    colorbar
    
    StFun = zeros(T1,1);

    for i = 1:T1
        StFun(i,1) = (I1(S1(i,1),1) + I1(S1(i,2),1) + I1(S1(i,3),1))/3; 
    end
    
    StFun = StFun - min(min(StFun));
    StFun = 0.6667*StFun/max(max(StFun));
    
    for i = 1:T1
        rgb(i,:) = hsv2rgb([StFun(i,1) 1 1]);
    end
end

%% ------------------------------------------------------------------------- % Отрисовка элементов первой катушки

figure('Name','Элементы первой катушки','NumberTitle','off'); 
movegui([0 30]);
x = [Node(S1(1,1),1) Node(S1(1,2),1) Node(S1(1,3),1)];
y = [Node(S1(1,1),2) Node(S1(1,2),2) Node(S1(1,3),2)];
z = [Node(S1(1,1),3) Node(S1(1,2),3) Node(S1(1,3),3)];
fill3(x, y, z,rgb(1,:))
hold on
%text(sum(x)/length(x),sum(y)/length(y),sum(z)/length(z),num2str(1),'Color','blue','FontSize',8) % Нумерация 1ого элемента
xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]')
axis equal
grid on
for i = 2:T1
    x = [Node(S1(i,1),1) Node(S1(i,2),1) Node(S1(i,3),1)];
    y = [Node(S1(i,1),2) Node(S1(i,2),2) Node(S1(i,3),2)];
    z = [Node(S1(i,1),3) Node(S1(i,2),3) Node(S1(i,3),3)];
    fill3(x, y, z,rgb(i,:)) % Отрисовка элементов
    %text(sum(x)/length(x),sum(y)/length(y),sum(z)/length(z),num2str(i),'Color','blue','FontSize',8) % Нумерация элементов
end

%% ------------------------------------------------------------------------- % Отрисовка потоковой функции второй катушки

rgb2 = ones(T2,3);
if Calc || OutData
    figure('Name','Потоковая функция второй катушки','NumberTitle','off'); 
    movegui([1145 560]);
    contourf(Fi2,Ll2,II2);
    %surf(fi,XX,II);
    box on
    hold on
    xlabel('fi'), ylabel('L'), title('Stream Function')
    colorbar
    
    StFun = zeros(T2,1);

    for i = 1:T2
        StFun(i,1) = (I2(S2(i,1)-N1,1) + I2(S2(i,2)-N1,1) + I2(S2(i,3)-N1,1))/3; 
    end
    
    StFun = StFun - min(min(StFun));
    StFun = 0.6667*StFun/max(max(StFun));
    
    for i = 1:T2
        rgb2(i,:) = hsv2rgb([StFun(i,1) 1 1]);
    end
end

%% ------------------------------------------------------------------------- % Отрисовка элементов второй катушки

figure('Name','Элементы второй катушки','NumberTitle','off'); 
movegui([0 30]);
x = [Node(S2(1,1),1) Node(S2(1,2),1) Node(S2(1,3),1)];
y = [Node(S2(1,1),2) Node(S2(1,2),2) Node(S2(1,3),2)];
z = [Node(S2(1,1),3) Node(S2(1,2),3) Node(S2(1,3),3)];
fill3(x, y, z,rgb2(1,:))
hold on
%text(sum(x)/length(x),sum(y)/length(y),sum(z)/length(z),num2str(1),'Color','blue','FontSize',8) % Нумерация 1ого элемента
xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]')
axis equal
grid on
for i = 2:T2
    x = [Node(S2(i,1),1) Node(S2(i,2),1) Node(S2(i,3),1)];
    y = [Node(S2(i,1),2) Node(S2(i,2),2) Node(S2(i,3),2)];
    z = [Node(S2(i,1),3) Node(S2(i,2),3) Node(S2(i,3),3)];
    fill3(x, y, z,rgb2(i,:)) % Отрисовка элементов
    %text(sum(x)/length(x),sum(y)/length(y),sum(z)/length(z),num2str(i),'Color','blue','FontSize',8) % Нумерация элементов
end

%% ------------------------------------------------------------------------- % Поле
 
if Calc || OutData
    
    Bz = zeros(K,1);
    Temp = 0;
    for j = 1:N
        b = Calcbz(j, S, Node, ROI);
        Temp = Temp + I(j,1)*b;
    end
    Bz(:,1) = Temp;

    [X1,Y1,Z1,Bz] = Transform3D(lx, ly, lz, nx, ny, nz, Bz);

    figure('Name','Полученное поле','NumberTitle','off'); 
    movegui([570 30]);
    if (plane == 0) % ---------------------------------------------------------- XY
        k = find(Z1(1,1,:)==level);
        surf(X1(:,:,k), Y1(:,:,k), Bz(:,:,k))
        xlabel ('x [m]'), ylabel ('y [m]'), title(strcat('Bz [T], z = ', num2str(level))) 
    elseif (plane == 1) % ------------------------------------------------------ XZ
        k = find(Y1(1,:,1)==level);
        XX1 = reshape(X1(:,k,:),[nx,nz,1]);
        ZZ1 = reshape(Z1(:,k,:),[nx,nz,1]);
        BB1 = reshape(Bz(:,k,:),[nx,nz,1]);
        surf(XX1, ZZ1, BB1)
        xlabel ('x [m]'), ylabel ('z [m]'), title(strcat('Bz [T], y = ', num2str(level)))
    elseif (plane == 2) % ------------------------------------------------------ YZ
        k = find(X1(:,1,1)==level);
        YY1 = reshape(Y1(k,:,:),[ny,nz,1]);
        ZZ1 = reshape(Z1(k,:,:),[ny,nz,1]);
        BB1 = reshape(Bz(k,:,:),[ny,nz,1]);
        surf(YY1, ZZ1, BB1)
        xlabel ('y [m]'), ylabel ('z [m]'), title(strcat('Bz [T], x = ', num2str(level)))
    else
        disp('___ОШИБКА ВЫБОРА ПЛОСКОСТИ___');
    end
    colorbar
    box on
    hold on
    grid on
end

%% ------------------------------------------------------------------------- % Интерполяция

if (Calc && Interp) || (OutData && Interp)
    nfi_int = 2001;
    nL_int = 1001;

    [Fi_int1, L_int1] = meshgrid(linspace(-pi,pi,nfi_int), linspace(-L1,L1,nL_int));

    I_int1 = interp2(Fi1,Ll1,II1,Fi_int1,L_int1);
end

if (Calc && Interp) || (OutData && Interp)

    [Fi_int2, L_int2] = meshgrid(linspace(-pi,pi,nfi_int), linspace(-L2,L2,nL_int));

    I_int2 = interp2(Fi2,Ll2,II2,Fi_int2,L_int2);
end

%% ------------------------------------------------------------------------- % Контур тока первой катушки

if (Calc && Interp) || (OutData && Interp)
    
    Max = max(max(I_int1));  % Максимальное значение функции потока
    Min = min(min(I_int1));  % Минимальное значение функции потока
    df = (Max-Min)/(Nc1+1);  % Шаг витков
    
    Redf = [];               % Контура с током против часовой
    Bluef = [];              % Контура с током по часовой 

    for k = 1:Nc1
        K1 = Min+df*k;       % Уровень
        
        if K1 >= (Max+Min)/2
            AA = find((I_int1 >= ((Min+df*k)-er1*(Min+df*k)) & I_int1 <= ((Min+df*k)+er1*(Min+df*k)) & I_int1 >= 0) | (I_int1 >= ((Min+df*k)+er1*(Min+df*k)) & I_int1 <= ((Min+df*k)-er1*(Min+df*k)) & I_int1 < 0));
        else
            AA = find((I_int1 >= ((Min+df*k)-eb1*(Min+df*k)) & I_int1 <= ((Min+df*k)+eb1*(Min+df*k)) & I_int1 >= 0) | (I_int1 >= ((Min+df*k)+eb1*(Min+df*k)) & I_int1 <= ((Min+df*k)-eb1*(Min+df*k)) & I_int1 < 0));
        end
       
        i = ceil(AA/(nL_int));   % Номера столбцов

        j = mod(AA, nL_int);     % Номера строк
        j(j==0) = nL_int;

        Temp = zeros(length(i),2);

        for a = 1:length(i)
            Temp(a,:) = [L_int1(j(a),i(a)) Fi_int1(j(a),i(a))];
        end

        if K1 >= (Max+Min)/2
            Redf = [Redf; Temp];
        else
            Bluef = [Bluef; Temp];
        end
    end
    
    Lred1 = length(Redf);
    Lblue1 = length(Bluef);
    
    if Lred1 ~= 0
        Red1f = zeros(length(Redf(:,1)),3);   % Массивы для записи контура в координатах
        Red2f = zeros(length(Redf(:,1)),2);   % Массивы для записи контура в плоскости
        Red1f(:,1) = Redf(:,1);               % X
        Red1f(:,2) = R1*sin(Redf(:,2));        % Y
        Red1f(:,3) = R1*cos(Redf(:,2));        % Z
        
        Red2f(:,1) = Redf(:,1);               % X
        Red2f(:,2) = Redf(:,2)*R1;             % R
    end
    
    if Lblue1 ~= 0
        Blue1f = zeros(length(Bluef(:,1)),3); % Массивы для записи контура в координатах 
        Blue2f = zeros(length(Bluef(:,1)),2); % Массивы для записи контура в плоскости
        Blue1f(:,1) = Bluef(:,1);             % X
        Blue1f(:,2) = R1*sin(Bluef(:,2));      % Y
        Blue1f(:,3) = R1*cos(Bluef(:,2));      % Z
        
        Blue2f(:,1) = Bluef(:,1);             % X
        Blue2f(:,2) = R1*Bluef(:,2);           % R
    end
     
end

%% ------------------------------------------------------------------------- % Контур тока второй катушки

if (Calc && Interp) || (OutData && Interp)
    
    Max = max(max(I_int2));  % Максимальное значение функции потока
    Min = min(min(I_int2));  % Минимальное значение функции потока
    df = (Max-Min)/(Nc2+1);  % Шаг витков
    
    Reds = [];               % Контура с током против часовой
    Blues = [];              % Контура с током по часовой 

    for k = 1:Nc2
        K1 = Min+df*k;       % Уровень
        
        if K1 >= (Max+Min)/2
            AA = find((I_int2 >= ((Min+df*k)-er2*(Min+df*k)) & I_int2 <= ((Min+df*k)+er2*(Min+df*k)) & I_int2 >= 0) | (I_int2 >= ((Min+df*k)+er2*(Min+df*k)) & I_int2 <= ((Min+df*k)-er2*(Min+df*k)) & I_int2 < 0));
        else
            AA = find((I_int2 >= ((Min+df*k)-eb2*(Min+df*k)) & I_int2 <= ((Min+df*k)+eb2*(Min+df*k)) & I_int2 >= 0) | (I_int2 >= ((Min+df*k)+eb2*(Min+df*k)) & I_int2 <= ((Min+df*k)-eb2*(Min+df*k)) & I_int2 < 0));
        end
       
        i = ceil(AA/(nL_int));   % Номера столбцов

        j = mod(AA, nL_int);     % Номера строк
        j(j==0) = nL_int;

        Temp = zeros(length(i),2);

        for a = 1:length(i)
            Temp(a,:) = [L_int2(j(a),i(a)) Fi_int2(j(a),i(a))];
        end

        if K1 >= (Max+Min)/2
            Reds = [Reds; Temp];
        else
            Blues = [Blues; Temp];
        end
    end
    
    Lred2 = length(Reds);
    Lblue2 = length(Blues);
    
    if Lred2 ~= 0
        Red1s = zeros(length(Reds(:,1)),3);   % Массивы для записи контура в координатах
        Red2s = zeros(length(Reds(:,1)),2);   % Массивы для записи контура в плоскости
        Red1s(:,1) = Reds(:,1);               % X
        Red1s(:,2) = R2*sin(Reds(:,2));        % Y
        Red1s(:,3) = R2*cos(Reds(:,2));        % Z
        
        Red2s(:,1) = Reds(:,1);               % X
        Red2s(:,2) = Reds(:,2)*R2;             % R
    end
    
    if Lblue2 ~= 0
        Blue1s = zeros(length(Blues(:,1)),3); % Массивы для записи контура в координатах 
        Blue2s = zeros(length(Blues(:,1)),2); % Массивы для записи контура в плоскости
        Blue1s(:,1) = Blues(:,1);             % X
        Blue1s(:,2) = R2*sin(Blues(:,2));      % Y
        Blue1s(:,3) = R2*cos(Blues(:,2));      % Z
        
        Blue2s(:,1) = Blues(:,1);             % X
        Blue2s(:,2) = R2*Blues(:,2);           % R
    end
     
end

%% ------------------------------------------------------------------------- % Отрисовка интерполированной фенкции и контуров

if (Calc && Interp) || (OutData && Interp)
    figure('Name','Интерполированная потоковая функция первой катушки','NumberTitle','off');    % 1ая Интерполированная функция 
    movegui([1255 560]);
    contourf(Fi_int1,L_int1,I_int1)
    box on
    hold on
    xlabel('fi'), ylabel('L'), title('Stream Function')
    if Lred1 ~= 0
        plot(Redf(:,2), Redf(:,1), '.r')
    end
    if Lblue1 ~= 0
        plot(Bluef(:,2), Bluef(:,1), '.b')
    end  
    colorbar

    figure('Name','Намотка 3D первой катушки','NumberTitle','off');    % Намотка 3D 
    movegui([1145 30]);
    if Lred1 ~= 0
        plot3(Red1f(:,1), Red1f(:,2), Red1f(:,3), '.r')
        hold on
    end
    if Lblue1 ~= 0
        plot3(Blue1f(:,1), Blue1f(:,2), Blue1f(:,3), '.b')
    end
    box on
    xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]')
    axis equal
    grid on 
    
    figure('Name','Намотка 2D первой катушки','NumberTitle','off');     % Намотка 2D 
    movegui([1255 30]);
    if Lred1 ~= 0
        plot(Red2f(:,2), Red2f(:,1), '.r')
        hold on
    end
    if Lblue1 ~= 0
        plot(Blue2f(:,2), Blue2f(:,1), '.b')
    end
    box on
    xlabel ('r [m]'), ylabel ('x [m]')
    axis equal
    grid on 
end

if (Calc && Interp) || (OutData && Interp)
    figure('Name','Интерполированная потоковая функция второй катушки','NumberTitle','off');    % 2ая Интерполированная функция 
    movegui([1255 560]);
    contourf(Fi_int2,L_int2,I_int2)
    box on
    hold on
    xlabel('fi'), ylabel('L'), title('Stream Function')
    if Lred2 ~= 0
        plot(Reds(:,2), Reds(:,1), '.r')
    end
    if Lblue2 ~= 0
        plot(Blues(:,2), Blues(:,1), '.b')
    end  
    colorbar

    figure('Name','Намотка 3D второй катушки','NumberTitle','off');    % Намотка 3D 
    movegui([1145 30]);
    if Lred2 ~= 0
        plot3(Red1s(:,1), Red1s(:,2), Red1s(:,3), '.r')
        hold on
    end
    if Lblue2 ~= 0
        plot3(Blue1s(:,1), Blue1s(:,2), Blue1s(:,3), '.b')
    end
    box on
    xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]')
    axis equal
    grid on 
    
    figure('Name','Намотка 2D второй катушки','NumberTitle','off');     % Намотка 2D 
    movegui([1255 30]);
    if Lred2 ~= 0
        plot(Red2s(:,2), Red2s(:,1), '.r')
        hold on
    end
    if Lblue2 ~= 0
        plot(Blue2s(:,2), Blue2s(:,1), '.b')
    end
    box on
    xlabel ('r [m]'), ylabel ('x [m]')
    axis equal
    grid on 
end
toc

%% ------------------------------------------------------------------------- % Функции

function Lmn = CalcLmn(n, m, S, Node) % Параметр Lmn 
    Sin = N2T(n, S);
    Sim = N2T(m, S);
    Sum = 0;
    for i=1:length(Sim)
        for j=1:length(Sin)
            Arn = Area(S(Sin(j),:),Node);
            Arm = Area(S(Sim(i),:),Node);
            Jn = CalcJn(n, S(Sin(j),:), Node);
            Jm = CalcJn(m, S(Sim(i),:), Node);
            rn = [(Node(S(Sin(j),1),1)+Node(S(Sin(j),2),1)+Node(S(Sin(j),3),1))/3 ...
                (Node(S(Sin(j),1),2)+Node(S(Sin(j),2),2)+Node(S(Sin(j),3),2))/3 ...
                (Node(S(Sin(j),1),3)+Node(S(Sin(j),2),3)+Node(S(Sin(j),3),3))/3];
            
            rm = [(Node(S(Sim(i),1),1)+Node(S(Sim(i),2),1)+Node(S(Sim(i),3),1))/3 ...
                (Node(S(Sim(i),1),2)+Node(S(Sim(i),2),2)+Node(S(Sim(i),3),2))/3 ...
                (Node(S(Sim(i),1),3)+Node(S(Sim(i),2),3)+Node(S(Sim(i),3),3))/3];
            d = Distance(rn, rm);
            if Sim(i)~=Sin(j)
                Sum = Sum + ScalarP(Jn(1,:),Jm(1,:))*Arn*Arm/d;
            else
                r1 = Node(S(Sim(i),1),:);
                r2 = Node(S(Sim(i),2),:);
                r3 = Node(S(Sim(i),3),:);
                aa = ScalarP((r3-r1),(r3-r1));
                bb = ScalarP((r3-r1),(r3-r2));
                cc = ScalarP((r3-r2),(r3-r2));
                p1 = sqrt(aa*cc);
                p2 = sqrt(aa-2*bb+cc);
                Sum = Sum + (4*ScalarP(Jn(1,:),Jm(1,:))*Arn*Arm)*((1/(6*sqrt(aa)))*log(((aa-bb+sqrt(aa)*p2)*(bb+p1))/((-bb+p1)*(-aa+bb+sqrt(aa)*p2)))+...
                    (1/(6*sqrt(cc)))*log(((bb+p1)*(-bb+cc+sqrt(cc)*p2))/((bb-cc+sqrt(cc)*p2)*(-bb+p1)))+...
                    (1/(6*p2))*log(((aa-bb+sqrt(aa)*p2)*(-bb+cc+sqrt(cc)*p2))/((bb-cc+sqrt(cc)*p2)*(-aa+bb+sqrt(aa)*p2))));
            end
        end
    end
    Lmn = Sum*10^(-7);
end

function bz = Calcbz(n, S, Node, ROI) % Параметр b 
    Si = N2T(n, S); % 
    Sum = zeros(length(ROI(:,1)),1);
    r1 = zeros(1,3); % Радиус вектор узлов
    r2 = zeros(1,3);
    r3 = zeros(1,3);
    point = zeros(4,3);
    w = [-9/32 25/96 25/96 25/96];
    for i=1:length(Si)
        integral = zeros(length(ROI(:,1)),1);  
        Jn = CalcJn(n, S(Si(i),:), Node);
        r1(:) = [Node(S(Si(i),1),1) Node(S(Si(i),1),2) Node(S(Si(i),1),3)]; 
        r2(:) = [Node(S(Si(i),2),1) Node(S(Si(i),2),2) Node(S(Si(i),2),3)];  
        r3(:) = [Node(S(Si(i),3),1) Node(S(Si(i),3),2) Node(S(Si(i),3),3)];
        point(1,:) = (r1 + r2 - 2*r3)/3 + r3;
        point(2,:) = 3*(r1 - r3)/5 + (r2 - r3)/5 + r3;
        point(3,:) = (r1 - r3)/5 + 3*(r2 - r3)/5 + r3;
        point(4,:) = (r1 + r2 - 2*r3)/5 + r3;      
        for j = 1:4
            d = DistanceROI(point(j,:),ROI);
            integral = integral + (w(j)*(Jn(1,2)*(ROI(:,1)-point(j,1))-Jn(1,1)*(ROI(:,2)-point(j,2))))./(d.^3);
        end      
        Sum = Sum + integral;
    end
    bz = Sum*10^(-7);
end

function Mx = CalcMx(n, S, Node, B0) % Параметр Mx 
    Sum = 0;
    Si = N2T(n, S); % 
    r1 = zeros(1,3); % Радиус вектор узлов
    r2 = zeros(1,3);
    r3 = zeros(1,3);
    point = zeros(4,3);
    w = [-9/32 25/96 25/96 25/96];
    for i=1:length(Si)
        integral = 0; 
        Jn = CalcJn(n, S(Si(i),:), Node);
        r1(:) = [Node(S(Si(i),1),1) Node(S(Si(i),1),2) Node(S(Si(i),1),3)]; 
        r2(:) = [Node(S(Si(i),2),1) Node(S(Si(i),2),2) Node(S(Si(i),2),3)];  
        r3(:) = [Node(S(Si(i),3),1) Node(S(Si(i),3),2) Node(S(Si(i),3),3)];
        point(1,:) = (r1 + r2 - 2*r3)/3 + r3;
        point(2,:) = 3*(r1 - r3)/5 + (r2 - r3)/5 + r3;
        point(3,:) = (r1 - r3)/5 + 3*(r2 - r3)/5 + r3;
        point(4,:) = (r1 + r2 - 2*r3)/5 + r3;      
        for j = 1:4
            integral = integral + (w(j)*(Jn(1,1)*point(j,3)));
        end      
        Sum = Sum + integral;
    end
    Mx = Sum*B0;
end

function My = CalcMy(n, S, Node, B0) % Параметр My 
    Sum = 0;
    Si = N2T(n, S); % 
    r1 = zeros(1,3); % Радиус вектор узлов
    r2 = zeros(1,3);
    r3 = zeros(1,3);
    point = zeros(4,3);
    w = [-9/32 25/96 25/96 25/96];
    for i=1:length(Si)
        integral = 0; 
        Jn = CalcJn(n, S(Si(i),:), Node);
        r1(:) = [Node(S(Si(i),1),1) Node(S(Si(i),1),2) Node(S(Si(i),1),3)]; 
        r2(:) = [Node(S(Si(i),2),1) Node(S(Si(i),2),2) Node(S(Si(i),2),3)];  
        r3(:) = [Node(S(Si(i),3),1) Node(S(Si(i),3),2) Node(S(Si(i),3),3)];
        point(1,:) = (r1 + r2 - 2*r3)/3 + r3;
        point(2,:) = 3*(r1 - r3)/5 + (r2 - r3)/5 + r3;
        point(3,:) = (r1 - r3)/5 + 3*(r2 - r3)/5 + r3;
        point(4,:) = (r1 + r2 - 2*r3)/5 + r3;      
        for j = 1:4
            integral = integral + (w(j)*(Jn(1,2)*point(j,3)));
        end      
        Sum = Sum + integral;
    end
    My = Sum*B0;
end

function Mz = CalcMz(n, S, Node, B0) % Параметр Mz 
    Sum = 0;
    Si = N2T(n, S); % 
    r1 = zeros(1,3); % Радиус вектор узлов
    r2 = zeros(1,3);
    r3 = zeros(1,3);
    point = zeros(4,3);
    w = [-9/32 25/96 25/96 25/96];
    for i=1:length(Si)
        integral = 0; 
        Jn = CalcJn(n, S(Si(i),:), Node);
        r1(:) = [Node(S(Si(i),1),1) Node(S(Si(i),1),2) Node(S(Si(i),1),3)]; 
        r2(:) = [Node(S(Si(i),2),1) Node(S(Si(i),2),2) Node(S(Si(i),2),3)];  
        r3(:) = [Node(S(Si(i),3),1) Node(S(Si(i),3),2) Node(S(Si(i),3),3)];
        point(1,:) = (r1 + r2 - 2*r3)/3 + r3;
        point(2,:) = 3*(r1 - r3)/5 + (r2 - r3)/5 + r3;
        point(3,:) = (r1 - r3)/5 + 3*(r2 - r3)/5 + r3;
        point(4,:) = (r1 + r2 - 2*r3)/5 + r3;      
        for j = 1:4
            integral = integral - (w(j)*(Jn(1,1)*point(j,1) + Jn(1,2)*point(j,2)));
        end      
        Sum = Sum + integral;
    end
    Mz = Sum*B0;
end

function pmn = CalcPmn(n, m, ro, t, S, Node) % Параметр p 
    Sin = N2T(n, S);
    Sim = N2T(m, S);
    Sum = 0;
    for i=1:length(Sim)
        for j=1:length(Sin)
            Arn = Area(S(Sin(j),:),Node);
            Arm = Area(S(Sim(i),:),Node);
            Jn = CalcJn(n, S(Sin(j),:), Node);
            Jm = CalcJn(m, S(Sim(i),:), Node);
            if Sim(i)==Sin(j)
                Sum = Sum + ScalarP(Jn(1,:),Jm(1,:))*Arn;
            end
        end
    end
    pmn = Sum*(ro/t);
end

function Jn = CalcJn(n, St, Node) % Базисный вектор тока
    t = find(St == n);
    Jn = zeros(2,3);
    Ar = Area(St, Node);
    if t == 1
        Jn(1,:) = (Node(St(3),:)-Node(St(2),:))/(2*Ar);
        Jn(2,:) = Node(St(2),:);
        %Jn(2,:) = -(Node(St(3),:)-Node(St(1),:))/(2*Area);
        %Jn(3,:) = -(Node(St(1),:)-Node(St(2),:))/(2*Area);
    elseif t == 2
        Jn(1,:) = (Node(St(1),:)-Node(St(3),:))/(2*Ar);
        Jn(2,:) = Node(St(3),:);
        %Jn(2,:) = -(Node(St(1),:)-Node(St(2),:))/(2*Area);
        %Jn(3,:) = -(Node(St(2),:)-Node(St(3),:))/(2*Area);
    elseif t == 3
        Jn(1,:) = (Node(St(2),:)-Node(St(1),:))/(2*Ar);
        Jn(2,:) = Node(St(1),:);
        %Jn(2,:) = -(Node(St(2),:)-Node(St(3),:))/(2*Area);
        %Jn(3,:) = -(Node(St(3),:)-Node(St(1),:))/(2*Area);
    else
        disp('___ОШИБКА ВЫБОРА УЗЛА/ЭЛЕМЕНТА___');
    end
end

function Ar = Area(St, Node) % Возвращает площадь элемента 
    S1 = det([1 Node(St(1),1) Node(St(1), 2); 1 Node(St(2),1) Node(St(2), 2); 1 Node(St(3),1) Node(St(3), 2)]);
    S2 = det([1 Node(St(1),1) Node(St(1), 3); 1 Node(St(2),1) Node(St(2), 3); 1 Node(St(3),1) Node(St(3), 3)]);
    S3 = det([1 Node(St(1),2) Node(St(1), 3); 1 Node(St(2),2) Node(St(2), 3); 1 Node(St(3),2) Node(St(3), 3)]);
    Ar = 0.5*sqrt(S1^2+S2^2+S3^2);
end

function Ar1 = Area1(St, Node) % Возвращает площадь элемента 
    a = sqrt((Node(St(1),1)-Node(St(2),1))^2 + (Node(St(1),2)-Node(St(2),2))^2 + (Node(St(1),3)-Node(St(2),3))^2);
    b = sqrt((Node(St(2),1)-Node(St(3),1))^2 + (Node(St(2),2)-Node(St(3),2))^2 + (Node(St(2),3)-Node(St(3),3))^2);
    c = sqrt((Node(St(3),1)-Node(St(1),1))^2 + (Node(St(3),2)-Node(St(1),2))^2 + (Node(St(3),3)-Node(St(1),3))^2);
    p = (a+b+c)/2;
    Ar1 = sqrt(p*(p-a)*(p-b)*(p-c));
end

function Si = N2T(n, S) % Возвращает номера элементов, содержихих указанный узел
    [rows, cols] = find(S == n);
    %Si = rows;
    Si = sort(rows);
end

function d = Distance(Node1, Node2) % Расстояние между двумя узлами
    d = sqrt((Node1(1)-Node2(1))^2+(Node1(2)-Node2(2))^2+(Node1(3)-Node2(3))^2);
end

function di = DistanceROI(Node, ROI) % Расстояние между двумя узлами
    di = sqrt((Node(1)-ROI(:,1)).^2+(Node(2)-ROI(:,2)).^2+(Node(3)-ROI(:,3)).^2);
end

function [Fi, L, II] = Transform2D(L, nfi, nL, I) % Преобразование массива вида Node или ROI в двумерный массив для построения
    [Fi, L] = meshgrid(linspace(-pi,pi,nfi+1), linspace(-L,L,nL));
    II = reshape(I,[nfi,nL]);
    II = [II(nfi,:); II];
    II = transpose(II);
end

function [X, Y, Z, B] = Transform3D(lx, ly, lz, nx, ny, nz, Bdes) % Преобразование массива вида Node или ROI в трехмерный массив для построения
    [X, Y, Z] = ndgrid(linspace(-lx/2,lx/2,nx),linspace(-ly/2,ly/2,ny),linspace(-lz/2,lz/2,nz));
    B = reshape(Bdes,[nx,ny,nz]);    
end

function p = ScalarP(Vector1, Vector2) % Радиус вектор точки
    p = Vector1(1)*Vector2(1) + Vector1(2)*Vector2(2) + Vector1(3)*Vector2(3);
end

function s = zum(Arr) % Радиус вектор точки
    s = 0;
    for i=1:length(Arr)
        s = s + Arr(i);
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



