function [BSmag] = BSmag_init()
%---------------------------------------------------

%  OUTPUTS:
%    BSmag   = Initialized BSmag data structure
%      BSmag.Nfilament      = Number of filaments
%---------------------------------------------------

BSmag.Nfilament = 0; %Number of source filament 

% Open default figure to plot source points and field points
figure(1), hold on, grid on, box on, axis equal
xlabel('x [m]'), ylabel('y [m]'), zlabel('z [m]')
view(3), axis tight