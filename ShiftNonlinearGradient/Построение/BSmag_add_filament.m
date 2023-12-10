function [BSmag] = BSmag_add_filament(BSmag,Gamma,I,dGamma)
%---------------------------------------------------

%    BSmag      = Updated BSmag data structure
%      BSmag.Nfilament              = Number of filaments
%      BSmag.filament(*).*          = Filament structure
%      BSmag.filament(*).Gamma      = Filament points coordinates (x,y,z), one point per line [m,m,m]
%      BSmag.filament(*).I          = Filament current (flows from first point towards last point) [A]
%      BSmag.filament(*).dGamma     = Filament max discretization step [m]
%----------------------------------------------------

n = BSmag.Nfilament+1;
BSmag.filament(n).Gamma = Gamma;
BSmag.filament(n).I = I;
BSmag.filament(n).dGamma = dGamma;
BSmag.Nfilament = n;

%Plot P (where there is a current source)
figure(1)
	plot3(Gamma(:,1),Gamma(:,2),Gamma(:,3),'.-r')
axis tight