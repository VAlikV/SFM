function [] = BSmag_plot_field_points(BSmag,X,Y,Z)

%Plot M (where we want to calculate the field)
figure(1)
	plot3(X(:),Y(:),Z(:),'kx')
axis  tight