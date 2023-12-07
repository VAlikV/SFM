function B = BSL(Coil, I, ROI)
    B = zeros(length(ROI(:,1)),3);
    for i = 1:(length(Coil(:,1))-1)
        dlx = Coil(i+1,1) - Coil(i,1);
        dly = Coil(i+1,2) - Coil(i,2);
        dlz = Coil(i+1,3) - Coil(i,3);
        [dx, dy, dz] = VectorDistanceROI(Coil(i,:), ROI);
        B(:,1) = B(:,1) + (10^(-7))*(I*(dly*dz-dlz*dy))./(sqrt(dx.^2 + dy.^2 + dz.^2)).^3;
        B(:,2) = B(:,2) + (10^(-7))*(I*(dlz*dx-dlx*dz))./(sqrt(dx.^2 + dy.^2 + dz.^2)).^3;
        B(:,3) = B(:,3) + (10^(-7))*(I*(dlx*dy-dly*dx))./(sqrt(dx.^2 + dy.^2 + dz.^2)).^3;
    end
end