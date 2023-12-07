function Coil1 = InterpCoil(Coil, n)
    Coil1 = zeros((length(Coil(:,1))+(length(Coil(:,1))-1)*(n-2)),3);
    for i=1:(length(Coil(:,1))-1)
        Coil1(((i-1)*(n-1))+1:(i*(n-1)+1),1) = linspace(Coil(i,1),Coil(i+1,1),n);
        Coil1(((i-1)*(n-1))+1:(i*(n-1)+1),2) = linspace(Coil(i,2),Coil(i+1,2),n);
        Coil1(((i-1)*(n-1))+1:(i*(n-1)+1),3) = linspace(Coil(i,3),Coil(i+1,3),n);
    end
end