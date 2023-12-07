function Coil = flex(Coil,R)
    Len = length(Coil(:,1));
    for i = 1:(Len)
        y = Coil(i,2);
        Coil(i,1)=Coil(i,1);
        Coil(i,2)=R*sin(y/R);
        Coil(i,3)=R*(cos(y/R));
    end
    Coil = real(Coil);
end