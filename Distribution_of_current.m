%% Constants:
H=0.5;
n=2;
R=1;
J0=1;

%%Representation of the current distribution:
hold on
%There is a loop for each element of a cylindrical coordinate system
for rop=linspace(0.000001,1,10)
    for phip=linspace(0,2*pi,10) %If it were to be particularized in pi/2, 
                                    %we would get the current distribution
                                    %for x=0
        for zp=linspace(-H/2,H/2,10) %If it were to be particularized in 0, 
                                        %we would get the current
                                        %distribution for z=0
            %The next three formulas represent each element (x,y,z) of the
            %current distribution
            Jx= J0*pi*((R-rop)./H).*(sin((pi.*zp)./H)).*cos(phip) - ((zp./H).^n).*sin(phip);
            Jy= J0*pi*((R-rop)./H).*(sin((pi.*zp)./H)).*sin(phip) + ((zp./H).^n).*cos(phip);
            Jz= J0*(R-2.*rop)./(rop).*(cos((pi.*zp)./H));
            %Vectorial representation:
            axis1=rop*cos(phip);
            axis2=rop*sin(phip);
            quiver3(axis1,axis2,zp,Jx,Jy,Jz), axis([-5 5 -5 5 -5 5]);  
        end
    end
end
title('Distribution of current (A/m^2)'), xlabel('x(m)'),ylabel('y(m)'),zlabel('z');
hold off
