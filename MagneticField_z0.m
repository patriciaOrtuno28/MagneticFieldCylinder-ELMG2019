%% Constants and observation points:
clear variables
vacuum = 4*pi*10^(-7); %Magnetic vacuum permeability
R=1; %Radius of the cylinder
H=0.5; %Height of the cylinder
%3 matrix to store each component of the magnetic field (x,y,z)
Bxz=zeros(17);
Byz=Bxz;
Bzz=Bxz;
%Indices for the matrix
i=1;
j=1;
warning('off'); 

%% Magnetic field (z=0)
%Declare the observation points with a loop for x and y, and z=0.001,
%almost 0 but not exactly, to avoid singularities

for x=linspace(-R,R,17) 
    for y=linspace(-R,R,17)
        z=0.001;
        [funx,funy,funz]=body(x,y,z); %Calling the function created at the end of this code
             
            %Storing in the matrix the corresponding value for each
            %observation point after doing integral3 for with integration
            %limits (ro', phi', z'), in that order
            Bxz(i,j) = (vacuum/(4*pi))*(integral3(funx,0,R,0,2*pi,-H/2,H/2));
            Byz(i,j) = (vacuum/(4*pi))*(integral3(funy,0,R,0,2*pi,-H/2,H/2));
            Bzz(i,j) = (vacuum/(4*pi))*(integral3(funz,0,R,0,2*pi,-H/2,H/2));
            
            i=i+1; %One new line for the matrix
            
    end
    i=1;
    j=j+1; %One new column for the matrix
end
 
%% Graphical representation of the magnetic field
%Creating the axis we are going to need
[Xaxis,Yaxis] = meshgrid(linspace(-5,5,17), linspace(-5,5,17));
Zaxis=zeros(size(Xaxis));
%Vector representation of the three components
subplot(2,2,1),quiver3(Xaxis,Yaxis,Zaxis,Bxz,Byz,Bzz),...
    xlabel('x(m)'),ylabel('y(m)'),zlabel('B (T)'),title('Magnetic Field Vector (z=0)');
%create_cylinder is a function designed to locate a cylinder of height H
%and radius R, in this case, which is the cylinder that generates our
%density current, to be able to see how the magnetic field is created
%inside and outside that volume more graphically, plotted in the same graph
%using the function hold on to overlap both of them
hold on,
subplot(2,2,1),create_cylinder(H,R),hidden off;
hold off

%Representation of the surface of each of the components to be able to see
%the exact value the maxima and minima of this function reaches, and its
%behavior
subplot(2,2,2),surf(Yaxis,Xaxis,Bxz),...
   xlabel('x(m)'),ylabel('y(m)'),zlabel('|B| (T)'),title('Magnetic Field Surface (z=0) - x component'),colorbar,hidden off;
subplot(2,2,3),surf(Yaxis,Xaxis,Byz),...
    xlabel('x(m)'),ylabel('y(m)'),zlabel('|B| (T)'),title('Magnetic Field Surface (z=0) - y component'),colorbar,hidden off;
subplot(2,2,4),surf(Yaxis,Xaxis,Bzz),...
    xlabel('x(m)'),ylabel('y(m)'),zlabel('|B| (T)'),title('Magnetic Field Surface (z=0) - z component'),colorbar,hidden off;


%% Body of the function that is going to be integrated

function[funx,funy,funz]=body(x,y,z)

%Definition of constants needed
R=1;
H=0.5;
n=2;
J0=1;
            %All functions will be represented in terms of ro',phi' and z',
            %for the cylindrical coordinate system
            
            %Position vector for x,y and z
            position_vector_x = @ (rop,phip,zp) x-rop.*cos(phip);
            position_vector_y = @ (rop,phip,zp) y-rop.*sin(phip);
            position_vector_z = @ (rop,phip,zp) z-zp;

            %Modulus of the position vector, with a small epsilon added in
            %the end to avoid possible singularities (divisions by zero)
            modulus_position_vector =@(rop,phip,zp) ((position_vector_x(rop,phip,zp)).^2+(position_vector_y(rop,phip,zp)).^2+(position_vector_z(rop,phip,zp)).^2)+0.01;
           
            %Components x,y and z of the current distribution
            Jx=@(rop,phip,zp) J0*pi*((R-rop)./H).*(sin((pi.*zp)./H)).*cos(phip) - ((zp./H).^n).*sin(phip);
            Jy=@(rop,phip,zp) J0*pi*((R-rop)./H).*(sin((pi.*zp)./H)).*sin(phip) + ((zp./H).^n).*cos(phip);
            Jz=@(rop,zp) J0*(R-2.*rop)./(rop).*(cos((pi.*zp)./H));
            
            %Crossproducts of x,y and z of the current distribution with
            %the position vector
            crossproduct_x = @(rop,phip,zp) Jy(rop,phip,zp).*position_vector_x(rop,phip,zp)-Jz(rop,zp).*position_vector_y(rop,phip,zp);
            crossproduct_y = @(rop,phip,zp) position_vector_x(rop,phip,zp).*Jz(rop,zp)-Jx(rop,phip,zp).*position_vector_z(rop,phip,zp);
            crossproduct_z = @(rop,phip,zp) Jx(rop,phip,zp).*position_vector_y(rop,phip,zp)-Jy(rop,phip,zp).*position_vector_x(rop,phip,zp);
            
            %The body of the integral for each component x,y and z
            funx = @(rop,phip,zp) (crossproduct_x(rop,phip,zp).*rop)./(modulus_position_vector(rop,phip,zp).^3);
            funy = @(rop,phip,zp) (crossproduct_y(rop,phip,zp).*rop)./(modulus_position_vector(rop,phip,zp).^3);
            funz = @(rop,phip,zp) (crossproduct_z(rop,phip,zp).*rop)./(modulus_position_vector(rop,phip,zp).^3);
            
end