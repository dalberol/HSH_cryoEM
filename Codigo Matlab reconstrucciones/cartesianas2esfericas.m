function [r,theta,phi] = cartesianas2esfericas(x,y,z)
%realiza el cambio de coordenadas cartesinanas a esféricas, todos los
%angulos están en radianes
%el radio siempre es >=0
r=sqrt(x.^2+y.^2+z.^2);
%theta va desde -pi hasta pi
theta=acos(z./r);
%hay que corregir para cuando el radio sea nulo
ind=find(r==0);
if r(ind)==0
    theta(ind)=0;
end

%phi va desde 0 hasta 2pi
phi=atan2(y,x)+pi;

end

