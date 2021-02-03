 function Z_nlm = armonico_hiper(n,l,m,chi,theta,phi)

%Z_nlm=N_nl*Y_nlm*C^(l+1)_n-l
%creamos el término de normalización
sinm = @(a,b) sin(a).^b;
N_nl = sinm(chi,l).*2^(l+1/2)*sqrt((n+1)*gamma(n-l+1)./(pi*gamma(n+l+2))).*gamma(l+1);

%3D HS Y_nlm=P_lm*Y_lm

%polinomios de legendre
L=legendre_bin(l,abs(m),cos(theta));

if m>=0
    P_lm=(1).^m.*L(1,:);
else 
    P_lm=(-1).^m.*(factorial(l-abs(m))./factorial(l+abs(m))).*L(1,:);
end
%display(P_lm);
%Definición de los HS
if m>0
    Y_lm=sqrt((2*l+1)/(2*pi)*factorial(l-m)./factorial(l+m))*cos(m.*phi);
end
if m<0
    Y_lm=-sqrt((2*l+1)/(2*pi)*factorial(l+m)./factorial(l-m))*sin((m).*phi).*((-1)^m*factorial(l-m)./factorial(l+m));
end
if m==0
    Y_lm=sqrt((2*l+1)/(4*pi));
end
%display(Y_lm);
Y_nlm=P_lm.*Y_lm;

%polinomios de Gegenbauer
cln = gamma(n +l +2)./(factorial(n-l).*gamma(2*l + 2));
a = l-n;
b = n+l+2;
c = l + 3/2;
x = cos(chi);
d = (1/2)*(1-x);
Fn = hypergeom2F1(a,b,c,d-eps);
C_nl=cln*Fn;

%Generamos los HSH
Z_nlm =N_nl.*C_nl.*Y_nlm;

%Hay que normalizarlos para que sea un conjunto ortonormal de vectores
if Z_nlm*Z_nlm'~=1
    Z_nlm=Z_nlm./sqrt(Z_nlm*Z_nlm');
end

end

