function P_lm = legendre_bin(n,m,x)

% particularizamos el caso especial de n=0
if n == 0 
    if isvector(x)
        P_lm = ones(1,length(x),class(x));
    else
        P_lm = ones(size(x),class(x));
    end
    return; %hace que se salga de esta función al ser igual a cero
end

% Generamos el polinomio
Pol=zeros(1,2*n+1); %generamos el vector donde vamos a guardar los coeficientes del polinomio
for k=1:n+1
    ind=2*(n-k+1)+1;
    Pol(ind)=((-1)^(k-1)*factorial(n))/(factorial(k-1)*factorial(n-k+1));
end

%Volteamos el polinomio para derivarlo
Pol=fliplr(Pol);

%Derivamos el polinomio
for i=1:n+m
    Pol_der=polyder(Pol);
    Pol=Pol_der;
end

%Calculamos las constantes del polimonio asociado y los valores de la
%derivada
Cte=((-1)^m*(1-x.^2).^(m/2))/(2^n*factorial(n));
Pol=fliplr(Pol);
val_der=0;

for i=1:length(Pol)
    der=Pol(i)*x.^(i-1);
    val_der=val_der+der;
end
A=val_der;

P_lm=Cte.*A;
end

