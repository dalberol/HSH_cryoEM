function [R,chi] = ProyeccionStereografica(sx,sy,sz,p_o)

s=sqrt(sx.^2+sy.^2+sz.^2);

u1=2*p_o^2*sx./(p_o^2+s.^2);
u2=2*p_o^2*sy./(p_o^2+s.^2);
u3=2*p_o^2*sz./(p_o^2+s.^2);
u4=p_o*(s.^2-p_o^2)./(p_o^2+s.^2);

R=sqrt(u1.^2+u2.^2+u3.^2+u4.^2); %radio de la hiperesfera = p_o

%ángulo
chi=acos(u4./R);

end