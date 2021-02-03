clear all:
close all;

%Vamos a leer el mapa de la muestra a estudiar
map = ReadMRC('emd_11997.map'); %Llamar al archivo que contiene el modelo
[Nx,Ny,Nz] = size(map); %tamaño del volumen

%vamos a reducir el número de voxeles para que no esceda la RAM
downsampling = round(Nx/67);
voxel=1.63;
%reducimos nuestras coordenadas
[x,y,z] = meshgrid(1:downsampling:Nx); 
%reducimos el mapa tras el downslamping
map = map(1:downsampling:Nx,1:downsampling:Nx,1:downsampling:Nx);
%escribimos el volumen pero reducido
WriteMRC(map,voxel/downsampling,'kk2.mrc');

%trasladamos las coordenadas para colocar el volumen en el origen
x = x-mean2(x);
y = y-mean2(y);
z = z-mean2(z);

%creamos una esfera en la que todos lo valores que sean mayores o igual que
%el máximo de las coordenadas se anulen, es decir creamos la esfera de
%proyección
mask = sqrt(x.^2 + y.^2 + z.^2) <= (max(x(:)));

%realizamos un cambio de coordenadas
[r,theta,phi]=cartesianas2esfericas(x(mask),y(mask),z(mask));

%vamos a realizar la proyección estereográfica de 3D a una hiperesfera de
%4D cuyo radio lo definimos como p_o
p_o=max(r(:));
[R,chi]=ProyeccionStereografica(x(mask),y(mask),z(mask),p_o);

%Grado n qu desaeamos nuestra reconstrucción
n=100;

mapRecon = x.*0;
temp = x.*0;

for i=1:n
    for l=0:i
        for m=-l:l
        
        %Cálculo de armónicos hiperesféricos
        P = armonico_hiper(i,l,m,chi',theta',phi');
        
        %Cálculo de los coeficientes
        coeff = (P*map(mask));
        
        %Reconstruction
        temp(mask) = P*coeff;
        mapRecon = temp+mapRecon;
        
        end 
    end
end

%El error del mapeo
RMSE=sqrt(mean((map(:)-mapRecon(:)).^2));

%representación del modelo original y la reconstrucción
figure,slice(x,y,z,(map),0,0,0); shading flat;
figure,slice(x,y,z,(mapRecon),0,0,0); shading flat;

%Guardado del modelo reconstruido
WriteMRC(mapRecon,voxel/downsampling,'kk.mrc');






