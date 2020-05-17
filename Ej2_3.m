%Script para Ejercicio 2.3
%Función de Dixmaana
clear all; close all; clc;
format long;

f=@Dixmaana;
maxiter=10000;


%Necesitamos hacer Line search con memoria limitada 
%n ∈{240,960} y para m ∈{1,3,5,17,29}. Para cada dimensión n hacer una tabla 
%que contenga el número de iteraciones, los últimos valores del gradiente de f en el punto
%que se aproximo, f(xk) y el tiempo.

%Para n=240;
n=240;
p=2; %El parametro para que la funcion Generarpunto genere el punto correcto 
x0=Generarpunto(p,n);
j=1;
for i=1:30
    if (i==1 || i==3 || i==5 || i==17 ||i==29)
        tI=tic;
        [x, iter] = lineLM_BFGS( f, x0, maxiter,i);
        tF=toc(tI);
        t1(j)=tF;
        s1(j)=iter;
        der1(j)=norm(apGrad(f,x));
        fx1(j)=f(x);
        j=j+1;
    end
end
%Para n=960;
n=960;
p=2; %El parametro para que la funcion Generarpunto genere el punto correcto 
x0=Generarpunto(p,n);
j=1;
for i=1:30
    if (i==1 || i==3 || i==5 || i==17 ||i==29)
        tI=tic;
        [x, iter] = lineLM_BFGS( f, x0, maxiter,i);
        tF=toc(tI);
        t2(j)=tF;
        s2(j)=iter;
        der2(j)=norm(apGrad(f,x));
        fx2(j)=f(x);
        j=j+1;
    end
end

%Genera la Tabla para LineSearch con Busqueda Limitada para cada n 
fprintf('\t\t\t\t\t\t%s' ,'Función Dixmaana')
fprintf('\n');
fprintf('Tabla para LineSearch con Busqueda Limitada para n=240 y n=960')
fprintf('\n\t%s \t\t%s \t\t%s \t\t%s \t\t\t%s \t\t\t%s \t\t\t%s\n','Valor de n','Valor de m' ,'Iteraciones','NormGrad','Evaluada','Tiempo');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n\t%s \t\t\t\t%s \t\t\t%d \t\t\t\t%d \t\t%d \t\t%.4d', '240','1',s1(1),der1(1),fx1(1),t1(1));
fprintf('\n\t%s \t\t\t\t%s \t\t\t%d \t\t\t\t%d \t\t%d \t\t%.4d', '240','3',s1(2),der1(2),fx1(2),t1(2));
fprintf('\n\t%s \t\t\t\t%s \t\t\t%d \t\t\t\t%d \t\t%d \t\t%.4d', '240','5',s1(3),der1(3),fx1(3),t1(3));
fprintf('\n\t%s \t\t\t\t%s \t\t\t%d \t\t\t\t%d \t\t%d \t\t%.4d', '240','17',s1(4),der1(4),fx1(4),t1(4));
fprintf('\n\t%s \t\t\t\t%s \t\t\t%d \t\t\t\t%d \t\t%d \t\t%.4d', '240','29',s1(5),der1(5),fx1(5),t1(5));
fprintf('\n\t%s \t\t\t\t%s \t\t\t%d \t\t\t%d \t\t%d \t\t%.4d', '960','1',s2(1),der2(1),fx2(1),t2(1));
fprintf('\n\t%s \t\t\t\t%s \t\t\t%d \t\t\t%d \t\t%d \t\t%.4d', '960','3',s2(2),der2(2),fx2(2),t2(2));
fprintf('\n\t%s \t\t\t\t%s \t\t\t%d \t\t\t%d \t\t%d \t\t%.4d', '960','5',s2(3),der2(3),fx2(3),t2(3));
fprintf('\n\t%s \t\t\t\t%s \t\t\t%d \t\t\t%d \t\t%d \t\t%.4d', '960','17',s2(4),der2(4),fx2(4),t2(4));
fprintf('\n\t%s \t\t\t\t%s \t\t\t%d \t\t\t%d \t\t%d \t\t%.4d', '960','29',s2(5),der2(5),fx2(5),t2(5));
fprintf('\n');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n');