%Script para Ejercicio 2.2 
%Función de Rosenbrock
clear all; close all; clc;
format long;

f=@extendedRosenbrock;
maxiter=10000;

%Necesitamos hacer Line search, Line search con memoria limitada y rcSR1
%para n={2,8,32,128} y para las últimas 5 iteraciones, calcular xk 
%,la norma del gradiente, f(xk), los errores entre el punto estacionario y la aproximación
%y el tiempo que tardo el algoritmo.

%Para n=2
n=2;
p=1; %El parametro para que la funcion generar punto, genere el punto correcto 
x0=Generarpunto(p,n);
[~, iter1] = rcSR1(f, x0,maxiter);
[~, iter2] = lineLM_BFGS( f, x0, maxiter);
[~, iter3] = lineBGFS( f, x0, maxiter );
%fprintf('El método rcSR1 aplicado a la funcion Rosenbrock, tomo %d iteraciones \n',iter1)
%fprintf('El método line search con memoria limitada aplicado a la funcion Rosenbrock, tomo %d iteraciones \n',iter2)
%fprintf('El método line search aplicado a la funcion Rosenbrock, tomo %d iteraciones \n',iter3)

%Tabla para ver iteraciones maximas de cada metodo para n=2
fprintf('\t\t\t\t\t\t%s' ,'Función Rosenbrock')
fprintf('\n\t%s \t\t\t%s \t\t%s\n' ,'Método', 'Iteraciones','Tamaño de n');
fprintf('\t--------------------------------------------');

fprintf('\n\t%s  \t\t\t\t%d \t\t\t\t%d', 'rcSR1', iter1,n);
fprintf('\n\t%s  \t\t%d \t\t\t\t%d', 'Line Seach LM', iter2,n);
fprintf('\n\t%s  \t\t%d \t\t\t\t%d', 'Line Seach', iter3,n);
fprintf('\n');
fprintf('\t--------------------------------------------');
fprintf('\n');

%Vamos a calcular cada valor mencionado anteriormente con cada metodo
%Definimos las matrices donde calcularemos cada valor
u=zeros(5,n);
y=zeros(5,n);
z=zeros(5,n);
der1=zeros(5,1);
der2=zeros(5,1);
der3=zeros(5,1);
fx1=zeros(5,1);
fx2=zeros(5,1);
fx3=zeros(5,1);
err1=zeros(5,1);
err2=zeros(5,1);
err3=zeros(5,1);
t1=zeros(5,1);
t2=zeros(5,1);
t3=zeros(5,1);
j=1;
for i=iter1-5:iter1
    t1 = tic;
        [x, iter1] = rcSR1(f, x0,maxiter);
    t2 = toc(t1);
    u(j,1)=x(1);
    u(j,2)=x(2);
    der1(j)=norm(apGrad(f,x));
    fx1(j)=f(x);
    %err1(j) = norm(x-xmin);
    j=j+1;
end
j=1;
for i=iter2-5:iter2
    t1 = tic;
        [x, iter2] = lineLM_BFGS( f, x0, maxiter);
    t2 = toc(t1);
    y(j,1)=x(1);
    y(j,2)=x(2);
    der2(j)=norm(apGrad(f,x));
    fx2(j)=f(x);
    %err2(j) = norm(x-xmin);
    j=j+1;
end
for i=iter2-5:iter2
    t1 = tic;
        [x, iter3] = lineBGFS( f, x0, maxiter );
    t2 = toc(t1);
    z(j,1)=x(1);
    z(j,2)=x(2);
    der3(j)=norm(apGrad(f,x));
    fx3(j)=f(x);
    %err3(j) = norm(x-xmin);
    j=j+1;
end
