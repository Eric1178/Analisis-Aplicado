%Script para Ejercicio 2.2 
%Función de Rosenbrock
clear all; close all; clc;
format long;

f=@extendedRosenbrock;
maxiter=10000;
m=3;
%Necesitamos hacer Line search, Line search con memoria limitada y rcSR1
%para n ∈{2,8,32,128} y para las últimas 5 iteraciones, calcular xk 
%,la norma del gradiente, f(xk), los errores entre el punto estacionario y la aproximación
%y el tiempo que tardo el algoritmo.

%Para n=2
n=2;
p=1; %El parametro para que la funcion Generarpunto genere el punto correcto 
x0=Generarpunto(p,n);
xmin=ones(n,1);
[~, iter1] = rcSR1(f, x0,maxiter);
[~, iter2] = lineLM_BFGS( f, x0, maxiter,m);
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
tF=zeros(5,1);
t3=zeros(5,1);
j=1;
for i=iter1-5:iter1
    tI = tic;
       [x, ~] = rcSR1(f, x0,i);
    tF = toc(tI);
    t1(j)=tF;
    der1(j)=norm(apGrad(f,x));
    fx1(j)=f(x);
    err1(j) = norm(x-xmin);
    j=j+1;
    s1(j)=i;
end
j=1;
for i=iter2-5:iter2
    tI = tic;
        [x, ~] = lineLM_BFGS( f, x0, i,m);
    tF = toc(tI);
    t2(j)=tF;
    der2(j)=norm(apGrad(f,x));
    fx2(j)=f(x);
    err2(j) = norm(x-xmin);
    j=j+1;
    s2(j)=i;
end
j=1;
for i=iter3-5:iter3
    tI = tic;
        [x, ~] = lineBGFS( f, x0, i);
    tF = toc(tI);
    t3(j)=tF;
    der3(j)=norm(apGrad(f,x));
    fx3(j)=f(x);
    err3(j) = norm(x-xmin);
    j=j+1;
    s3(j)=i;
end
%Genera la Tabla para rcSR1
fprintf('Tabla para rcSR1 para n=2')
fprintf('\n\t%s \t\t%s \t\t\t%s \t\t\t%s \t\t\t%s\n' ,'Iteración','NormGrad','Evaluada','Error','Tiempo');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s1(3),der1(2),fx1(2),err1(2),t1(2));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s1(4),der1(3),fx1(3),err1(3),t1(3));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s1(5),der1(4),fx1(4),err1(4),t1(4));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s1(6),der1(5),fx1(5),err1(5),t1(5));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s1(7),der1(6),fx1(6),err1(6),t1(6));
fprintf('\n');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n');
%Genera la Tabla para Line Search Memoria Limitada
fprintf('Tabla para Line Search Memoria Limitada para n=2')
fprintf('\n\t%s \t\t%s \t\t\t%s \t\t\t%s \t\t\t%s\n' ,'Iteración','NormGrad','Evaluada','Error','Tiempo');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s2(3),der2(2),fx2(2),err2(2),t2(2));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s2(4),der2(3),fx2(3),err2(3),t2(3));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s2(5),der2(4),fx2(4),err2(4),t2(4));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s2(6),der2(5),fx2(5),err2(5),t2(5));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s2(7),der2(6),fx2(6),err2(6),t2(6));
fprintf('\n');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n');
%Genera la Tabla para Line Search 
fprintf('Tabla para Line Search para n=2 ')
fprintf('\n\t%s \t\t%s \t\t\t%s \t\t\t%s \t\t\t%s\n' ,'Iteración','NormGrad','Evaluada','Error','Tiempo');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s3(3),der3(2),fx3(2),err3(2),t3(2));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s3(4),der3(3),fx3(3),err3(3),t3(3));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s3(5),der3(4),fx3(4),err3(4),t3(4));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s3(6),der3(5),fx3(5),err3(5),t3(5));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s3(7),der3(6),fx3(6),err3(6),t3(6));
fprintf('\n')
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n');

%Para n=8
n=8;
p=1; %El parametro para que la funcion Generarpunto genere el punto correcto 
x0=Generarpunto(p,n);
xmin=ones(n,1);
[~, iter1] = rcSR1(f, x0,maxiter);
[~, iter2] = lineLM_BFGS( f, x0, maxiter,m);
[~, iter3] = lineBGFS( f, x0, maxiter );
%Tabla para ver iteraciones maximas de cada metodo para n=8
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
tF=zeros(5,1);
t3=zeros(5,1);
j=1;
for i=iter1-5:iter1
    tI = tic;
       [x, ~] = rcSR1(f, x0,i);
    tF = toc(tI);
    t1(j)=tF;
    der1(j)=norm(apGrad(f,x));
    fx1(j)=f(x);
    err1(j) = norm(x-xmin);
    j=j+1;
    s1(j)=i;
end
j=1;
for i=iter2-5:iter2
    tI = tic;
        [x, ~] = lineLM_BFGS( f, x0, i,m);
    tF = toc(tI);
    t2(j)=tF;
    der2(j)=norm(apGrad(f,x));
    fx2(j)=f(x);
    err2(j) = norm(x-xmin);
    j=j+1;
    s2(j)=i;
end
j=1;
for i=iter3-5:iter3
    tI = tic;
        [x, ~] = lineBGFS( f, x0, i);
    tF = toc(tI);
    t3(j)=tF;
    der3(j)=norm(apGrad(f,x));
    fx3(j)=f(x);
    err3(j) = norm(x-xmin);
    j=j+1;
    s3(j)=i;
end
%Genera la Tabla para rcSR1
fprintf('Tabla para rcSR1 para n=8')
fprintf('\n\t%s \t\t%s \t\t\t%s \t\t\t%s \t\t\t%s\n' ,'Iteración','NormGrad','Evaluada','Error','Tiempo');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s1(3),der1(2),fx1(2),err1(2),t1(2));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s1(4),der1(3),fx1(3),err1(3),t1(3));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s1(5),der1(4),fx1(4),err1(4),t1(4));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s1(6),der1(5),fx1(5),err1(5),t1(5));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s1(7),der1(6),fx1(6),err1(6),t1(6));
fprintf('\n');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n');
%Genera la Tabla para Line Search Memoria Limitada
fprintf('Tabla para Line Search Memoria Limitada para n=8')
fprintf('\n\t%s \t\t%s \t\t\t%s \t\t\t%s \t\t\t%s\n' ,'Iteración','NormGrad','Evaluada','Error','Tiempo');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s2(3),der2(2),fx2(2),err2(2),t2(2));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s2(4),der2(3),fx2(3),err2(3),t2(3));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s2(5),der2(4),fx2(4),err2(4),t2(4));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s2(6),der2(5),fx2(5),err2(5),t2(5));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s2(7),der2(6),fx2(6),err2(6),t2(6));
fprintf('\n');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n');
%Genera la Tabla para Line Search 
fprintf('Tabla para Line Search para n=8 ')
fprintf('\n\t%s \t\t%s \t\t\t%s \t\t\t%s \t\t\t%s\n' ,'Iteración','NormGrad','Evaluada','Error','Tiempo');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s3(3),der3(2),fx3(2),err3(2),t3(2));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s3(4),der3(3),fx3(3),err3(3),t3(3));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s3(5),der3(4),fx3(4),err3(4),t3(4));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s3(6),der3(5),fx3(5),err3(5),t3(5));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s3(7),der3(6),fx3(6),err3(6),t3(6));
fprintf('\n')
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n');

%Para n=32
n=32;
p=1; %El parametro para que la funcion Generarpunto genere el punto correcto 
x0=Generarpunto(p,n);
xmin=ones(n,1);
[~, iter1] = rcSR1(f, x0,maxiter);
[~, iter2] = lineLM_BFGS( f, x0, maxiter,m);
[~, iter3] = lineBGFS( f, x0, maxiter );
%Tabla para ver iteraciones maximas de cada metodo para n=32
fprintf('\t\t\t\t\t\t%s' ,'Función Rosenbrock')
fprintf('\n\t%s \t\t\t%s \t\t%s\n' ,'Método', 'Iteraciones','Tamaño de n');
fprintf('\t--------------------------------------------');
fprintf('\n\t%s  \t\t\t\t%d \t\t\t\t%d', 'rcSR1', iter1,n);
fprintf('\n\t%s  \t\t%d \t\t\t\t\t%d', 'Line Seach LM', iter2,n);
fprintf('\n\t%s  \t\t%d \t\t\t\t%d', 'Line Seach', iter3,n);
fprintf('\n');
fprintf('\t--------------------------------------------');
fprintf('\n');

%Vamos a calcular cada valor mencionado anteriormente con cada metodo
%Definimos las matrices donde calcularemos cada valor
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
tF=zeros(5,1);
t3=zeros(5,1);
j=1;
for i=iter1-5:iter1
    tI = tic;
       [x, ~] = rcSR1(f, x0,i);
    tF = toc(tI);
    t1(j)=tF;
    der1(j)=norm(apGrad(f,x));
    fx1(j)=f(x);
    err1(j) = norm(x-xmin);
    j=j+1;
    s1(j)=i;
end
j=1;
for i=iter2-5:iter2
    tI = tic;
        [x, ~] = lineLM_BFGS( f, x0, i,m);
    tF = toc(tI);
    t2(j)=tF;
    der2(j)=norm(apGrad(f,x));
    fx2(j)=f(x);
    err2(j) = norm(x-xmin);
    j=j+1;
    s2(j)=i;
end
j=1;
for i=iter3-5:iter3
    tI = tic;
        [x, ~] = lineBGFS( f, x0, i);
    tF = toc(tI);
    t3(j)=tF;
    der3(j)=norm(apGrad(f,x));
    fx3(j)=f(x);
    err3(j) = norm(x-xmin);
    j=j+1;
    s3(j)=i;
end
%Genera la Tabla para rcSR1
fprintf('Tabla para rcSR1 para n=32')
fprintf('\n\t%s \t\t%s \t\t\t%s \t\t\t%s \t\t\t%s\n' ,'Iteración','NormGrad','Evaluada','Error','Tiempo');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s1(3),der1(2),fx1(2),err1(2),t1(2));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s1(4),der1(3),fx1(3),err1(3),t1(3));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s1(5),der1(4),fx1(4),err1(4),t1(4));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s1(6),der1(5),fx1(5),err1(5),t1(5));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s1(7),der1(6),fx1(6),err1(6),t1(6));
fprintf('\n');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n');
%Genera la Tabla para Line Search Memoria Limitada
fprintf('Tabla para Line Search Memoria Limitada para n=32')
fprintf('\n\t%s \t\t%s \t\t\t%s \t\t\t%s \t\t\t%s\n' ,'Iteración','NormGrad','Evaluada','Error','Tiempo');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s2(3),der2(2),fx2(2),err2(2),t2(2));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s2(4),der2(3),fx2(3),err2(3),t2(3));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s2(5),der2(4),fx2(4),err2(4),t2(4));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s2(6),der2(5),fx2(5),err2(5),t2(5));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s2(7),der2(6),fx2(6),err2(6),t2(6));
fprintf('\n');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n');
%Genera la Tabla para Line Search 
fprintf('Tabla para Line Search para n=32 ')
fprintf('\n\t%s \t\t%s \t\t\t%s \t\t\t%s \t\t\t%s\n' ,'Iteración','NormGrad','Evaluada','Error','Tiempo');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s3(3),der3(2),fx3(2),err3(2),t3(2));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s3(4),der3(3),fx3(3),err3(3),t3(3));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s3(5),der3(4),fx3(4),err3(4),t3(4));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s3(6),der3(5),fx3(5),err3(5),t3(5));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s3(7),der3(6),fx3(6),err3(6),t3(6));
fprintf('\n')
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n');

%Para n=128
n=128;
p=1; %El parametro para que la funcion Generarpunto genere el punto correcto 
x0=Generarpunto(p,n);
xmin=ones(n,1);
[~, iter1] = rcSR1(f, x0,maxiter);
[~, iter2] = lineLM_BFGS( f, x0, maxiter,m);
[~, iter3] = lineBGFS( f, x0, maxiter );
%Tabla para ver iteraciones maximas de cada metodo para n=128
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
tF=zeros(5,1);
t3=zeros(5,1);
j=1;
for i=iter1-5:iter1
    tI = tic;
       [x, ~] = rcSR1(f, x0,i);
    tF = toc(tI);
    t1(j)=tF;
    der1(j)=norm(apGrad(f,x));
    fx1(j)=f(x);
    err1(j) = norm(x-xmin);
    j=j+1;
    s1(j)=i;
end
j=1;
for i=iter2-5:iter2
    tI = tic;
        [x, ~] = lineLM_BFGS( f, x0, i,m);
    tF = toc(tI);
    t2(j)=tF;
    der2(j)=norm(apGrad(f,x));
    fx2(j)=f(x);
    err2(j) = norm(x-xmin);
    j=j+1;
    s2(j)=i;
end
j=1;
for i=iter3-5:iter3
    tI = tic;
        [x, ~] = lineBGFS( f, x0, i);
    tF = toc(tI);
    t3(j)=tF;
    der3(j)=norm(apGrad(f,x));
    fx3(j)=f(x);
    err3(j) = norm(x-xmin);
    j=j+1;
    s3(j)=i;
end
%Genera la Tabla para rcSR1
fprintf('Tabla para rcSR1 para n=128')
fprintf('\n\t%s \t\t%s \t\t\t%s \t\t\t%s \t\t\t%s\n' ,'Iteración','NormGrad','Evaluada','Error','Tiempo');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s1(3),der1(2),fx1(2),err1(2),t1(2));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s1(4),der1(3),fx1(3),err1(3),t1(3));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s1(5),der1(4),fx1(4),err1(4),t1(4));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s1(6),der1(5),fx1(5),err1(5),t1(5));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s1(7),der1(6),fx1(6),err1(6),t1(6));
fprintf('\n');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n');
%Genera la Tabla para Line Search Memoria Limitada
fprintf('Tabla para Line Search Memoria Limitada para n=128')
fprintf('\n\t%s \t\t%s \t\t\t%s \t\t\t%s \t\t\t%s\n' ,'Iteración','NormGrad','Evaluada','Error','Tiempo');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s2(3),der2(2),fx2(2),err2(2),t2(2));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s2(4),der2(3),fx2(3),err2(3),t2(3));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s2(5),der2(4),fx2(4),err2(4),t2(4));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s2(6),der2(5),fx2(5),err2(5),t2(5));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s2(7),der2(6),fx2(6),err2(6),t2(6));
fprintf('\n');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n');
%Genera la Tabla para Line Search 
fprintf('Tabla para Line Search para n=128 ')
fprintf('\n\t%s \t\t%s \t\t\t%s \t\t\t%s \t\t\t%s\n' ,'Iteración','NormGrad','Evaluada','Error','Tiempo');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s3(3),der3(2),fx3(2),err3(2),t3(2));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s3(4),der3(3),fx3(3),err3(3),t3(3));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s3(5),der3(4),fx3(4),err3(4),t3(4));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s3(6),der3(5),fx3(5),err3(5),t3(5));
fprintf('\n\t%d \t\t\t%d \t\t%d \t\t%d \t\t%.4d', s3(7),der3(6),fx3(6),err3(6),t3(6));
fprintf('\n')
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n');