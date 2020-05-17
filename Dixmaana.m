function [f] = Dixmaana(x)

n = length(x);
m = floor(n/3);

%Parametros para nuestra selección de parametros: G
alpha = 1;
beta = 0.125;
gamma = 0.125;
delta = 0.125;
k1 = 1;
k2 = 0;
k3 = 0;
k4 = 1;

%Partimos la funcion en cuatro sumas
sum1 = 0;
sum2 = 0;
sum3 = 0;
sum4 = 0;

%
for i = 1:n
    sum1 = sum1 + alpha*(x(i)^2)*(i/n)^k1;
    if i <= n-1
        sum2 = sum2 + beta*(x(i)^2)*((x(i+1) + x(i+1)^2)^2)*(i/n)^k2;
    end
    if i <= 2*m
        sum3 = sum3 + gamma*(x(i)^2)*(x(i+m)^4)*(i/n)^k3;
    end
    if i <= m
        sum4 = sum4 + delta*x(i)*x(i+2*m)*(i/n)^k4;
    end     
end

f = 1 + sum1 + sum2 + sum3 + sum4;

end
