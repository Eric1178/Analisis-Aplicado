function [f] = extendedRosenbrock(x)
n=length(x);
c = 100;
f = 0;
     for i = 1:n/2
         f = f+c*(x(2*i)-x(2*i-1)^2)^2 + (1-x(2*i-1))^2;
     end
end