function [x]=Generarpunto(p,n)
%Genera el punto inicial para la función de Rosenbrock o Dixmaana 
%para un vector de longitud n
%p=1 es para Rosenbrock 
%p=2 es para Dixmaana
x=zeros(n,1);
if p==1
    for i=1:n
       if mod(i,2)==0
          x(i)=1; 
       else
           x(i)=-1.2;
       end
    end
else
    for i=1:n
        x(i)=2;
    end
end
end
