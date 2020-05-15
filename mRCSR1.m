function [xk, iter] = mRCSR1(f, x0,maxiter)
% Metodo de Región de Confianza
%
% Entradas : f ... función a optimizar 
% x0 ... (vector) punto inicial 
% itmax ... (numero natural) iteraciones maximas que realizara el programa
%
% Salida: x ... (vector) ultima aproximacion a un punto estacionario 
% iter... (numero natural) iteraciones realizadas 

    eta = 0.1;
    r = 1e-6;
    n = length(x0);
    radioMax=1.25;
    radio = radioMax;
    iter = 0;
    xk = x0;
    g = apGrad(f, xk);
    H = speye(n);
    B = H;
    tol=1e-5;
    
    while norm(g, 'inf') > tol
        %Paso 1 
        s = -H*g;
        if dot(s,g) < 0
            if norm(s) > radio
                s = radio*s/norm(s);
            end
        else
            s = pCauchy( B, g, radio);
        end
         
        %Paso 2 
        h = -(f(xk) - f(xk+s))/(dot(g,s) + 0.5*dot(s,B*s));
        gnew = apGrad(f, xk + s);
        gamma = gnew -g;
        
        %Paso 3
        if(h > eta)
            xk = xk + s;
            g = gnew;
            iter = iter + 1;
            if iter == maxiter
                break
            end
        end
        
        if h > 0.75 %P4
            if norm(s) > 0.8*radio
                radio = min(2*radio, radioMax);
            end
        elseif h < 0.1 %P5
            radio = 0.5*radio;
        end
            
        v = gamma -B*s;
        %P6
        if abs(dot(v,s)) >= r*norm(v)*norm(s)
            B = B + v/dot(v,s)*v';
            u = s-H*gamma;
            H = H + u/dot(u,gamma)*u';
        end
    end
end 
