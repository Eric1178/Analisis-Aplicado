function [alphaOp, gnew] = lineSearch(f, xk, dk, gk)
%% Parametros
alphaO = 0;
alphaN = 1;
alphaMax = 500;
c1 = 1e-4;
c2 = 0.99;

%%
%Definimos funciones que usaremos para realizar el algoritmo
PhiD0 = dot(gk,dk);
Phi = @(x) f(xk + x*dk);
Phi0 = @(y) f(xk) + c1*y*PhiD0; 
PhiD = @(z) dot(apGrad(f, xk + z*dk),dk);

%%
while alphaN > 0 && alphaN < alphaMax
    if Phi(alphaN) > Phi0(alphaN) || Phi(alphaN) >= Phi(alphaO)
        alphaAux = zoom(alphaO, alphaN, f, xk, gk, dk);
        break
        elseif abs(PhiD(alphaN)) <= -c2*PhiD0
        alphaAux = alphaNew;
        break
        elseif PhiD(alphaN) >= 0
            alphaAux = zoom(alphaN, alphaO, f, xk, gk, dk);
        break
    else
        alphaO = alphaN;
        alphaN = 2*alphaN;
    end
end
    alphaOp = alphaAux;
    gnew = apGrad(f, xk + alphaOp*dk);
end