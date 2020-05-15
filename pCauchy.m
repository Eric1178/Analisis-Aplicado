function [pC] = pCauchy(B, g, delta)
    
    if dot(g, B * g) <= 0
        tau = 1;
    else
        tau = min([norm(g)^3 / (delta * dot(g, B * g)), 1]);
    end
    
    pC = -tau * delta / norm(g) * g;
end