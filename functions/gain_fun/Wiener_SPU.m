function gain = Wiener_SPU(xi,gamma_k,qk)
% - Calculate the gain factor for the Wiener filter with SPU weighting
% - Gaussian speech priori/Gaussian noise

    qkr = qk/(1-qk);

    wiener_gain = (xi./(xi+1));
    
    v_k = wiener_gain.*gamma_k;
    
    w = 1+qkr*exp(-v_k/2).*sqrt(1+xi);
    
    gain = wiener_gain./w;
end


 