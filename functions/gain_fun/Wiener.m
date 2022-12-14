function gain = Wiener(Xi)
% - Calculate the gain factor for the Wiener filter
% - Gaussian speech priori/Gaussian noise
% - double-precision  
    gain = (Xi./(Xi+1));
end


 