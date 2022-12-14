function gain = CSA_n(xi,gamma_k)
% - Calculate the gain factor for the 
%   short-time discrete Cosine transform 
%   Spectral Amplitude estimator (CSA)
% - Gaussian speech priori/Gaussian noise
% - double-precision 
    vk = gamma_k.*xi./(1+xi); 
    gain = (exp(-vk/2).*sqrt(2/pi./vk)+erf(sqrt(vk/2))).*xi./(1+xi);       
%     gain = (erf(sqrt(vk/2)) ...
%         + sqrt(2/pi)./sqrt(vk)./exp(vk/2)).*xi./(1+xi);
end

