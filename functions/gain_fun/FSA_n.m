function gain = FSA_n(xi,gammaK)
% Calculate the gain factor for the 
% short-time Discrete Fourier Transform (DFT)
% Spectral Amplitude MMSE Estimator(FSAE)
% Gaussian speech priori/Gaussian noise
% double-precision 
%--------------------------------------------------------------------------
% global Fn
% M = [xi.',gammaK.'];
% gain = Fn(M);
%--------------------------------------------------------------------------
% for n = 1:length(xi)
%     gain(n) = Fn({xi(n),gammaK(n)});
% end
% gain1 = diag(gain2D);
%--------------------------------------------------------------------------
    vk = gammaK.*xi./(1+xi);
    A=(((sqrt(pi)/2)*(sqrt(vk))).*(exp(-vk/2)))./gammaK;
    B=(1+vk).*(besseli(0,vk/2))+vk.*(besseli(1,vk/2));
    gain = A.*B;
end