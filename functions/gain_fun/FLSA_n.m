function gain = FLSA_n(Xi,gamma_k)
% - Calculate the gain factor for the 
%   short-time discrete Fourier transform 
%   Log Spectral Amplitude MMSE Estimator (FLSA)
% - Gaussian speech priori/Gaussian noise
% - double-precision 

    vk = gamma_k.*Xi./(Xi+1);   
    gain=((Xi./(Xi+1))).*real(exp(expint(vk)./2));
end
