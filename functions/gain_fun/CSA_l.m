function gain = CSA_l(xi_hat,gamma_k)
% Calculate the gain factor for the 
% short-time discrete Cosine transform 
% Spectral Amplitude estimator (CSA)
% Laplacian speech priori/Gaussian Noise
% double-precision 

% xi_hat = 1./sqrt(Xi);
xi_hat = 1./xi_hat;
delta_hat = sqrt(gamma_k/2);

E_plus = xi_hat+delta_hat;
E_minus= xi_hat-delta_hat;

C_hat = exp(-E_minus.^2);

c1 =2./sqrt(pi);


D_plus = exp(4.*(xi_hat.*delta_hat));

I_plus = D_plus.*erfc(E_plus);
I_minus = erfc(E_minus);

gain = (C_hat.*c1...
    -(E_plus.*I_plus+E_minus.*I_minus))./(delta_hat.*(I_plus+I_minus));

end
