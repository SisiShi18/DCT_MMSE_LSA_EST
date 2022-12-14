function [xi_hat, gamma_k] = est_xi_lap_ref(noisy_pow ,noise_pow, xi_min, alpha2 ,n,Xk_prev1)
% reference :
% 'Speech enhancement using an MMSE short time DCT coefficients estimator
% with supergaussian speech modeling', 2007, Zou and Zhang
%
%     gamma_k=min(noisy_pow./noise_pow,40);
    gamma_k= max(min(noisy_pow./(noise_pow),1000),0.001);

    if n==1 
        xi_hat=alpha2+(1-alpha2)*sqrt(max(gamma_k-1,0));
    else
        xi_hat=alpha2*Xk_prev1./sqrt(noise_pow) + (1-alpha2)*sqrt(max(gamma_k-1,0));
        xi_hat=max(sqrt(xi_min),xi_hat);
    end
end