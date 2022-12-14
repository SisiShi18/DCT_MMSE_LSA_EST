function [xi_hat, gamma_k] = est_xi_lap(noisy_pow ,noise_pow, xi_min, alpha2 ,n,Xk_prev1)
%% ========================Calculate posteriori SNR========================
    % posteriori SNR : gamma_k = R_k^2/lambda_D;
    % where R_k is the noisy magnitude, lambda_D is the noise variance
%     gamma_k=min(noisy_pow./noise_pow,40);
    gamma_k= max(min(noisy_pow./(noise_pow),1000),0.001);
    if n==1  % initial condition for recursive computation
        xi_hat=alpha2+(1-alpha2)*sqrt(max(gamma_k-1,0));
    else
        xi_hat=alpha2*Xk_prev1./sqrt(noise_pow) + (1-alpha2)*sqrt(max(gamma_k-1,0));
%         xi_hat=alpha2*Xk_prev1./sqrt(2*noise_pow) + (1-alpha2)*sqrt(max(gamma_k-1,0)); %1
%         xi_hat=alpha2*Xk_prev1./sqrt(noise_pow) + (1-alpha2)*sqrt(max(gamma_k-1,0)*2);
%         xi_hat=max(sqrt(xi_min*2),xi_hat);
        xi_hat=max(sqrt(xi_min),xi_hat);
    end
end