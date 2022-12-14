function [xi_hat, gamma_k] = est_xi_gamma(noisy_pow ,noise_pow, xi_min, alpha3 ,n,Xk_prev)
%% ========================Calculate posteriori SNR========================
    % posteriori SNR : gamma_k = R_k^2/lambda_D;
    % where R_k is the noisy magnitude, lambda_D is the noise variance
%     gamma_k=min(noisy_pow./noise_pow,40);
    gamma_k= max(min(noisy_pow./(noise_pow),1000),0.001);
    if n==1  % initial condition for recursive computation
        xi_hat=alpha3+(1-alpha3)*sqrt(max(gamma_k-1,0));
    else
        xi_hat=alpha3*Xk_prev./sqrt(noise_pow) + (1-alpha3)*sqrt(max(gamma_k-1,0));
%         xi_hat=alpha3.*Xk_prev*sqrt(3)./sqrt(noise_pow)/2+ (1-alpha3)*sqrt(max(gamma_k-1,0)); %1
%         xi_hat=alpha3.*Xk_prev./sqrt(noise_pow)+ (1-alpha3)*sqrt(max(gamma_k-1,0)/3)*2;
%         xi_hat=max(sqrt(xi_min*3)/2,xi_hat);
        xi_hat=max(sqrt(xi_min),xi_hat);
    end
end