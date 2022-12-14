function [xi, gamma_k] = est_xi_norm(noisy_pow ,noise_pow, xi_min, alpha ,idxFr,Xk_prev)
%% ========================Calculate posteriori SNR========================
% posteriori SNR : gamma_k = R_k^2/lambda_D;
% where R_k is the noisy magnitude, lambda_D is the noise variance
%     gamma_k=min(noisy_pow./noise_pow,40);
    gamma_k= max(min(noisy_pow./(noise_pow),1000),0.001);
%% =============decision-direct estimate of a priori SNR===================
    if idxFr==1  % initial condition for recursive computation
        xi=alpha+(1-alpha)*max(gamma_k-1,0);
    else
        xi=alpha.*Xk_prev./noise_pow + (1-alpha).*max(gamma_k-1,0);
        xi=max(xi_min,xi);  % limit ksi to -25 dB, but in the text book is -15dB. p219
    end
%     if idxFr==1  % initial condition for recursive computation
%         xi_hat=alpha+(1-alpha)*sqrt(max(gamma_k-1,0));
%     else
% %         xi_hat=alpha*Xk_prev./sqrt(noise_pow) + (1-alpha)*sqrt(max(gamma_k-1,0));
%         xi_hat=alpha*Xk_prev*sqrt(2)./sqrt(pi*noise_pow) + (1-alpha)*sqrt(max(gamma_k-1,0));
%         xi_hat=max(sqrt(xi_min),xi_hat);
%     end
%     
%     xi = xi_hat.^2;
end