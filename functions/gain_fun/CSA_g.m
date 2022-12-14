function gain = CSA_g(xi_hat,gamma_k)
% Calculate the gain factor for the 
% short-time discrete Cosine transform 
% Spectral Amplitude estimator (CSA)
% Gamma speech priori/Gaussian noise
% double-precision 
    c = pi/sqrt(2);
%      xi_hat = sqrt(3)./sqrt(Xi)./2; %faster
     xi_hat = sqrt(3)./xi_hat/2;
    
%      xi_hat = 1./xi_hat;
     delta_hat = sqrt(gamma_k); % note this  is different with Laplacian priori model
     
     M_pls = xi_hat+delta_hat;
     M_mns= abs(xi_hat-delta_hat);
     Msgn = sign(xi_hat-delta_hat);
%      
%      M_p2 = M_pls.^2/4;
%      M_m2 = M_mns.^2/4;
%     
%      S_pls = sqrt(M_pls);
%      S_mns = sqrt(M_mns);
%       
%      C_pls = sqrt(M_pls.^3);
%      C_mns = sqrt(M_mns.^3);

        M_p2 = (M_pls/2).^2;
        M_m2 = (M_mns/2).^2;

        S_pls = sqrt(M_pls);
        S_mns = sqrt(M_mns);

        C_pls = sqrt(M_pls).^3;
        C_mns = sqrt(M_mns).^3;
        
     D_mns = exp(-xi_hat.*delta_hat);
     
     I_m_n_quat = besseli(-1/4,M_m2);
     I_m_p_quat = besseli(1/4,M_m2);
     
     Ik_p_quat = besselk(1/4,M_p2);
     
     N = C_pls.*(besselk(3/4,M_p2)+...
         -Ik_p_quat)+...
         c.*D_mns.*C_mns.*(besseli(-3/4,M_m2)+...
         +I_m_p_quat...
         -Msgn.*(besseli(3/4,M_m2)+I_m_n_quat));
     
%      D = 2*delta_hat.*(S_pls.*Ik_p_quat)+...
%          2*delta_hat.*(c.*D_mns.*S_mns.*(I_m_n_quat...
%          -Msgn.*I_m_p_quat));
     
%      D(D == 0) = eps*1e-9;
     D = 2*delta_hat.*(S_pls.*Ik_p_quat+...
         c.*D_mns.*S_mns.*(I_m_n_quat...
         -Msgn.*I_m_p_quat));
%      
     D(D == 0) = eps;  
     gain = N./D;     
%      gain(isnan(gain)) =eps;
%      D = S_pls.*Ik_p_quat+...
%          c.*D_mns.*S_mns.*(I_m_n_quat...
%          -Msgn.*I_m_p_quat);
%      
%      gain = N./(2*delta_hat.*D);
end