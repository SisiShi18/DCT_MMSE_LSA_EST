function [CLSA_l_Gain] = Tabulate_CLSA_l_gain(Xi_dB,gamma_k_dB)
% - Calculate the gain factor for the 
%   short-time discrete Cosine transform 
%   Log Spectral Amplitude MMSE Estimator (CLSA)
% - Laplacian speech priori/Gaussian noise
% - double-precision 

    Xi = 10.^(Xi_dB./10);
    gamma_k= 10.^(gamma_k_dB./10);
    
    g = 0.5772156649015328606065120900824;
    one = 1.0; zero=0.0;
    l4 = log(4);
    
    my_2F2 = @(a)integral(@(x)((expm1(a.*x))./(a.*x))...
    ./realsqrt(one-x),zero,one,'ArrayValued',true);

    my_int = @(a,b)integral(@(x)((expm1(a.*x))./(a.*x))...
        .*erf(b.*realsqrt(one-x)),zero,one,'ArrayValued',true);

    xi_hat = 1./sqrt(Xi);
    delta_hat = sqrt(gamma_k/2);

    E_plus = (xi_hat'+delta_hat);
    E_minus= (xi_hat'-delta_hat);
    
    E_plus2 = realpow(E_plus,2);
    E_minus2 = realpow(E_minus,2);
    
    D_minus = expm1(-realsqrt((8./Xi)'*gamma_k))+one;
    
    I_plus= erfc(E_plus);
    I_minus=D_minus.*erfc(E_minus);
    
    A1 = (E_plus2).*my_2F2(-(E_plus2));
    B1 = (E_minus2).*my_2F2(-(E_minus2));
    
    if sum(isnan(B1(:)))
       [r,c] = find(isnan(B1));
       B1(r,c) = (B1(r,c-1)+B1(r,c+1))/2;
       B1(r,c) = (B1(r-1,c)+B1(r+1,c))/2;
    end
    
    A2 = E_plus2.*my_int(-E_plus2,E_plus);
    B2 = E_minus2.*my_int(-E_minus2,E_minus);
    
    if sum(isnan(B2(:)))
        [r,c] = find(isnan(B2));
        B2(r,c) = (B2(r,c-1)+B2(r,c+1))/2;
    end

    P_plus =  (A1-A2-l4);
    P_minus = D_minus.*(B1-B2-l4);

Q1 = (expm1(P_plus+P_minus)+1).^(1./(I_plus+I_minus));
Q2 = (gamma_k/2).*(expm1(g)+1);
 
CLSA_l_Gain = sqrt(Q1./Q2); 

ind=find(E_plus > 5.36);

if ~isempty(ind)
    vk = (Xi./(1+Xi))'*gamma_k;
    gain = realsqrt((Xi./(Xi+one))'./(2*exp(g).*gamma_k))...
                            .* exp(my_2F2(-vk/2).*vk./4);
    CLSA_l_Gain(ind) = gain(ind);               
else
    
end
