function gain = CLSA_l(Xi,gamma_k)
% function gain = CLSA_l(xi_hat,gamma_k)
% - Calculate the gain factor for the 
%   short-time discrete Cosine transform 
%   Log Spectral Amplitude MMSE Estimator (CLSA)
% - Laplacian speech priori/Gaussian noise
% - double-precision 
    g = 0.5772156649015328606065120900824;
    
    my_2F2 = @(a)integral(@(x)((expm1(a.*x))./(a.*x))...
        ./realsqrt(1-x),0,1,'ArrayValued',true);

    my_int = @(a,b)integral(@(x)((expm1(a.*x))./(a.*x))...
        .*erf(b.*realsqrt(1-x)),0,1,'ArrayValued',true);
    
    l4 = log(4);
    xi_hat = 1./sqrt(Xi);
    delta_hat = sqrt(gamma_k/2);

    E_plus = min(xi_hat+delta_hat,37.5);
    E_minus= min(xi_hat-delta_hat,37.5);
    
%     E_plus = (xi_hat+delta_hat);
%     E_minus= (xi_hat-delta_hat);

    E_plus2 = realpow(E_plus,2);
    E_minus2 = realpow(E_minus,2);

%     D_plus = min(expm1(realsqrt(2.*gamma_k./Xi))+1,1e40);
    D_minus = expm1(-realsqrt(4.*gamma_k./Xi))+1;
%     D_minus = expm1(-realsqrt(4.*gamma_k)./xi_hat)+1; % not working???
%     D_plus = min(expm1(realsqrt(4.*gamma_k./Xi))+1,1e40);
    
%     I_plus = D_plus.*erfc(E_plus);
%     I_minus = erfc(E_minus);
    
    I_plus = erfc(E_plus);
    I_minus = D_minus.*erfc(E_minus);

    A1 = (E_plus2).*my_2F2(-(E_plus2));
    B1 = (E_minus2).*my_2F2(-(E_minus2));
    B1(isnan(B1)) = 0;

    A2 = E_plus2.*my_int(-E_plus2,E_plus);
    B2 = E_minus2.*my_int(-E_minus2,E_minus);

    B2(isnan(B2)) = 0;
%     P_plus =  D_plus.*(A1-A2-l4);
%     P_minus = (B1-B2-l4);
    
    P_plus =  (A1-A2-l4);
    P_minus = D_minus.*(B1-B2-l4);

    gain = expm1((((P_plus+P_minus)./(I_plus+I_minus))-log(gamma_k./2)-g)./2)+1;
    gain(isinf(gain)) = 0;  
%--------------------------------------------------------------------------
% if(~strcmpi(class(Xi),'mp'))
%         mp.Digits(34);
%         Xi = mp(Xi);
%         gamma_k = mp(gamma_k);
% end
%     g = mp('euler');
%     one = mp('1');
%     zero = mp('0');
%     l4 = mp('log(4)');
%     
%     my_2F2 = @(a)integral(@(x)((expm1(a.*x))./(a.*x))...
%     ./realsqrt(one-x),zero,one,'ArrayValued',true);
% 
%     my_int = @(a,b)integral(@(x)((expm1(a.*x))./(a.*x))...
%         .*erf(b.*realsqrt(one-x)),zero,one,'ArrayValued',true);
%     
%     xi_hat = 1./sqrt(Xi);
%     delta_hat = sqrt(gamma_k/2);
% 
%     E_plus = min(xi_hat+delta_hat,mp('38'));
%     E_minus= min(xi_hat-delta_hat,mp('38'));
%     
%     E_plus2 = realpow(E_plus,2);
%     E_minus2 = realpow(E_minus,2);
%     
% %     D_plus = min(expm1(realsqrt(2.*gamma_k./Xi))+one,mp('1e30'));
%     D_minus = expm1(-realsqrt(4.*gamma_k./Xi))+one;
%     
% %     I_plus= D_plus.*erfc(E_plus);
%     I_plus= erfc(E_plus);
%     I_minus=D_minus.*erfc(E_minus);
%     
%     A1 = (E_plus2).*my_2F2(-(E_plus2));
%     B1 = (E_minus2).*my_2F2(-(E_minus2));
%     B1(isnan(B1)) = zero;
% 
%     A2 = E_plus2.*my_int(-E_plus2,E_plus);
%     B2 = E_minus2.*my_int(-E_minus2,E_minus);
%     
%     B2(isnan(B2)) = zero;
% %     P_plus =  D_plus.*(A1-A2-l4);
%     P_plus =  (A1-A2-l4);
%     P_minus = D_minus.*(B1-B2-l4);
%     
%     gain = exp((((P_plus+P_minus)./(I_plus+I_minus+eps('mp')))-log(gamma_k./2)-g)./2);
%     gain(isinf(gain)) = zero;
end
