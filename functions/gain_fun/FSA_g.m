function gain = FSA_g(xi,gamma_k)
% Calculate the gain factor for the 
% short-time Discrete Fourier Transform (DFT)
% Spectral Amplitude MMSE Estimator(FSAE)
% Gamma speech priori/Gaussian noise
% double-precision 
%--------------------------------------------------------------------------
% Jan S.Erkelens, Richard C.Hendriks (2007) 
% Minimum Mean-Square Error Estimation of Discrete Fourier Coefficients
% With Generalized Gamma Priors
% eq.(18): for large SNRs ; gamma = 1; v = 1; for Laplacian speech priori
%--------------------------------------------------------------------------
global Fg
gain2D = Fg({xi,gamma_k});
gain = diag(gain2D);
%--------------------------------------------------------------------------
% global parabolic_x_neg3_2 parabolic_x_neg1_2 ...
%         parabolic_y_neg3_2 parabolic_y_neg1_2
%     v = 1/2;
%     h = v-1/2;
%     mu = sqrt(v*(v+1));
%     
%     q = sqrt(2*gamma_k);
%     p = mu./sqrt(2.*xi)-q;
%     
%     nu = interp1(parabolic_x_neg3_2,parabolic_y_neg3_2,p,'spline');
%     de = interp1(parabolic_x_neg1_2,parabolic_y_neg1_2,p,'spline');    
%     gain = h*nu./de./q;
%--------------------------------------------------------------------------
%     v = 1/2;
%     mu = sqrt(v*(v+1));
%     f = @(y,x,g,v,mu)(y.^v).*exp(-y.^2./(4*g)-mu*y./(2*sqrt(x.*g))).*besseli(0,y); %eq.(14)
%     n = integral(@(y)f(y,xi,gamma_k,v,mu),0,40,'ArrayValued',true);
%     d = integral(@(y)f(y,xi,gamma_k,v-1,mu),0,40,'ArrayValued',true);
%     gain = n./d./gamma_k/2;
%--------------------------------------------------------------------------    
end    
