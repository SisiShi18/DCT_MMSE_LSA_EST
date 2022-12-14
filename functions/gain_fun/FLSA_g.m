function gain = FLSA_g(xi,gamma_k,zeta)
% - Calculate the gain factor for the 
%   short-time discrete Fourier transform 
%   Log Spectral Amplitude MMSE Estimator (FLSA)
% - generalized Gamma speech priori/Gaussian noise
% - double-precision 
% - from paper 'Log-spectral amplitude estimation with generalized Gamma
%   distributions for speech enhancement
%   hg=ConflHyperGeomFun(a,b,x)
%     M(a,b,x)
%     Input  : a  --- Parameter
%     b  --- Parameter ( b <> 0,-1,-2,... )
%     x  --- (Vector) argument
%     Output:  HG --- M(a,b,x)

tau = 10^(-4);%best value

if zeta == 2
    
    v = 0.2 ; %best value
    alpha = xi.*gamma_k./(v+xi);
    N = ConflHyperGeomFun(-v+1-tau/2,1,-alpha).*gamma(v+tau/2);
    D = ConflHyperGeomFun(-v+1,1,-alpha).*gamma(v);
    c = sqrt(alpha)./gamma_k;   
else
    v = 0.6;
    mu = sqrt(2*gamma_k)-sqrt(v*(v+1)./(2*xi));
    N = ParCylFun(-v-tau+1/2,-mu).*gamma(v+tau-1/2);
    D = ParCylFun(-v+1/2,-mu).*gamma(v-1/2);
    c = 1./sqrt(2*gamma_k);

end

    
gain = N./D;
gain = gain.^(1./tau);
gain = gain.*c;
    
end
