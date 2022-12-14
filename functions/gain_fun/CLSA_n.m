function gain = CLSA_n(Xi,gamma_k)
% - Calculate the gain factor for the 
%   short-time discrete Cosine transform 
%   Log Spectral Amplitude MMSE Estimator (CLSA)
% - Gaussian speech priori/Gaussian noise
% - double-precisi
g = 0.5772156649015328606065120900824;
vk = gamma_k.*Xi./(1+Xi);
% 
my_2F2 = @(a)integral(@(x)((expm1(a.*x))./(a.*x))...
    ./realsqrt(1-x),0,1,'ArrayValued',true);

% my_2F2 = @(a)integral(@(x)hypergeom(mp('1'),mp('2'),a.*x)./sqrt(1-x),mp(0),mp(1),'ArrayValued',true);

gain = realsqrt((Xi./(Xi+1))./(2*exp(g).*gamma_k))...
                            .* exp(my_2F2(-vk/2).*vk./4);
end

