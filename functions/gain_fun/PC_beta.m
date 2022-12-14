function gain = PC_beta(xi, gamma_k, beta, p)
% Perceptually weighted beta-order spectrual amplitude Bayecian estimator for
% phase compensated speech enhancement 2021, Applied Acoustics (PC-beta)
% eq. (21) 

v = (xi./(1+xi)).*gamma_k;

N = ConflHyperGeomFun(-(beta+p)/2,1,-v).*gamma((p+beta)/2+1);
D = ConflHyperGeomFun(-p/2,1,-v).*gamma(p/2+1);

c = sqrt(v)./gamma_k;

gain = real(N./D);
gain = gain.^(1./beta);
gain = gain.*c;

W=(xi./(1+xi));
I=find(v>700); %700?
gain(I)=W(I);
%limit G for problematic cases
I=find(~isfinite(gain));
gain(I)=1;
I=find(gain>1e4);
gain(I)=1e4;

end
