function [G1]=TabulateGainGamma2logmmse(Rprior,Rpost,nu);
% function [G1]=TabulateGainGamma2logmmse(Rprior,Rpost,nu);
% This function tabulates the gain functions   for log  MMSE estimation
% of speech spectral magnitudes, for a
% Generalized Gamma speech amplitude prior with gamma=2. For explanations,
% derivations and motivations, see:
% R. C. Hendriks ,
% "Log-spectral Magnitude MMSE Estimators under Super-Gaussian Densities",
% Interspeech, 2009
%
% INPUT variables:
% Rprior: Array of "a priori" SNR (SNRprior) values for which the gains are
% computed. NOTE: The values must be in dBs.
% Rpost: Array of "a posteriori" SNR (SNRpost) values for which the gains
% are computed. NOTE: The values must be in dBs.
% 'nu': parameter of the Generalized Gamma distribution.
% 
% OUTPUT variables:
% G1: Matrix with gain values for speech spectral amplitude estimation,
% evaluated at all combinations of a priori and a posteriori SNR in the
% input variables Rprior and Rpost. To be multiplied with the noisy amplitude.
% 
%
% Copyright 2008: Delft University of Technology, Information and
% Communication Theory Group. 
%
% Last modified: 21-9-2008.

% Convert to Rpost and Rprior to a linear scale.
Rpost=10.^(Rpost/10);
Rprior=10.^(Rprior(:)'/10);% transformed into a row vector.
G1=zeros(length(Rprior),length(Rpost));
for k=1:length(Rprior)
    SNRprior=Rprior(k);
    [g1]=gam2gains(nu,Rpost,SNRprior);
g1=    g1(:);
    G1(k,:)=g1.';
end

function [G]=gam2gains(nu,zeta,xi);

NN=15000;
Rvec=[1:NN].';
R=NN;
psi_nu=psi(nu);
gamma_nu=gamma(nu);
Z= xi.*zeta./(xi+nu);
vec=(gammaln(nu+Rvec))  + log((psi(nu+Rvec)-psi_nu))-log(2*gamma_nu)-2*((gammaln(Rvec+1)));
vec= vec*ones(1,length(zeta))+(Rvec)*log(Z);
vec2=ConflHyperGeomFun(nu,1,Z);

G=sqrt((xi./(zeta*((xi+nu))))).*exp(0.5*psi_nu).*exp(sum(exp(vec-ones(R,1)*log(vec2)),1));
ind=find((xi.*zeta./(xi+nu))>700);
if ~isempty(ind)
vec2_ind=log((gamma(1)/gamma(nu)))+((Z(ind))) +log(Z(ind).^(nu-1));
 G_hulp =sqrt((xi./(zeta(ind)*((xi+nu))))).*exp(0.5*psi_nu).*exp(sum(exp(vec(:,ind)-ones(R,1)*(vec2_ind)),1));  
G(ind)=G_hulp;
end

