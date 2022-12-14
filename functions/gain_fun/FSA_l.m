function gain = FSA_l(xi,gammaK)
% Calculate the gain factor for the 
% short-time Discrete Fourier Transform (DFT)
% Spectral Amplitude MMSE Estimator(FSAE)
% Laplacian speech priori/Gaussian noise
% double-precision 
%--------------------------------------------------------------------------
% beta order, 2015: eq.(28)
% Speech enhancement based on beta-order MMSE estimation of STSA and
% Laplacian speech modeling
%--------------------------------------------------------------------------
% global Fl
% gain2D = Fl({xi,gammaK});
% gain = diag(gain2D);
%--------------------------------------------------------------------------
%     global parabolic_x_neg3_2 parabolic_x_neg5_2 ...
%             parabolic_y_neg3_2 parabolic_y_neg5_2
% 
%     q = sqrt(2*gamma_k);
%     p = 2*sqrt(2./xi)./pi;
%     z = p-q;
% 
%     nu = interp1(parabolic_x_neg5_2,parabolic_y_neg5_2,z,'spline');
%     de = interp1(parabolic_x_neg3_2,parabolic_y_neg3_2,z,'spline');    
%     gain = (3/2).*nu./de./q;
%--------------------------------------------------------------------------
% global parabolic_x_neg3_2 parabolic_x_neg1_2 ...
%         parabolic_y_neg3_2 parabolic_y_neg1_2
%     v = 1;
%     h = v-1/2;
%     mu = sqrt(v*(v+1));
%     q = sqrt(2*gamma_k);
%     p = mu./sqrt(2.*xi)-q;
%     
%     nu = interp1(parabolic_x_neg3_2,parabolic_y_neg3_2,p,'spline');
%     de = interp1(parabolic_x_neg1_2,parabolic_y_neg1_2,p,'spline');   
%     
%     gain1 = h*nu./de./q;
%     
%     f = @(p,z) 2.^(p/2).*exp(-z.^2/4).*...
%         (sqrt(pi)./gamma((1-p)/2).*hypergeom(-p/2,1/2,(z.^2)/2)-...
%         sqrt(2*pi).*z./gamma(-p/2).*hypergeom((1-p)/2,3/2,(z.^2)/2));
%     
%     p1 = mu./sqrt(2.*xi);
%         
%     n1 = f(-2,p1) + ...
%         gamma_k.*3.*f(-4,p1)+...
%         gamma_k.^2.*30.*f(-6,p1)./4+...
%         (gamma_k/2).^3.*gamma(8).*f(-8,p1)./(factorial(3).^2)+...
%         (gamma_k/2).^4.*gamma(10).*f(-10,p1)./(factorial(4).^2)+...
%         (gamma_k/2).^5.*gamma(12).*f(-12,p1)./(factorial(5).^2);
%     
%     de1 = f(-1,p1) + ...
%         gamma_k.*D1(-3,p1)+...
%         gamma_k.^2.*6.*D1(-5,p1)./4+...
%         (gamma_k/2).^3.*gamma(7).*f(-7,p1)./(factorial(3).^2)+...
%         (gamma_k/2).^4.*gamma(9).*f(-9,p1)./(factorial(4).^2)+...
%         (gamma_k/2).^5.*gamma(11).*f(-11,p1)./(factorial(5).^2);
%     gain2 = nu1./de1./q;
%--------------------------------------------------------------------------   
    v = 1;
    mu = sqrt(v*(v+1));
    f = @(y,x,g,v,mu)(y.^v).*exp(-y.^2./(4*g)-mu*y./(2*sqrt(x.*g))).*besseli(0,y); %eq.(14)
    n = integral(@(y)f(y,xi,gamma_k,v,mu),0,80,'ArrayValued',true);
    d = integral(@(y)f(y,xi,gamma_k,v-1,mu),0,80,'ArrayValued',true);

    gain = n./d./gamma_k/2;
%     gain = max(gain1,gain2);
%--------------------------------------------------------------------------
end