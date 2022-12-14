function [ output,coef ] = LSF( order, x,y_input )
%
% LSF : least square fit
%	LSF() is a function for first, second and third order least square fit.
%
%	inputs :
%           	order : "the order of the fitting" ;
%           	'1' - linear
%           	'2' - quadratic
%           	'3' - cubic
%               x : independent variable ;
%               y : dependent variable ; i.e., y = a0+a1*x+a2*x^2+..+an*x^n;
%	outputs :	
%   		coef : coefficients of the fitted function: a0,a1,...
%   		output : fitted output
%===========================================================================
    output = [];

    switch order
        case {'1','linear'}
            order = 1;
        case {'2','quad','quadratic'}
            order = 2;
        case {'3','cube','cubic'}
            order = 3;
    end

    x=x(:)';
    x_orders = [];
    
    if size(y_input,1) == 1
        y_input=y_input(:);
    end
    
    X= fliplr(vander(x));
    
    M = X'*X;

    y_input = X'*y_input;

    M = M(1:(order+1),1:(order+1));

    y = y_input(1:(order+1),:);

    coef = M\y;
    
    for i = 1:(order+1)
        x_orders = [x_orders; x.^(i-1)];
    end
    
    output = coef'*x_orders;
end





%     if size(y_input,1) ~= 1
%         X= fliplr(vander(x));
% 
%         M = X'*X;
% %         y_input = y_input(:,11);
%         y_input = X'*y_input;
%         M = M(1:(order+1),1:(order+1));
%         y = y_input(1:(order+1),:);
%         coef = M\y;
% 
%         for i = 1:(order+1)
%             x_orders = [x_orders; x.^(i-1)];
%         end
%         output = coef'*x_orders;
%         % output = coef(1,:)'+coef(2,:)'.*x+coef(3,:)'.*x.^2+coef(4,:)'.*x.^3;
% 
%     elseif size(y_input,1) == 1
%         y_input=y_input(:);
% 
%         X= fliplr(vander(x));
% 
%         M = X'*X;
%         y_input = X'*y_input;
%         M = M(1:(order+1),1:(order+1));
%         y = y_input(1:(order+1));
%         coef = M\y;
% 
%         for i = 1:(order+1)
%             x_orders = [x_orders; x.^(i-1)];
%         end
%         
%         output = coef'*x_orders;
%         %output = coef(1)+coef(2).*x+coef(3).*x.^2+coef(4).*x.^3;
%     end
