function [D] = my_idct(X,varargin)
% Description: 
% Create DCT (1-4) forward coefficients and inverse DCT to reconstruct the
% signal back to the original
% Inputs:
%   X: input data
%   type : string, '1'/'2'/'3'/'4' 
% Outputs:
%   D: output data (frames), dct spectrum
    in = inputParser;
    addRequired(in,'X',@(x) isnumeric(x));
    
    default_opt = '2';
    valid_opt = {'1','2','3','4','one','two','three','four'};
    check_opt = @(x) any(validatestring(x,valid_opt));
    addParameter(in,'type',default_opt,check_opt);


    in.parse(X,varargin{:});

    X = in.Results.X;
    type = in.Results.type;
    
    N = size(X,2);
    
    switch type
        case {'1','one'}
            C_mat = DCT_mat(N-1,'type','1');
        case {'2','two'}
            C_mat = DCT_mat(N,'type','2');
        case {'3','three'}
            C_mat = DCT_mat(N,'type','3');
        case {'4','four'}
            C_mat = DCT_mat(N,'type','4');
        otherwise
            fprintf('Invalid inputs.\n');
            return
    end
    D = X*C_mat';
end

