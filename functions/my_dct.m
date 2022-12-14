function [D] = my_dct(X,N,varargin)
% Description: 
% Create DCT (1-4) forward coefficients and inverse DCT to reconstruct the
% signal back to the original
% Inputs:
%   X: input data (frames)
%   type : string, '1'/'2'/'3'/'4' 
%   N : Dimension 
%   direction : 'forward' ,'inverse'
% Outputs:
%   D: output data (frames), dct spectrum
% spec_bench = dct(speech_frames',f_clean.N_analysis/2)';
    in = inputParser;
    addRequired(in,'X',@(x) isnumeric(x));
    addRequired(in,'N',@(x) isnumeric(x) && x>0);
    
    default_opt = '2';
    valid_opt = {'1','2','3','4','one','two','three','four'};
    check_opt = @(x) any(validatestring(x,valid_opt));
    addParameter(in,'type',default_opt,check_opt);

    default_opt = 'f';
    valid_opt_str = {'f','i','forward','invert','inverse','F','I','b','backward'};
    check_opt = @(x) any(validatestring(x,valid_opt_str));
    addParameter(in,'direction',default_opt,check_opt);

    in.parse(X,N,varargin{:});

    X = in.Results.X;
    N = in.Results.N;
    type = in.Results.type;
    direction = in.Results.direction;
   
%     N = size(X,2);
    L = size(X,2);

    if L < N
        Y = zeros([size(X,1),N]);
        Y(1:size(X,1),1:size(X,2)) = X;
        X = Y;
    end
    
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
    
    switch direction
        case {'f','forward','F'}
            D = X*C_mat;
        case {'i','invert','inverse','I','b','backward'}
            D = X*C_mat';
    end
end

