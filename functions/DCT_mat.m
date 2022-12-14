function [ C ] = DCT_mat(N,varargin )
    in = inputParser;

    addRequired(in,'N',@(x) isnumeric(x));

    default_opt = '2';
    valid_opt = {'1','2','3','4','one','two','three','four'};
    check_opt = @(x) any(validatestring(x,valid_opt));
    addParameter(in,'type',default_opt,check_opt);

    in.parse(N,varargin{:});

    N = in.Results.N;
    type = in.Results.type;
    
    m = 0:N-1;  % transform domain m_th component (R^N)
    m = m';
    n = 0:N-1;  % input vector range
    switch type
        case {'1','one'}           
            C = sqrt(2/N)*(cos((m*n*pi)/N));
            C(:,1) = C(:,1)*1/sqrt(2);
            C(:,end) = C(:,1)*1/sqrt(2);
            C(1,:) = C(1,:)*1/sqrt(2);
            C(end,:) = C(end,:)*1/sqrt(2);
        case {'2','two'}
            C= sqrt(2/N)*cos((m*(2*n+1)*pi)/(2*N));  
            C(1,:) = C(1,:)*1/sqrt(2);
        case {'3','three'}
            C= sqrt(2/N)*(cos((2*m+1)*n*pi/(2*N)));
            C(:,1) = C(:,1)*1/sqrt(2);
        case {'4','four'}
            C=  sqrt(2/N)*cos((2*m+1)*(2*n+1)*pi/(4*N));
        otherwise
            fprintf('Invalid inputs,\n');
            return
    end
    C = C';
end

