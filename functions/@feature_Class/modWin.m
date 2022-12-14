function winfunc = modWin(Tw,Ts,varargin)
% DESCRIPTION: given the length of the window return the modified discrete window function
% INPPUT PARAMETERS:
%   Tw : frame length in samples, time domain
%   Ts : frame shift in samples, time domain
    in = inputParser;
    
    addRequired(in,'Tw',@(x) isnumeric(x) && x>0);
    addRequired(in,'Ts',@(x) isnumeric(x) && x>0);

    default_opt = 'hanning';
    valid_opt = {'hanning','hann','hamming','rect'};
    check_opt = @(x) any(validatestring(x,valid_opt));
    addParameter(in,'type',default_opt,check_opt);
    
    in.parse(Tw, Ts, varargin{:});
%     in.parse(Tw, varargin{:});
    Tw = in.Results.Tw;
    Ts = in.Results.Ts;
    type = in.Results.type;  

    switch type 
        case {'hanning','hann'}
            a = 0.5;
            b = -0.5;
        case 'hamming'
            a = 0.54;
            b = -0.46;           
        case 'rect'
            winfunc = rec_window(Tw,Ts);   
            return
        otherwise
            a = 0.54;
            b = -0.46; 
    end
    
    l = (2*pi*(0:Tw-1)')/Tw;
    phi = pi/Tw;
    wr= sqrt(Ts/Tw);
    coef = (2*wr)/(sqrt(4*a^2+2*b^2));
    winfunc = coef*(a+b*cos(l+phi));
end

function winfunc = rec_window(L,S)
     m = ones(L,1);
     winfunc = m*sqrt(S/L);
end

