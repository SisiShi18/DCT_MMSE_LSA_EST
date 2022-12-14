function mfcc = MFCC(obj,varargin)
    in = inputParser;

    addParameter(in,'uselifter',true,...
        @(x) validateattributes(x,{'logical'},{'scalar'}));

    addParameter(in,'liftercoeff',22,@(x) isnumeric(x) && x>0);
    in.parse(varargin{:});

    addParameter(in,'n_coef',12,@(x) isnumeric(x) && x>0);
    in.parse(varargin{:});

    addParameter(in,'en',true,...
        @(x) validateattributes(x,{'logical'},{'scalar'}));

    default_opt = 'raw';
    valid_opt = {'log','raw'};
    check_opt = @(x) any(validatestring(x,valid_opt));

    addOptional(in,'en_opt',default_opt,check_opt);

    addParameter(in,'delta',true,...
        @(x) validateattributes(x,{'logical'},{'scalar'}));

    addParameter(in,'acc',true,...
        @(x) validateattributes(x,{'logical'},{'scalar'}));

    in.parse(varargin{:});
    uselifter = in.Results.uselifter;
    liftercoeff = in.Results.liftercoeff;
    n_coef = in.Results.n_coef;
    fre_energy = in.Results.en;
    en_opt = in.Results.en_opt;
    delta = in.Results.delta;
    acc = in.Results.acc;

    if (uselifter == 1 || strcmp(uselifter,'true'))
        %liftered MFCC applied DCT to LSSE
        mfcc = obj.lifter(dct(obj.LSSE')',liftercoeff);
        %only kept 12 MFCC
        mfcc = mfcc(:,2:n_coef+1);
    else
        %without lifter, apply DCT to LSSE
        mfcc = dct(obj.LSSE')';
        mfcc = mfcc(:,2:n_coef+1);
    end

    if (delta == 1 || strcmp(delta,'true'))
        % MFCC delta coefficients
        delta = obj.time_diff([mfcc,obj.rawE],2);
    else
        delta = [];
    end

    if (acc == 1 || strcmp(acc,'true'))
        % MFCC acceleration coefficients
        acc = obj.time_diff(delta,2);
    else
        acc = [];
    end

    if (fre_energy == 1 || strcmp(fre_energy,'true'))
        % frame energy
        if (strcmp(en_opt,'raw'))
            % energy for each frame is calculated before any windowing or pre-emphesis
            en = obj.rawE; 
        else
            % log energy coefficient(applied pre-emphesis and windowing)
            en = obj.logE;
        end
    else
        en = [];
    end

    mfcc = [mfcc,en,delta,acc];
end