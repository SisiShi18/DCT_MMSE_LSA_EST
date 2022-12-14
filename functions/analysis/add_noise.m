function [ noisy, noise_new ] = addnoise( waveform, snr, varargin )
% inline function for SNR calculation
%     SNR = @(signal,noisy)( 20*log10(norm(signal)/norm(signal-noisy)) );
% needed for older realases of MATLAB
%     randi = @(n)( round(1+(n-1)*rand) );
    in = inputParser;

    addRequired(in,'waveform',@(x) isnumeric(x));
    addRequired(in,'snr',@(x) isnumeric(x));
  
    addParameter(in,'noise_init',[],@(x) isnumeric(x));

    in.parse(waveform,snr,varargin{:});

    waveform = in.Results.waveform(:);
    noise_init = in.Results.noise_init;
    snr = in.Results.snr;

    if isempty(noise_init)
        % generate initial noise with zero mean and unit variance
        N = size(waveform,1);
        noise_init = normrnd (0 , 1 , N ,1) ; 
    end
            
    % inputs are vectors
    % ensure masker is at least as long as the target
    S = length( waveform );
    N = length( noise_init );
    if( S>N ), error( 'Error: length(signal)>length(noise)' ); end
    
    % generate a random start location in the masker signal
    R = randi(1+N-S);
    
    % extract random section of the masker signal
    noise_init = noise_init(R:R+S-1);
    % scale the masker w.r.t. to target at a desired SNR level
    % noise_new = noise_init / norm(noise_init) * norm(clean) / 10.0^(0.05*snr);
    
    speech_pow = sum(waveform.^2)/length(waveform);
    noise_init_pow = sum(noise_init.^2)/length(noise_init);
    ratio =sqrt(speech_pow/(noise_init_pow*10^(snr/10)));
    noise_new = ratio*noise_init;
    % generate the mixture signal
    noisy = waveform + noise_new;  
end
