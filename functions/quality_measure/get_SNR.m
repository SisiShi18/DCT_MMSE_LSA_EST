function snr = get_SNR(sig,noise,mod,varargin)
% Description: 
%   returns the signal-to-noise ratio (SNR) of a signal in linear scale 
%   or in decibels, by computing the ratio of its summed squared magnitude 
%   to that of the noise (or noise estimate). 
%   Noise must have the same dimensions as the signal. 
%  
% Inputs:
%           SIGNAL is a target signal as vector.
%
%           NOISE is a masker signal as vector, such that 
%                 length(NOISE)>=length(SIGNAL). Note that 
%                 in the case that length(NOISE)>length(SIGNAL),
%                 a vector of length length(SIGNAL) is selected 
%                 from NOISE starting at a random sample number.
%           
%           SNR is the desired signal-to-noise ratio level (dB).
%
% Outputs: 
%           NOISY is a mixture signal of SIGNAL and NOISE at given SNR.
%
%           NOISE is a scaled masker signal, such that the mixture
%                 NOISY=SIGNAL+NOISE has the desired SNR.
%
% Example:
%           % inline function for SNR calculation
%           SNR = @(signal,noisy)( 20*log10(norm(signal)/norm(signal-noisy)) );
%   
%           fs   = 16000;                       % sampling frequency (Hz)
%           freq = 1000;                        % sinusoid frequency (Hz)
%           time = [ 0:1/fs:2 ];                % time vector (noises)
%           signal = sin( 2*pi*freq*time );     % signal vector (s)
%           noise = randn( size(signal) );      % noise vector (s)
%           snr = -5;                           % desired SNR level (dB)
% 
%           % generate mixture signal: noisy = signal + noise
%           [ noisy, noise ] = addnoise( signal, noise, snr ); 
% 
%           % check the resulting signal-to-noise ratio
%           fprintf( 'SNR: %0.2f dB\n', SNR(signal,noisy) );
    in = inputParser;

    addRequired(in,'sig',@(x) isnumeric(x));
    addOptional(in,'noise',[],@(x) isnumeric(x));

    default_opt = 'SNR';
    valid_opt = {'SSNR','SNR','seg_S','seg_N','seg','global','SegSNR','freq','freq_dB'};
    check_opt = @(x) any(validatestring(x,valid_opt));

    addParameter(in,'type',default_opt,check_opt);

    in.parse(sig,noise,varargin{:});

    sig = in.Results.sig;
    noise = in.Results.noise;
    type = in.Results.type;
    
    if size(sig,1) == 1
        sig = sig';
    end
    
    if size(noise,1) == 1
        noise = noise';
    end
    
    global min_db max_db
    
    switch type
        case {'seg_S'}%segmental Speech SNR
%             sig_E =10*log10(sum(sig.^2,2));
%             r = max(sig_E);
%             idx = sig_E>=r-45;
%             sig_new = sig(idx,:);
%             mod_new = mod(idx,:);
%             snr=10*log10(sum(sig_new,2).^2./sum(mod_new,2).^2);
%             snr = mean(snr);
            
            fun = @(ref,noi)max(min(10*log10(sum(ref.^2,2)./sum(noi.^2,2)), 35), -10);
            snr = mean(fun(sig,sig-mod));
        case {'SSNR','seg','SegSNR'} % segmental SNR 
            
            fun = @(ref,noi)max(min(10*log10(sum(ref.^2,2)./sum(noi.^2,2)), max_db), min_db); 
        %snr_sum = sum(fun(sig.^2,noise.^2),2);
%             fprintf('*SNR.m* :: *seg* \n');
            % waveform as inputs (i.e., x(n),s(n))
            % segmental SNR measured in time domain
%             snr_min = -35; %23 (20 seams too small when max is 30)
%             snr_max = 10; %15
%             
%             r = max(10*log10(norm(sig)))-35;
%             snr = mean(min(max(10*log10(norm(sig)), r), max_db)-min(max(10*log10(norm(noise)), r), max_db));     
%% Non-speech frames were excluded from the calculation
         %   snr_new = min(max(snr, snr_min), snr_max);
%             gamma = max(snr)-35; % 45 was used in the paper
%             snr_new = snr(snr > gamma);
%             
%             snr_new = snr(snr >= snr_min & snr <= snr_max);
%             fprintf('no. of snr(frame wise) %i - %i = %i, %4.2f %% removed .\n', size(snr,1),size(snr_new,1),...
%                 size(snr,1)-size(snr_new,1),(size(snr,1)-size(snr_new,1))/(size(snr,1))*100);
%             disp(snr_new)
%             % C.Loiszou, p503
%             snr = min(snr , snr_max);
%             snr = max(snr, snr_min);
            snr = mean(fun(sig,sig-mod)); % average over all the frames (the mean value)
        case {'seg_N'}%segmental noise reduction
%             sig_E =10*log10(sum(sig.^2,2));
%             r = max(sig_E);
%             idx = sig_E>=r-45;
%             n_new = noise(idx,:);
%             nMod_new = mod(idx,:);
%             snr=10*log10(sum(n_new,2).^2./sum(nMod_new,2).^2);
            
            snr=10*log10(sum(noise,2).^2./sum(mod,2).^2);
            snr = mean(snr);
        case {'SNR','global'} % global SNR
            % waveform as inputs
            % gives a single value
%             fprintf('*SNR.m* :: *global* \n');
            snr = 10*log10(sum(sig.^2,2)./(sum(noise.^2,2)+eps));
        case 'freq'
            % segmental SNR(local spectral SNR) in frequency domain
            % measured with filter bank energy
            % Sub-band energy (SSE) as inputs
            % 
            % SSE = PSD*(filter_banks'); % Spectral filterband(subband) energy
            % SNR are computed at each time-frequency point
            snr = sig./max(noise,eps);
        case 'freq_dB'
            % segmental SNR(local spectral SNR) in dB
            % Sub-band energy (SSE) as inputs
%           % SNR are computed at each time-frequency point
            snr = 10*log10(sig./max(noise,eps));
    end
end