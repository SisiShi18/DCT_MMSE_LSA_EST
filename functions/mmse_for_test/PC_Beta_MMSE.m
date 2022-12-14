function [xfinal,sig_mod,n_mod]= PC_Beta_MMSE(varargin)
% Description:
%   short-time (Fourier) Spectral Amplitude estimator with
%   complex Gaussian speech/noise priori.
%--------------------------------------------------------------------------
%-------------------------- Input validation ------------------------------
%--------------------------------------------------------------------------
in = inputParser;
addParameter(in,'noisy',@(x) isnumeric(x)); % noisy speech signal
addOptional(in,'noise',@(x) isnumeric(x));  % noise signal for ideal variance estimate
addOptional(in,'clean',@(x) isnumeric(x));  % clean signal for ideal variance estimate

% define signal parameters and short-time parameters
addParameter(in,'Fs',@(x) isnumeric(x) && x>0 ); % Sampling frequency (samples/sec)
addParameter(in,'frame_len',32,@(x) isnumeric(x) && x>0); % frame length in samples
addParameter(in,'shift_len',8,@(x) isnumeric(x) && x>0); % frame shift in samples
addParameter(in,'xi_min',10^(-30/10),@(x) isnumeric(x) && x>0); %in the text book is -15dB. p219
addParameter(in,'gain_min',eps,@(x) isnumeric(x) && x>0); %in the text book is -15dB. p219

% synthesis window type
default_opt = 'modified';
valid_opt = {'modified','hamming','rect','hanning'};
check_opt = @(x) any(validatestring(x,valid_opt));
addParameter(in,'synWinType',default_opt,check_opt);

% a-priori SNR estimation method
default_opt = 'DD'; % Decision-directed
valid_opt = {'ideal','DD','noise_ideal'};
check_opt = @(x) any(validatestring(x,valid_opt));
addParameter(in,'xi_type',default_opt,check_opt);

default_opt = 'normal';
valid_opt = {'norm','normal','gaussian'};
check_opt = @(x) any(validatestring(x,valid_opt));
addParameter(in,'speech_prior',default_opt,check_opt);

in.parse(varargin{:});

noisy = in.Results.noisy;
noise = in.Results.noise;
clean = in.Results.clean;
xi_min = in.Results.xi_min;
gain_min = in.Results.gain_min;

Fs = in.Results.Fs;
frame_len = in.Results.frame_len;
shift_len = in.Results.shift_len;
inputs.synWinType = in.Results.synWinType;
speech_prior = in.Results.speech_prior;
xi_type = in.Results.xi_type;
%%
f_inputs.frameLen = frame_len; % frame length in ms
f_inputs.shiftLen = shift_len; % frame shift in ms
f_noisy = feature_Class(noisy,Fs,f_inputs);
f_noise = feature_Class(noise,Fs,f_inputs);
f_clean = feature_Class(clean,Fs,f_inputs);
%% -------------------Apply Phase Compensation procedure-------------------
lambda = 3.74; % scale of compensation 
[psc_sig] = psc(noisy, Fs, frame_len, shift_len, 'G&L', lambda); 
f_psc = feature_Class(psc_sig,Fs,f_inputs);
%--------------------------------------------------------------------------
clean_mag = f_clean.mag;
noise_mag = f_noise.mag;
%--------------------------------------------------------------------------
% get noisy STDFT spectra
% noisy_mag = f_psc.mag; % short-time magnitude spectrum
noisy_mag = f_noisy.mag; % short-time magnitude spectrum
noisy_pow = noisy_mag.^2;
%--------------------------------------------------------------------------
pha_noise = f_noise.pha;
pha_sig = f_clean.pha;
%--------------------------------------------------------------------------
%%
inputs.Ts = f_noisy.Ts; %window-shift length in discrete time
beta = 0.5; p = 0.1;

switch xi_type
    case 'ideal'
        noise_pow = noise_mag.^2;
        clean_pow = clean_mag.^2;
        xi_frames = clean_pow./noise_pow;
        xi_frames = max(xi_min,xi_frames);

        gamma_k=min(noisy_pow./noise_pow,40);     
       
        switch speech_prior
            case {'norm','normal','gaussian'}
                gain = PC_beta(xi_frames, gamma_k, beta, p);
        end
        
        gain = max(gain,gain_min);
        noisy_mag = noisy_mag.*gain;
    case {'DD','noise_ideal'}
        %Set parameters for decistion-directed xi-estimate in STDFT domain
        alpha=0.98; 
        
        [Nframes,~] = size(noisy_mag);
        %% Initialize noise sample mean vector using first 6 frames- which assumed to be noise/silence.
        % initialize noise variance (approx) for decision-directed method
        if strcmp(xi_type ,'noise_ideal')
            noise_pow = noise_mag.^2;
        else       
            noise_pow = init_noise(noisy_mag,frame_len,shift_len);
            qp=struct;
            [noise_pow]=estnoiseg_dft(noisy_pow,shift_len.*1e-3,qp,noise_pow);
        end
        
        for n=1:Nframes
            %% ========================Calculate posteriori SNR========================
            % posteriori SNR : gamma_k = R_k^2/lambda_D;
            % where R_k is the noisy magnitude, lambda_D is the noise variance
            gamma_k= max(min(noisy_pow(n,:)./noise_pow(n,:),1000),0.001);
            %% =============decision-direct estimate of a priori SNR===================
            % alpha control the trade-off between speech distortion and residual noise.
            if n==1  % initial condition for recursive computation
                xi=alpha+(1-alpha)*max(gamma_k-1,0);
            else           
                xi=alpha*Xk_prev./noise_pow(n,:) + (1-alpha)*max(gamma_k-1,0);
                xi=max(xi_min,xi);  % limit ksi to -25 dB, but in the text book is -15dB. p219
            end   
            %% update/estimate the noise spectrum using VAD algorithm or noise-estimation alghorithm
            switch speech_prior
                case {'norm','normal','gaussian'}
                    gain = PC_beta(xi, gamma_k, beta, p);
            end
            
            gain = max(gain,gain_min);
            
            noisy_mag(n,:) = noisy_mag(n,:).*gain;
            
            Xk_prev = noisy_mag(n,:).^2;
            %--------------------------------------------------------------------------
            noise_mag(n,:) = noise_mag(n,:).*gain;
            clean_mag(n,:) = clean_mag(n,:).*gain;
            %--------------------------------------------------------------------------
        end   
end
%-----------------------------------------------------------------------------------------------------------------
xfinal = f_noisy.ISTDFT(noisy_mag,f_psc.pha,f_noisy.idx_mat,f_noisy.N_analysis,inputs);
%-----------------------------------------------------------------------------------------------------------------
sig_mod = f_noisy.ISTDFT(clean_mag,pha_sig,f_noisy.idx_mat,f_noisy.N_analysis,inputs);
n_mod = f_noisy.ISTDFT(noise_mag,pha_noise,f_noisy.idx_mat,f_noisy.N_analysis,inputs);
end

