function [xfinal,sig_mod,n_mod]= CLSA_MMSE_PE(varargin)
% Description:
%   short-time (Cosine) Logarithmic Spectral Amplitude estimator with
%   Gaussian/Laplacian speech prior + Gaussian noise prior + Polarity
%   Estimation
%--------------------------------------------------------------------------
%-------------------------- Input validation ------------------------------
%--------------------------------------------------------------------------
in = inputParser;
addParameter(in,'noisy',@(x) isnumeric(x)); % noisy speech signal
addOptional(in,'noise',@(x) isnumeric(x));  % noise signal for ideal variance estimate
addOptional(in,'clean',@(x) isnumeric(x));  % clean signal for ideal variance estimate

% define signal parameters and short-time parameters
addParameter(in,'Fs',@(x) isnumeric(x) && x>0 ); % Sampling frequency (samples/sec)
addParameter(in,'frame_len',20,@(x) isnumeric(x) && x>0); % frame length in samples
addParameter(in,'shift_len',10,@(x) isnumeric(x) && x>0); % frame shift in samples
addParameter(in,'xi_min',10^(-35/10),@(x) isnumeric(x) && x>0); %in the text book is -15dB. p219
addParameter(in,'gain_min',eps,@(x) isnumeric(x) && x>0); %in the text book is -15dB. p219

% synthesis window type
default_opt = 'modified';
valid_opt = {'modified','hamming','rect'};
check_opt = @(x) any(validatestring(x,valid_opt));
addParameter(in,'synWinType',default_opt,check_opt);

% a-priori SNR estimation method
default_opt = 'DD'; % Decision-directed
valid_opt = {'ideal','DD','noise_ideal'};
check_opt = @(x) any(validatestring(x,valid_opt));
addParameter(in,'xi_type',default_opt,check_opt);

default_opt = 'normal';
valid_opt = {'norm','lap', 'generalized-gamma','gamma','laplace','normal','gauss'};
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
%--------------------------------------------------------------------------
f_noise = feature_Class(noise,Fs,f_inputs);
f_clean = feature_Class(clean,Fs,f_inputs);
f_noisy = feature_Class(noisy,Fs,f_inputs);
polar_sig = f_clean.polar;
inputs.Ts = f_clean.Ts; %window-shift length in discrete time
% [acc,pe_sig] = polarity_est(noisy,f_inputs,Fs,xi_min,gain_min,speech_prior,polar_sig,inputs);
% f_noisy = feature_Class(pe_sig,Fs,f_inputs);
%%  constant
noisy_mag = f_noisy.abs; % get noisy STDCT magnitude spectrum
noisy_pow = noisy_mag.^2;
%--------------------------------------------------------------------------
noise_mag = f_noise.abs;
clean_mag = f_clean.abs;

polar_noisy = f_noisy.polar;
polar_noise = f_noise.polar;
[r,c] = size(noisy_pow); % row, column size , dimension of the matrix
% number of elements
%--------------------------------------------------------------------------

[Nframes,~] = size(noisy_mag);
global CLSA_l_gain acc
%%
switch xi_type
    case 'ideal' % get ideal xi       
        noise_pow = noise_mag.^2;
        clean_pow = clean_mag.^2;
        
        xi_frames = clean_pow./noise_pow;
        xi_frames = max(xi_min,xi_frames);
        
        gamma_k = min(noisy_pow./noise_pow,40);     
        
        N = num2cell((gamma_k),2);
        
        switch speech_prior
            case {'norm','normal'}
%--------------------------------------------------------------------------                
                gain = CLSA_n(xi_frames,gamma_k);
%--------------------------------------------------------------------------       
            case {'lap','laplace'}
%--------------------------------------------------------------------------                
%                 gain = CLSA_l((xi_frames),gamma_k);
%--------------------------------------------------------------------------
                M = num2cell((xi_frames),2);
                L = cat(2,M,N);

                for i = 1:length(L)
                    gain(i,:) = lookup_gain_in_table(CLSA_l_gain,L{i,2}.',L{i,1}.',-40:1:40,-40:1:40,1);
                end
        end
        
        gain = max(gain,gain_min);
        noisy_mag= noisy_mag.*gain;
             
    case {'DD','noise_ideal'}
        
        switch speech_prior
            case {'norm','normal'}
                alpha1 = 0.98;
            case {'lap','laplace'}
                alpha1 = 0.98;
                alpha2 = 0.92;
        end
%         clean_est_frame_lap = [];
        if strcmp(xi_type , 'noise_ideal')
            noise_pow = noise_mag.^2;
            
            clean_est_frame = [];
%             T = (0:Nframes-1)/(Nframes-1)*length(noisy)/Fs;
%             
%             [~,t,pv]=fxpefac(noisy,Fs,shift_len*0.001);
%             voicedYN = interp1(t,pv,T,'linear','extrap'); % ==0: unvoiced / ==1: voiced
%             voicedYN(voicedYN<0.5)=0;
%             voicedYN(voicedYN>=0.5)=1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             nc = 0;
            acc_tot = 0;
            
            for n=1:Nframes
                %% ====================== Compute the MMSE Gain ==========================
                switch speech_prior
                    case {'norm','normal'}
                        
                        [xi,gamma_k] = est_xi_norm(noisy_pow(n,:),noise_pow(n,:),xi_min,alpha1,n,clean_est_frame);
%--------------------------------------------------------------------------                        
                        gain = CLSA_n(xi,gamma_k);
%--------------------------------------------------------------------------
%                         gain = lookup_gain_in_table(CLSA_n_gain,gamma_k,xi,-40:1:40,-40:1:40,1);
%--------------------------------------------------------------------------
                        gain = max(gain,gain_min);
                        
                        noisy_mag(n,:) = noisy_mag(n,:).*gain;
                        
                        clean_est_frame = noisy_mag(n,:).^2;

                    case {'lap','laplace'}
%                         
%                         [xi_hat,~] = est_xi_lap(noisy_pow(n,:),noise_pow(n,:), xi_min, alpha2 ,n,clean_est_frame_lap);
%--------------------------------------------------------------------------                        
                        [xi,gamma_k] = est_xi_norm(noisy_pow(n,:),noise_pow(n,:),xi_min,alpha1,n,clean_est_frame);
%--------------------------------------------------------------------------                                         
%                         gain = CLSA_l_double(xi,gamma_k);% works well but slow
%-------------------------------------------------------------------------- 
                        gain =lookup_gain_in_table(CLSA_l_gain,gamma_k,xi,...
                            -40:1:40,-40:1:40,1);
%                         % works well and fast
%--------------------------------------------------------------------------
                        gain = max(gain,gain_min);

                        noisy_mag(n,:) = noisy_mag(n,:).*gain;
                        
                        clean_est_frame = noisy_mag(n,:).^2;  
%                         clean_est_frame_lap = noisy_mag(n,:);
                end
%                 
%                 if voicedYN(n) && (n>1) %&& ~isOnsetYN

%                     gamma_logic = sqrt(xi+1) > 1;
                    gamma_logic = sqrt(gamma_k) > 1;
%                     gamma_logic = sqrt(noisy_pow(n,:)./noise_pow(n,:)) > 1;
                    
                    % gamma_logic = sqrt(noisy_pow(n,:)./estimate) > 1;
                    first = gamma_logic == 0;
                    
                    xi_logic = sqrt(xi) > 1;
%                     xi_logic = xi_hat > 1;
                    %                     xi_logic = sqrt(gamma_k-1) > 1;
                    %                     xi_logic = noisy_mag(n,:)./sqrt(estimate) > 1;
                    
                    second = first & xi_logic;
                    third = gamma_logic| second;
                    
                    sx_est = polar_noisy(n,:);
                    sx_est(~third) = sx_est(~third)*(-1);
                    polar_noisy(n,:) = sx_est;
                    temp = (polar_sig(n,:) == sx_est);
                    acc_tot = acc_tot+(sum(sum(temp))/c)*100; % accuracy
%                     nc = nc + 1;
%                 end
                %--------------------------------------------------------------------------
                        noise_mag(n,:) = noise_mag(n,:).*gain;
                        clean_mag(n,:) = clean_mag(n,:).*gain;
                %--------------------------------------------------------------------------
            end
%             acc = acc_tot/nc;
            acc = acc_tot/Nframes;
        else
            %Initialize noise sample mean vector using first 6 frames- which assumed to be noise/silence.
            noise_pow = init_noise(noisy_mag,frame_len,shift_len);
%             noise_pow = max(eps, noise_pow);
            %             noise_pow = repmat(noise_pow,size(noisy_mag,1),1);
%             qp=struct;
%             [noise_pow]=estnoiseg_dct(noisy_pow,shift_len.*1e-3,qp,noise_pow);

            clean_est_frame = [];
            %------------------------------------------------------------------
            q = 0.5; % a priori probability of speech presence:
            PH1mean = 0.5;
            alphaPH1mean = 0.97; %0.95
            alphaPSD = 0.97; %0.95
            priorFact  = q./(1-q);
            xiOptDb    = 15.6; % optimal fixed a priori SNR for SPP estimation
            xiOpt      = 10.^(xiOptDb./10);
            logGLRFact = log(1./sqrt(1+xiOpt));
            GLRexp     = xiOpt./(1+xiOpt)/2;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % You may put your fundamental frequency estimator and voiced/unvoiced classificator here: (or within the for loop if you like)
%             addpath('voicebox')

%             GET SPECTROGRAM DATA
%             [S,F,T] = spectrogram(s,w(Nw).',Nw-Ns,nfft,fs);            % MATLAB's new spectrogram function
%             [S,F,T] = specgram(s,nfft,fs,w(Nw).',Nw-Ns);               % MATLAB's depreciated spectrogram function
%             [S,F,T] = toframes(s,w,T,fs,nfft,type);

%             T = (0:Nframes-1)/(Nframes-1)*length(noisy)/Fs;
%             
%             [~,t,pv]=fxpefac(noisy,Fs,shift_len*0.001);
%             voicedYN = interp1(t,pv,T,'linear','extrap'); % ==0: unvoiced / ==1: voiced
%             voicedYN(voicedYN<0.5)=0;
%             voicedYN(voicedYN>=0.5)=1;
% % %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             nc = 0;
            acc_tot = 0;
            
%             xi_mat = zeros(size(noisy_pow));
            for n=1:Nframes
                %% ====================== Compute the MMSE Gain ==========================
                switch speech_prior
                    case {'norm','normal'}
                        
                        GLR  = priorFact .* exp(min(logGLRFact + GLRexp.*(noisy_pow(n,:)./noise_pow),100));
                        PH1   = GLR./(1+GLR); % a posteriori speech presence probability
                        
                        PH1mean = alphaPH1mean * PH1mean + (1-alphaPH1mean) * PH1;
                        
                        stuckInd = PH1mean > 0.99;
                        
                        PH1(stuckInd) = min(PH1(stuckInd),0.99);
                        
                        estimate =  PH1 .* noise_pow + (1-PH1) .* noisy_pow(n,:) ;
                        
                        noise_pow = alphaPSD *noise_pow+(1-alphaPSD)*estimate;
%                         
                        [xi,gamma_k] = est_xi_norm(noisy_pow(n,:),noise_pow,xi_min,alpha1,n,clean_est_frame);  
%                         [xi,gamma_k] = est_xi_norm(noisy_pow(n,:),noise_pow(n,:),xi_min,alpha1,n,clean_est_frame);  
                        
                        gain = CLSA_n(xi,gamma_k);
                        gain = max(gain,gain_min); 
                        
                        noisy_mag(n,:) = noisy_mag(n,:).*gain;
                        
                        clean_est_frame = noisy_mag(n,:).^2;
                        
%                         xi_mat(n,:)= xi;
                    case {'lap','laplace'}          
                        GLR  = priorFact .* exp(min(logGLRFact + GLRexp.*(noisy_pow(n,:)./noise_pow),200));
                        PH1   = GLR./(1+GLR); % a posteriori speech presence probability
                        
                        PH1mean = alphaPH1mean * PH1mean + (1-alphaPH1mean) * PH1;
                        stuckInd = PH1mean > 0.99;
                        PH1(stuckInd) = min(PH1(stuckInd),0.99);
                        estimate =  PH1 .* noise_pow + (1-PH1) .* noisy_pow(n,:) ;
                        
                        noise_pow = alphaPSD *noise_pow+(1-alphaPSD)*estimate;
% --------------------------------------------------------------------------
                        [xi,gamma_k] = est_xi_norm(noisy_pow(n,:),noise_pow,xi_min,alpha1,n,clean_est_frame);
%                         [xi,gamma_k] = est_xi_norm(noisy_pow(n,:),noise_pow(n,:),xi_min,alpha1,n,clean_est_frame);  
%                         [xi_hat,~] = est_xi_lap(noisy_pow(n,:),noise_pow, xi_min, alpha2 ,n,clean_est_frame_lap);
%--------------------------------------------------------------------------   
                        gain =lookup_gain_in_table(CLSA_l_gain,gamma_k,xi,...
                                            -40:1:40,-40:1:40,1);
%--------------------------------------------------------------------------
                        gain = max(gain,gain_min);

                        noisy_mag(n,:) = noisy_mag(n,:).*gain;
                        
                        clean_est_frame = noisy_mag(n,:).^2;
%                         clean_est_frame_lap = noisy_mag(n,:);
                        
%                         xi_mat(n,:)= xi;
%--------------------------------------------------------------------------
                end
                %--------------------------------------------------------------------------
                        noise_mag(n,:) = noise_mag(n,:).*gain;
                        clean_mag(n,:) = clean_mag(n,:).*gain;
%                 %--------------------------------------------------------------------------
%                 if voicedYN(n) && (n>1) %&& ~isOnsetYN
%                     gamma_logic = sqrt(xi+1) > 1;
%                     gamma_logic = sqrt(gamma_k) > 1;
%                       gamma_logic = sqrt(xi_hat.^2+1) > 1;
                    gamma_logic = sqrt(noisy_pow(n,:)./noise_pow) > 1;
%                     gamma_logic = sqrt(noisy_pow(n,:)./estimate) > 1;
                    first = gamma_logic == false;
                    
%                     xi_logic = sqrt(xi) > 1;
%                     xi_logic = xi_hat > 1;
%                     xi_logic = sqrt(gamma_k-1) > 1;
%                     xi_logic = sqrt(noisy_pow(n,:)./estimate-1) > 1;
%                     xi_logic = sqrt(noisy_pow(n,:)./noise_pow-1) > 1;
                    xi_logic = noisy_mag(n,:)./sqrt(noise_pow) > 1;
%                     xi_logic = noisy_mag(n,:)./sqrt(estimate) > 1;
                    
                    second = first & xi_logic;
                    third = gamma_logic| second;
%                     
                    sx_est = polar_noisy(n,:);
                    sx_est(~third) = sx_est(~third)*(-1);
                    polar_noisy(n,:) = sx_est;
                    temp = (polar_sig(n,:) == sx_est);
                    acc_tot = acc_tot+(sum(sum(temp))/c)*100; % accuracy
%                     nc = nc + 1;
%                 end
            end
%--------------------------------------------------------------------------            
%             acc = acc_tot/nc;
            acc = acc_tot/Nframes;
%--------------------------------------------------------------------------
%             gamma_logic = sqrt(xi_mat+1) > 1;

%             gamma_logic = sqrt(noisy_pow./noise_pow) > 1;
%             first = gamma_logic == 0;
% %             xi_logic = sqrt(xi_mat) > 1;
%             xi_logic = sqrt(noisy_pow./noise_pow-1) > 1;
%             second = first & xi_logic;
%             third = gamma_logic| second;
% 
%             sx_est = polar_noisy;
%             sx_est(~third) = sx_est(~third)*(-1);
%             polar_noisy= sx_est;
%             temp = (polar_sig == sx_est);
%             acc = (sum(sum(temp))/(r*c))*100; % accuracy
           
        end

end
xfinal = f_noisy.ISTDCT(noisy_mag,polar_noisy,f_noisy.idx_mat,inputs);
sig_mod = f_noisy.ISTDCT(clean_mag,polar_sig,f_noisy.idx_mat,inputs);
n_mod = f_noisy.ISTDCT(noise_mag,polar_noise,f_noisy.idx_mat,inputs);
end

function [acc,xfinal] = polarity_est(noisy,f_inputs,Fs,xi_min,gain_min,speech_prior,polar_sig,inputs)
global acc CLSA_l_gain
f_noisy = feature_Class(noisy,Fs,f_inputs);
noisy_mag = f_noisy.abs;
noisy_pow = noisy_mag.^2;
noisy_polar = f_noisy.polar;
Nframes = size(noisy_mag,1);
[r,c] = size(noisy_pow); 

noise_pow = init_noise(noisy_mag,f_inputs.frameLen,f_inputs.shiftLen);
noise_pow = max(eps, noise_pow);
%             noise_pow = repmat(noise_pow,size(noisy_mag,1),1);
qp=struct;
[noise_pow]=estnoiseg_dct(noisy_pow,f_inputs.shiftLen.*1e-3,qp,noise_pow);

switch speech_prior
    case {'norm','normal'}
        alpha1 = 0.98;
    case {'lap','laplace'}
        alpha1 = 0.98;
end
        
clean_est_frame = [];
            for n=1:Nframes
                %% ====================== Compute the MMSE Gain ==========================
                switch speech_prior
                    case {'norm','normal'}
                        
                     
                        [xi,gamma_k] = est_xi_norm(noisy_pow(n,:),noise_pow(n,:),xi_min,alpha1,n,clean_est_frame);  
                        
                        gain = CLSA_n(xi,gamma_k);
                        gain = max(gain,gain_min); 
                        
                        noisy_mag(n,:) = noisy_mag(n,:).*gain;
                        
                        clean_est_frame = noisy_mag(n,:).^2;
                        
                        xi_mat(n,:)= xi;
                    case {'lap','laplace'}
                               
                        [xi,gamma_k] = est_xi_norm(noisy_pow(n,:),noise_pow(n,:),xi_min,alpha1,n,clean_est_frame);  
   
                        gain =lookup_gain_in_table(CLSA_l_gain,gamma_k,xi,...
                                            -40:1:40,-40:1:40,1);

                        gain = max(gain,gain_min);

                        noisy_mag(n,:) = noisy_mag(n,:).*gain;
                        
                        clean_est_frame = noisy_mag(n,:).^2;
                        
                        xi_mat(n,:)= xi;

                end
            end

            gamma_logic = sqrt(gamma_k) > 1;
            first = gamma_logic == 0;
            xi_logic = sqrt(xi_mat) > 1;
            second = first & xi_logic;
            third = gamma_logic| second;

            sx_est = noisy_polar;
            sx_est(~third) = sx_est(~third)*(-1);
            noisy_polar= sx_est;
            temp = (polar_sig == sx_est);
            acc = (sum(sum(temp))/(r*c))*100; % accuracy
            
            xfinal = f_noisy.ISTDCT(noisy_mag,noisy_polar,f_noisy.idx_mat,inputs);
end
%================================E.O.F.====================================