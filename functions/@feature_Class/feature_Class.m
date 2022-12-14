classdef feature_Class 
    % DESCRIPTION: Class for creating features from given signals.
    % type " doc feature_Class " in command window for documentation of the Class.
    %
    % Constructor takes the signal_Class as required input, extra parameters
    % are passed as a structure, if not specified, uses default values.
    % To change certain parameters:
    % i.e.,
    %   paraStruc.N_analysis = 512
    %   paraStruct.N_fbanks = 22
    %   fea = feature_Class(signal,paraStruc);
    %
    % PARAMETERS: 
    %  signal :      signal waveform
    %  alpha :       pre-emphasis (high pass filter) coefficient
    %                (default : 0.97)
    %  N_analysis :        length of FFT
    %                (default : 1024)
    %  mag :         signal magnitude spectrum
    %  pha :         signal phase spectrum
    %  N_fbanks :  number of filter banks for feature extraction
    %                (default : 26)
    %  (default values are from Chapper's paper on phase distrotion)
    %  frameLen :   frame length in ms, time domain
    %                (default : 32)
    %  shiftLen :   shift length between frames in ms, time domain
    %                (default : 4)
    %  Tw :          frame length in samples, time domain
    %  Ts :          frame shift in samples, time domain
    %  idx_mat :     index matrix indicating the location of the frame samples in the signal
    %  pre_emphesis: signal after apply pre-emphesis
    %
    % USAGE: Compute features from the given signal
    %  Features:
    %     mag:     Signal magnitude spectrum
    %     pha:     Signal phase spectrum
    %     PSD:     Power spectral density
    %     SSC:     Spectral subband centroids (SSCS)
    %     SSE:     Spectral filterband(subband) energy.
    %     LSSE:    Spectral filterband(subband) energy in log scale
    %
    % EXAMPLES:
    %
%% PROPERTIES
    properties
        signal;         % signal which is processed by silence removal
        Fs;
        alpha;          % pre-emphasis (high pass filter) coefficient (default : 0.97)
        N_fbanks;       % number of filter banks for feature extraction (default : 26)
        normFbanks;     % if ture, area under each filter bank is the same. otherwise has the same amplitudes
        frameLen;       % frame length in ms, time domain (default : 25)
        shiftLen;       % shift length between frames in ms, time domain (default : 10)
        winType;
        preE;           % use pre-emphasis
        rmSil;            % silience removal
    end
    
    properties(SetAccess = protected)
        waveform;       % original waveform
        speech_frames;
        idx_mat;    % index matrix indicating the location
        % of the frame samples in the signal
    end
    
    properties (SetAccess = protected,Dependent)
        mag;        % signal magnitude spectrum
        pha;        % signal phase spectrum
        abs;        % Absolute spectrum
        polar;      % Polarity spectrum
        spec_dft;
        spec_dct;
        PSD;        % power spectral density
        SSC;        % spectral subband centroids (SSCS)
        SSE;        % Spectral filterband(subband) energy.
        LSSE;       %taking the log
    end 
    
    properties (SetAccess = protected,Dependent,Hidden) % Only computed when required
        Tw;         % frame length in samples (Time domain)
        Ts;         % frame shift in samples  (Time domain)
        N_analysis; % length of FFT/DCT (default : 2*Tw)
        preEmphasis; % signal after apply pre-emphesis
        E;
        logE;
        rawE;
        rawElog;
    end
%% METHODS  
    methods
        function obj =  feature_Class(waveform,Fs,varargin) % Class object constructor
            in = inputParser;
            
            addRequired(in,'waveform',@(x) isnumeric(x));
            addRequired(in,'Fs',@(x) isnumeric(x) && x>0);
    
            addParameter(in,'frameLen',32,@(x) isnumeric(x) && x>0);
            addParameter(in,'shiftLen',4,@(x) isnumeric(x) && x>0);
            addParameter(in,'alpha',0.97,@(x) isnumeric(x));
            addParameter(in,'N_fbanks',26,@(x) isnumeric(x) && x>0);
            addParameter(in,'normFbanks',false,...
                @(x) validateattributes(x,{'logical'},{'scalar'}));
            
            default_opt = 'hamming';
            valid_opt = {'hamming','hanning','rect','mod_Hamming','mod_Hanning','mod_rect'};
            check_opt = @(x) any(validatestring(x,valid_opt));
            
            addParameter(in,'winType',default_opt,check_opt);
            
            addParameter(in,'preE',false,...
                @(x) validateattributes(x,{'logical'},{'scalar'}));
            
            addParameter(in,'rmSil',false,...
                @(x) validateattributes(x,{'logical'},{'scalar'}));
            
            in.parse(waveform,Fs,varargin{:});
            
            obj.waveform = in.Results.waveform;
            obj.Fs       = in.Results.Fs;
            obj.frameLen = in.Results.frameLen;
            obj.shiftLen = in.Results.shiftLen;
            obj.alpha    = in.Results.alpha;
            obj.N_fbanks = in.Results.N_fbanks;
            obj.normFbanks = in.Results.normFbanks;
            obj.winType   = in.Results.winType;
            obj.preE = in.Results.preE;
            obj.rmSil  = in.Results.rmSil;
            
            if obj.rmSil
                % default window:'hamming', default input signal:'obj.signal'
                frame_input.signal = obj.waveform; 
                frame_input.winType = 'rect';
                [speech_frames_old,idx_mat_old] = obj.framing(frame_input);
                [obj.speech_frames,obj.idx_mat,obj.signal] = obj.rmSilence(idx_mat_old,speech_frames_old);
            else
                obj.signal = obj.waveform; 
                [obj.speech_frames,obj.idx_mat] = obj.framing; % default window:'hamming', default input signal:'obj.signal'
            end
        end
        
        function mag = get.mag (obj)
            mag = obj.STDFT(obj.speech_frames,obj.N_analysis);
        end
        
        function pha = get.pha (obj)  
            [~,pha] = obj.STDFT(obj.speech_frames,obj.N_analysis);
        end
        
        function spec = get.spec_dft (obj)
            [~,~,spec] = obj.STDFT(obj.speech_frames,obj.N_analysis);
        end
        
        function abs = get.abs (obj)            
            abs = obj.STDCT(obj.speech_frames,obj.N_analysis);
        end
        
        function polar = get.polar (obj)
            [~,polar] = obj.STDCT(obj.speech_frames,obj.N_analysis);
        end
        
        function spec = get.spec_dct (obj)
            [~,~,spec] = obj.STDCT(obj.speech_frames,obj.N_analysis);
        end
        
        function N_analysis = get.N_analysis(obj)
            % length of FFT/DCT (default : 2*Tw)
%             N_analysis = obj.Tw;
            N_analysis = 2*obj.Tw;
        end
        
        function pre = get.preEmphasis(obj)
            pre = filter([1,-obj.alpha],1,obj.signal);
        end
        
        function Tw = get.Tw (obj)
            Tw = round(obj.Fs*obj.frameLen*0.001);
        end
        
        function Ts = get.Ts (obj)
            Ts = round(obj.Fs*obj.shiftLen*0.001);
        end
        
        function PSD = get.PSD(obj)
%             PSD = obj.mag.^2/obj.Tw;
            PSD = obj.abs.^2/obj.Tw;
        end 

        function SSC = get.SSC (obj)
            % spectral subband centroids (SSCS)
            SSC = zeros (size(obj.PSD,1),obj.N_fbanks);
            w = obj.filter_banks; % create filter bank masks/weight
            f = (obj.Fs/2)*(1 : (obj.N_analysis/2)+1)/(obj.N_analysis/2+1);
            %             f = (1 : (obj.N_analysis/2))/(obj.N_analysis/2); % normalized frequency
            
            for m = 1:obj.N_fbanks
                wP = bsxfun(@times,w(m,:),obj.PSD);
                fwP = bsxfun(@times,f,wP);
                
                SSC(:,m) = sum(fwP,2)./sum(wP,2);
            end
        end
        
        function SSE = get.SSE(obj)
            % Spectral filterband(subband) energy.
            SSE = obj.PSD*(obj.filter_banks');
        end
        
        function LSSE = get.LSSE(obj)
            % Log spectral filterband energy
            sse = obj.SSE;
            LSSE = log(sse+eps); %taking the log
        end
         
        function logE = get.logE(obj)
            % log energy coefficient(applied pre-emphesis and windowing)
            logE = log(sum(obj.PSD,2));
        end 
        
        function E = get.E(obj)
            E = sum(obj.PSD,2);
        end
        
        function rawE = get.rawE(obj)
            % energy for each frame is calculated before any windowing or pre-emphesis
            rawE = sum(obj.framing('winType','rect').^2,2);
        end
        
        function rawElog = get.rawElog(obj)
            % energy for each frame is calculated before any windowing or pre-emphesis
            rawElog = log(sum(obj.framing('winType','rect').^2,2));
        end
             
        % External methods declarations
        mfcc = MFCC(obj,varargin)
        [frame_data,idx_mat]= framing(obj,varargin)
        fbanks = filter_banks(obj)
        h = imagesc(obj,varargin)
    end
    
%%  STATIC METHODS
    methods (Static)
        function [ speech_frames_new,idx_mat_new,new_speech,speech_frames_rect] = rmSilence(idx_mat_old,speech_frames_old,varargin)
            in = inputParser;
            % Only use modified version when window shift divides the
            % window length evenly
            default_opt = 'hamming';
            valid_opt = {'hamming','hanning','rect','mod_hamming','mod_hanning'};
            check_opt = @(x) any(validatestring(x,valid_opt));
            addParameter(in,'winType',default_opt,check_opt);      
            
            in.parse(varargin{:});
            winType = in.Results.winType;
% N = size(frame_clean,2);
% PSD = frame_clean.^2/N;
% E = sum(PSD,2);
% E = sum((speech_frames_old.^2),2);
            Tw = size(speech_frames_old,2);
            E = sum((speech_frames_old.^2)/Tw,2);
            perTrain = 0.2; % percentage of the data that identified as silience.
            ESort = sort(E);
            Thresh = ESort(round(perTrain * size(E,1)));
            speech_frames_rect= speech_frames_old(E >= Thresh,:);
            idx_mat_new = idx_mat_old(1:size(speech_frames_rect,1),:);
            new_speech = overlap_add(speech_frames_rect,idx_mat_new);
            
            switch winType
                % Only use modified version when window shift divides the
                % window length evenly
                case 'hamming'
                    speech_frames_new = bsxfun(@times,speech_frames_rect,hamming(Tw)'); % apply window
                case 'hanning'
                    speech_frames_new = bsxfun(@times,speech_frames_rect,hanning(Tw)');
                case 'rect'
                    %frame_data = bsxfun(@times,frame_data,rectwin(obj.Tw)');
                    speech_frames_new = bsxfun(@times,speech_frames_rect,ones(1,Tw)); % apply window
                case 'mod_hamming'
                    speech_frames_new = bsxfun(@times,speech_frames_rect,obj.modWin(Tw,Ts,'type','hamming')'); % apply window
                case 'mod_hanning'
                    speech_frames_new = bsxfun(@times,speech_frames_rect,obj.modWin(Tw,Ts,'type','hann')');
%                 case 'mod_rect'
%                     speech_frames_new = bsxfun(@times,speech_frames_rect,obj.modified_window(obj.Tw,'type','rect')');

            end
            
        end
%% Analysis STDCT      
        function [abs_spec, polar_spec, spec] = STDCT(frames,N_analysis)
            % STDCT : short-time discrete cosine transform
            %---- matlab dct function ( DCT-II)
%             spec = dct(frame_data',N_analysis)';      
            %---- my dct function (DCT-I,II,IV,VI)
            spec = my_dct(frames,N_analysis,'type','3');
            %             spec = my_dct(frame_data,N_analysis,'type','2');
            abs_spec = abs(spec);
            polar_spec = sign(spec);
        end
%% Analysis STDFT         
        function [mag, pha , spec] = STDFT(frames,N_analysis)
            % apply STDFT : short-time discrete FFT
            spec  = fft(frames, N_analysis, 2);
            % 1. use the DC component and the Nyquist frequency component
            spec = spec (:, 1:(N_analysis/2+1));
%           % or 2. discard the DC component
%             spec = spec (:, 2:(0.5*N_analysis+1));
%             % avoid dividing zeros by adding small values
%             spec(:,1) =ones(size(spec,1),1)*eps;
            
            mag = abs(spec);
            pha = angle(spec);
        end
%% Synthesis STDCT        
        function outsig = ISTDCT(abs_spec,polar,idx_mat,varargin)

            in = inputParser;
            
            addRequired(in,'abs_spec',@(x) isnumeric(x));
            addRequired(in,'polar',@(x) isnumeric(x));
            addRequired(in,'idx_mat',@(x) isnumeric(x));
            
            
            default_opt = 'modified';
            valid_opt = {'modified','hamming','rect'};
            check_opt = @(x) any(validatestring(x,valid_opt));
            addParameter(in,'synWinType',default_opt,check_opt);
            addParameter(in,'Ts',32,@(x) isnumeric(x) && x>0);
            
            in.parse(abs_spec,polar,idx_mat,varargin{:});
            
            abs_spec = in.Results.abs_spec;
            polar = in.Results.polar;
            synWinType = in.Results.synWinType;
            Ts = in.Results.Ts;  
       
            [~,Tw] = size(idx_mat);
            
            % Reconstruct frames
            Y = abs_spec.*polar;
            %% matlab dct invert function
%             frames_syn = real(idct(Y')');
            %% my dct invert function
            frames_syn = my_idct(Y,'type','3'); 
            frames_syn = frames_syn(:,1:Tw);
            
            switch synWinType
                case 'modified'
                    winfunc = modWin(Tw,Ts,'type','hann')';
                    frames_syn = bsxfun(@times,frames_syn,winfunc./(sum(winfunc.^2)));
%                      frames_syn = bsxfun(@times,frames_syn,winfunc);
                    
%                      overadd_count = zeros(1, max(idx_mat(:)));
%                      wnt_sum = zeros(1, max(idx_mat(:))); % double check with the paper
                    outsig = zeros(1, max(idx_mat(:)));
                    
                    [num_frames] = size(idx_mat);
                    
                    for n = 1:num_frames
                        outsig(idx_mat(n,:)) = outsig(idx_mat(n,:)) + frames_syn(n,:);
%                         overadd_count(idx_mat(n,:)) = overadd_count(idx_mat(n,:)) + 1;
%                         wnt_sum(idx_mat(n,:)) = wnt_sum(idx_mat(n,:))+ winfunc.^2;
                    end
                    
                    % averaging the overlapped estimates for each sample
%                      outsig = outsig./overadd_count;
%                      outsig = outsig./wnt_sum;
                    
                case 'hamming'
                    winfunc = hamming(Tw)';
                    
                    outsig = zeros(1, max(idx_mat(:)));
                    wnt_sum = zeros(1, max(idx_mat(:)));
                    [num_frames] = size(idx_mat);
                    
                    for n = 1:num_frames
                        outsig(idx_mat(n,:)) = outsig(idx_mat(n,:)) + frames_syn(n,:);
                        wnt_sum(idx_mat(n,:)) = wnt_sum(idx_mat(n,:))+ winfunc;
                    end
                    
                    outsig = outsig./wnt_sum;
                case 'rect'
                    winfunc = ones(size(Tw))';
                    outsig = zeros(1, max(idx_mat(:)));
                    wnt_sum = zeros(1, max(idx_mat(:)));
                    [num_frames] = size(idx_mat);
                    
                    for n = 1:num_frames
                        outsig(idx_mat(n,:)) = outsig(idx_mat(n,:)) + frames_syn(n,:);
                        wnt_sum(idx_mat(n,:)) = wnt_sum(idx_mat(n,:))+ winfunc;
                    end
                    outsig = outsig./wnt_sum;
            end
        end
%% Synthesis STDFT
        function outsig = ISTDFT(mag,pha,idx_mat,N_analysis,varargin)
            
            in = inputParser;
            
            addRequired(in,'mag',@(x) isnumeric(x)); % or abs spec
            addRequired(in,'pha',@(x) isnumeric(x)); % or polar spec
            addRequired(in,'idx_mat',@(x) isnumeric(x));
            addRequired(in,'N_analysis',@(x) isnumeric(x) && x>0);
            
            
            default_opt = 'modified';
            valid_opt = {'modified','hamming','rect'};
            check_opt = @(x) any(validatestring(x,valid_opt));
            addParameter(in,'synWinType',default_opt,check_opt);        
            addParameter(in,'Ts',32,@(x) isnumeric(x) && x>0);
            
            in.parse(mag,pha,idx_mat,N_analysis,varargin{:});
            
            mag = in.Results.mag;
            pha = in.Results.pha;
            synWinType = in.Results.synWinType;
            Ts = in.Results.Ts;
            
            [~,Tw] = size(idx_mat);
            % Mirror magnitude spectrum
            % Reconstruct frames
            mag = [mag mag(:,(end-1):-1:2)];
            pha = [pha -pha(:,(end-1):-1:2)];

%             % According to the analysis stage,method 2,the DC components were discarded
%             mag = [zeros(size(idx_mat,1),1) mag mag(:,(end-1):-1:1)];
%             pha = [zeros(size(idx_mat,1),1) pha -pha(:,(end-1):-1:1)];
            Y = mag.*exp(1i*pha);
            frames_syn = real(ifft(Y,N_analysis,2));
            frames_syn = frames_syn(:,1:Tw);
           
            switch synWinType
                case 'modified'
                    winfunc = modWin(Tw,Ts,'type','hann')';
                    frames_syn = bsxfun(@times,frames_syn,winfunc./(sum(winfunc.^2)));
%                     frames_syn = bsxfun(@times,frames_syn,winfunc);
                    % double check with the paper 
%                     overadd_count = zeros(1, max(idx_mat(:)));
%                     wnt_sum = zeros(1, max(idx_mat(:)));
                    
                    outsig = zeros(1, max(idx_mat(:)));
                   
                    [num_frames] = size(idx_mat);

                    for n = 1:num_frames
                        outsig(idx_mat(n,:)) = outsig(idx_mat(n,:)) + frames_syn(n,:);
%                         wnt_sum(idx_mat(n,:)) = wnt_sum(idx_mat(n,:))+ winfunc.^2;
%                         overadd_count(idx_mat(n,:)) = overadd_count(idx_mat(n,:)) + 1;
                    end

                    % averaging the overlapped estimates for each sample
%                      outsig = outsig./overadd_count;
%                      outsig = outsig./wnt_sum;
                case 'hamming'
                    winfunc = hamming(Tw)';
                    
%                     overadd_count = zeros(1, max(idx_mat(:)));
                    outsig = zeros(1, max(idx_mat(:)));
                    wnt_sum = zeros(1, max(idx_mat(:)));
                    [num_frames] = size(idx_mat);

                    for n = 1:num_frames
                        outsig(idx_mat(n,:)) = outsig(idx_mat(n,:)) + frames_syn(n,:);
                        wnt_sum(idx_mat(n,:)) = wnt_sum(idx_mat(n,:))+ winfunc;
%                         overadd_count(idx_mat(n,:)) = overadd_count(idx_mat(n,:)) + 1;
                    end

                    % averaging the overlapped estimates for each sample
%                     outsig = outsig./overadd_count;
%                     wnt_sum = wnt_sum/overadd_count;
                    outsig = outsig./wnt_sum;
                case 'rect'
                    winfunc = ones(size(Tw))';
                    
%                     overadd_count = zeros(1, max(idx_mat(:)));
                    outsig = zeros(1, max(idx_mat(:)));
                    wnt_sum = zeros(1, max(idx_mat(:)));
                    [num_frames] = size(idx_mat);
                    
                    for n = 1:num_frames
                        outsig(idx_mat(n,:)) = outsig(idx_mat(n,:)) + frames_syn(n,:);
                        wnt_sum(idx_mat(n,:)) = wnt_sum(idx_mat(n,:))+ winfunc;
%                         overadd_count(idx_mat(n,:)) = overadd_count(idx_mat(n,:)) + 1;
                    end
                    
                    % averaging the overlapped estimates for each sample
%                     outsig = outsig./overadd_count;
%                     wnt_sum = wnt_sum/overadd_count;
                    outsig = outsig./wnt_sum;
            end
        end

        function [mod_pha,theta_max,I_SNR] = modPhase(phase,alpha)
            % alpha range from 0.1-0.28, i.e., in steps of 0.02
            theta_max = alpha * pi;
            eta = 1/(tan(theta_max)^2);
            I_SNR = 10*log10(eta);
%             phi = -pi+(pi+pi)*rand(size(phase,1),size(phase,2));
%             %             mod_pha = bsxfun(@plus,phase,alpha*phi);
%             mod_pha = phase + alpha*phi;
            
            dist = alpha*unifrnd(-pi,pi,size(phase));
            mod_pha = phase + dist;
        end 
        
        outsig = overlap_add(frames,idx_mat);
        
        winfunc = modWin(Tw,Ts,varargin);
        
        [spec] = my_dct(frames,N,varargin);
        
        lcep = lifter(cep,L);

        delta = time_diff(coef, theta);
        
        function mel = f2mel(f)
            mel = 2595 * log10(1+f/700);
        end
        
        function f = mel2f(mel)
            f = (10.^(mel/2595)-1)*700;
        end
    end
end


