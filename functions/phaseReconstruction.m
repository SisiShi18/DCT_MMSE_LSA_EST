function [yhat,phaseRecon] = phaseReconstruction(x,fs,freqMax,estimationMethod)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstructs the STFT-phase of voiced speech sounds given the noisy 
% time-domain signal x according to:
%
% [1] Martin Krawczyk, Timo Gerkmann, "STFT Phase Reconstruction in Voiced Speech
% for an Improved Single-Channel Speech Enhancement", Transaction on Audio,
% Speech, and Language Processing, UNDER REVIEW
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Attention: For this code to run, you need to provide an estimate of the 
% fundamental frequency "f0" for each STFT signal segment together with a 
% binary decision if the segment contains voiced speech (1) or not (0) in 
% "voicedYN". For more instructions you may also have a look at line 108. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:  
%         x:
%             - noisy time domain input
%         fs:
%             - sampling frequency in [Hz]
%         freqMax (optional):
%             - frequency in [Hz] up to which the phase should be reconstructed
%         estimationMethod (optional):
%             - if == "freqOnly", then the phase is only
%               reconstructed along frequency between the spectral harmonics. The
%               phase of bands that directly contain harmonic components is not reconstructed
%               (instead the noisy phase is used). (DEFAULT)
%             - if == "combi", then the phase is estimated along time on the
%               harmonics and along frequency between the harmonics as outlined
%               in [1, Fig.4]
%
% Output: 
%         phaseRecon: 
%             - matrix containing the reconstructed STFT-phase for voiced speech. 
%               During unvoiced speech, "phaseRecon" contains the noisy phase.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Martin Krawczyk, martin.krawczyk@uni-oldenburg.de
% Timo Gerkmann, timo.gerkmann@uni-oldenburg.de
% Speech Signal Processing Group
% Department of Medical Physics and Acoustics and Cluster of Excellence "Hearing4all"
% University of Oldenburg, Germany
% (c) 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set estimation method (if not provided):
if ~exist('estimationMethod','var') || isempty(estimationMethod)
    estimationMethod = 'freqOnly';
end
switch estimationMethod
    case 'combi'
        phaseEstimationAlongTimeYN = 1; % Reconstruct the phase along time in bands which directly contain harmonics
        phaseEstimationAlongFreqYN = 1; % AND between the harmonics across frequencies
    case 'freqOnly'
        phaseEstimationAlongTimeYN = 0; % Only reconstruct the phase across frequencies between harmonics
        phaseEstimationAlongFreqYN = 1;
    case 'timeOnly'
        phaseEstimationAlongTimeYN = 1; % Only reconstruct the phase across frequencies between harmonics        
        phaseEstimationAlongFreqYN = 0;
    otherwise
        error(['Unknown value for "estimationMethod": ' estimationMethod]);
end

% Setup:
fLengthTime = (32e-3);                        % frame length in seconds
fShiftTime  = fLengthTime/8;                  % frame shift in seconds
fLength     = fLengthTime*fs;                 % frame length in samples
fShift      = fShiftTime*fs;                  % frame shift in samples
fOverlap    = fLength - fShift;               % frame overlap in samples
zeroPadFact = 1;                              % amount of zero padding (fftLength/fLength == zeroPadFact)
fftLength   = fLength*zeroPadFact;            % fft length (zero-padding if fftLength>fLength)
nBins       = fftLength/2+1;                  % number of relevant FFT bins
fftBA       = fs/fftLength;                   % distance between center-frequencies of two adjacent STFT-bands in [Hz]

% Spectral analysis window function:
% anaWin      = sqrt(hann(fLength,'periodic')); % spectral analysis window
% anaWin      = ones(fLength,1);
% anaWin      = hamming(fLength,'periodic');  
anaWin      = hann(fLength,'periodic');

% Set maximum frequency up to which we reconstruct the phase (if not provided):
if ~exist('freqMax','var') || isempty(freqMax)
    freqMax     = min(fs/2-fftBA,4000-fftBA);
end

% Check which window is used:
if max(abs(anaWin-hann(fLength,'periodic')))<1e-13
    windowType = 'hann'; a = 0.5;
elseif max(abs(anaWin-ones(fLength,1)))<1e-13
    windowType = 'rect'; a = 1;
elseif max(abs(anaWin-sqrt(hann(fLength,'periodic'))))<1e-13
    windowType = 'sqrtHann'; a = 0.5; fftBANoZP = fs/fLength;  
elseif max(abs(anaWin-hamming(fLength,'periodic')))<1e-13
    windowType = 'hamming'; a = 0.54;
else
    windowType = 'other';
end

% Compute STFT of x:
[X,fax,tax] = spectrogram(x,anaWin,fOverlap,fftLength,fs);
nFrames     = length(tax);                     % number of frames  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% You may put your fundamental frequency estimator and voiced/unvoiced classificator here: (or within the for loop if you like)
addpath('voicebox')
[f0,t,pv]=fxpefac(x,fs,fShift/fs);
f0       = interp1(t,f0,tax,'linear','extrap'); % Estimates must be in [Hz]
f0(f0<50)=50;f0(f0>500)=500;
voicedYN = interp1(t,pv,tax,'linear','extrap'); % ==0: unvoiced / ==1: voiced
voicedYN(voicedYN<0.5)=0;
voicedYN(voicedYN>=0.5)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialization:
phaseRecon = zeros(nBins,nFrames);            % Reconstructed phase
if strcmp(windowType,'other')                 % Initialize look-up table if necessary
    [~] = getWindowPhaseViaZP(anaWin,[],[],fs,1); 
end

% Framewise processing:
for frameIdx = 1:nFrames
    
    % Set the enhanced phase to the noisy phase as a starting point for enhancement:
    phaseRecon(:,frameIdx) = angle(X(:,frameIdx));
    
    % Check if current frame is a voiced speech onset:
    isOnsetYN =  voicedYN(frameIdx) > voicedYN(max(1,frameIdx-1));
    
    % Fundamental frequency dependent bandwidth:
    deltaFreq = f0(frameIdx)/2;
    deltaBin  = floor(deltaFreq/fftBA) + 1;
    
    if voicedYN(frameIdx) && (frameIdx>1) %&& ~isOnsetYN
        
        % Compute harmonic frequencies and corresponding frequency bins:
        nHarmonics  = floor(freqMax  / f0(frameIdx));  % Number of harmonics
        frequencies = [1:nHarmonics].' * f0(frameIdx); % Frequencies of the harmonics
        freqBin     = 1 + round(frequencies/fftBA);    % STFT-bands containing the harmonic frequencies
        
        % a) Phase reconstruction along TIME on the harmonics (if activated):
        if ~isOnsetYN && (phaseEstimationAlongTimeYN~=0)
            phaseRecon(freqBin,frameIdx) = phaseRecon(freqBin,frameIdx-1) + 2.*pi.*frequencies.*fShiftTime;
        end
        
        % b) Phase reconstruction along FREQUENCY:
        if phaseEstimationAlongFreqYN
        if any(strcmp(windowType,({'rect','hamming','hann'})))
            phiW0 = getAnalyticWindowPhase(fax(freqBin),frequencies,fLength,fs,a);
        elseif strcmp(windowType,'other')
            phiW0 = getWindowPhaseViaZP(anaWin,fax(freqBin),frequencies,fs,0);
        end
        end
        if f0(frameIdx)>fftBA % Exclude extremely low fundamental frequency estimates    
            % Reconstruct the phase along FREQUENCY starting from the phase of
            % the bands which contain harmonic components [1,eq. (13)]: 
            for binCntr = -deltaBin:+deltaBin
                binIdx = max(2,min(freqBin+binCntr,nBins-1)); % Omit the first and the last frequency index
                if (binCntr~=0)&&(~any(binIdx==freqBin))   
                    switch windowType
                        case {'rect','hamming','hann'} 
                            phiW = getAnalyticWindowPhase(fax(binIdx),frequencies,fLength,fs,a);
                            if phaseEstimationAlongFreqYN
                            phaseRecon(binIdx,frameIdx) = phaseRecon(freqBin,frameIdx) - phiW0 + phiW; % [1, eq(13)]
                            end
                        case 'sqrtHann' % Specialized implementation for sqrtHann (faster than zero-padding):
                            phaseRecon(binIdx,frameIdx) = phaseRecon(freqBin,frameIdx) - binCntr*pi/zeroPadFact;
                            myI = abs(frequencies-fax(binIdx))>=fftBANoZP*1.5;
                            if any(myI==true)
                                phaseRecon(binIdx(myI),frameIdx) = phaseRecon(freqBin(myI),frameIdx) - binCntr*pi/zeroPadFact + pi*floor((abs(frequencies(myI)-fax(binIdx(myI)))-fftBANoZP/2)./fftBANoZP);
                            end
                        case {'other'}
                            phiW = getWindowPhaseViaZP(anaWin,fax(binIdx),frequencies,fs,0);
                            phaseRecon(binIdx,frameIdx) = phaseRecon(freqBin,frameIdx) - phiW0 + phiW; % [1, eq(13)]
                    end
                end
            end
            % Take care of lowest frequency bands:
            if ((freqBin(1)-deltaBin) >2) && (freqMax==fs/2-fftBA)
                for toFixIdx = freqBin(1)-deltaBin-1:-1:2
                    binCntr = toFixIdx-freqBin(1);
                    switch windowType
                        case {'rect','hamming','hann'} 
                            phiW(1) = getAnalyticWindowPhase(fax(toFixIdx),frequencies(1),fLength,fs,a);
                            if phaseEstimationAlongFreqYN
                            phaseRecon(toFixIdx,frameIdx) = phaseRecon(freqBin(1),frameIdx) - phiW0(1) + phiW(1);                
                            end
                        case 'sqrtHann' % Specialized implementation for sqrtHann (faster than zero-padding):
                            phaseRecon(toFixIdx,frameIdx) = phaseRecon(freqBin(1),frameIdx) - binCntr*pi/zeroPadFact;
                            if abs(frequencies(1)-fax(toFixIdx))>=fftBANoZP*1.5
                                phaseRecon(toFixIdx,frameIdx) = phaseRecon(freqBin(1),frameIdx) - binCntr*pi/zeroPadFact + pi*floor((abs(frequencies(1)-fax(toFixIdx))-fftBANoZP/2)./fftBANoZP);       
                            end
                        case 'other'
                             phiW(1) = getWindowPhaseViaZP(anaWin,fax(toFixIdx),frequencies(1),fs,0);
                             phaseRecon(toFixIdx,frameIdx) = phaseRecon(freqBin(1),frameIdx) - phiW0(1) + phiW(1);
                    end
                end
            end
            % Take care of highest frequency bands:
            if ((freqBin(end)+deltaBin) < nBins-1) && (freqMax==fs/2-fftBA)
                for toFixIdx = freqBin(end)+deltaBin+1:nBins-1
                    binCntr = toFixIdx-freqBin(end);
                    switch windowType
                        case {'rect','hamming','hann'} 
                            phiW(end) = getAnalyticWindowPhase(fax(toFixIdx),frequencies(end),fLength,fs,a);
                            if phaseEstimationAlongFreqYN
                            phaseRecon(toFixIdx,frameIdx) = phaseRecon(freqBin(end),frameIdx) - phiW0(end) + phiW(end);
                            end
                        case 'sqrtHann' % Specialized implementation for sqrtHann (faster than zero-padding):
                            phaseRecon(toFixIdx,frameIdx) = phaseRecon(freqBin(end),frameIdx) - binCntr*pi/zeroPadFact;
                            if abs(frequencies(end)-fax(toFixIdx))>=fftBANoZP*1.5
                                phaseRecon(toFixIdx,frameIdx) = phaseRecon(freqBin(end),frameIdx) - binCntr*pi/zeroPadFact + pi*floor((abs(frequencies(end)-fax(toFixIdx))-fftBANoZP/2)./fftBANoZP);
                            end
                        case 'other'
                            phiW(end) = getWindowPhaseViaZP(anaWin,fax(toFixIdx),frequencies(end),fs,0);
                            phaseRecon(toFixIdx,frameIdx) = phaseRecon(freqBin(end),frameIdx) - phiW0(end) + phiW(end);
                    end
                end
            end
        end      
    end % end of "isVoiced"
   
end % end of frameIdx

% Wrap the phase:
phaseRecon = angle(exp(1i*phaseRecon));
yhat=overlapadd(irfft((abs(X).*exp(1i*phaseRecon)))',anaWin,fShift);
end % end of main function



function phiW = getAnalyticWindowPhase(fax,frequencies,M,fs,a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the spectral phase of symmetric half-cosine-based window functions
% like the rectangular, Hamming, and Hann windows using [1, eq. (16)].
% Input:  fax:         Frequency [Hz] where the phase should be estimated
%         frequencies: Frequency [Hz] of the harmonic (around which the
%                      spectral representation of the window is centered).
%         M:           length of the window in samples
%         fs:          sampling rate [Hz]
%         a:           shape coefficient as defined in [1, eq. (15)]
% Output: phiW:        phase at frequency "fax" of the spectral representation of the window 
%                      centered around "frequencies"
%
% Please note that "fax" and "frequencies" can also be vectors of the same size.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Omega = 2*pi*(fax-frequencies)/fs;
    if a~=1
        phiW = angle(sin(M/2*Omega) .* exp(-1i*(M-1)/2*Omega) .* (a./sin(0.5*Omega) - (1-a)/2 * ( exp(-1i*pi/M)./sin(0.5*(Omega - 2*pi/M)) + exp(1i*pi/M)./sin(0.5*(Omega + 2*pi/M)) ) ));
    else
        phiW = angle(sin(M/2*Omega) .* exp(-1i*(M-1)/2*Omega) .* (a./sin(0.5*Omega) ));
    end
    % Handle NaNs at "0/0":
    if any(isnan(phiW))
        phiW(sin(0.5*Omega)==0) = 0;
        if a~=1
            phiW(sin(0.5*(Omega - 2*pi/M)) == 0) = -pi;
            phiW(sin(0.5*(Omega + 2*pi/M)) == 0) = +pi;
        end
    end
end



function phiW = getWindowPhaseViaZP(anaWin,fax,frequencies,fs,initYN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the spectral phase of arbitrary windows by using a large amount 
% of zero padding, achieving a quasi-continous sampling.
% Input:  anaWin:      column vector containing the time-domain window function
%         fax:         Frequency [Hz] where the phase should be estimated
%         frequencies: Frequency [Hz] of the harmonic (around which the
%                      spectral representation of the window is centered).
%         fs:          sampling rate [Hz]
%         initYN:      if initYN==1(true), then zero-padding is performed
%                      and a large array containing all phase values is
%                      created (persistent look-up table). if initYN==0(false),
%                      then the previously created look-up table is used directly.
% Output: phiW:        phase at frequency "fax" of the spectral representation of the window 
%                      centered around "frequencies" 
%
% Please note that "fax" and "frequencies" can also be vectors of the same size.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    persistent anaPhase; persistent anaPhaseFax;
    
    if initYN ==1
        anaPhase    = angle(fft(anaWin,2^16));
        anaPhaseFax = (0:length(anaPhase)-1)./length(anaPhase) .*fs;
        anaPhase    = anaPhase((1:2^15)+1);anaPhaseFax = anaPhaseFax((1:2^15)+1);
        phiW        = [];
    else
        nHarmonics = length(frequencies);
        phiW       = zeros(nHarmonics,1);
        for h = 1:nHarmonics
            [~,myidx] = min(abs( anaPhaseFax-(abs(frequencies(h)-fax(h))) ));
            phiW(h)   = -sign(frequencies(h)-fax(h)).*anaPhase(myidx);
        end
    end
end
