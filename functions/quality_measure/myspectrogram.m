function [handle] = myspectrogram(s, fs, T, w, nfft, Slim, alpha, cmap, cbar, type)

    %__________________________________________________________________________________________________________________
    % VALIDATE INPUTS, SET DEFAULTS
    switch nargin
    case 1, type='per'; cbar=false; cmap='default'; alpha=false; Slim=[-59,-1]; w=@hamming; T=[32 32/8]; nfft=1024; fs=8000;
    case 2, type='per'; cbar=false; cmap='default'; alpha=false; Slim=[-59,-1]; w=@hamming; T=[32 32/8]; nfft=1024; 
    case 3, type='per'; cbar=false; cmap='default'; alpha=false; Slim=[-59,-1]; w=@hamming; T=[32 32/8]; 
    case 4, type='per'; cbar=true; cmap='default'; alpha=false; Slim=[-59,-1]; w=@hamming; 
    case 5, type='per'; cbar=true; cmap='default'; alpha=false; Slim=[-59,-1]; 
    case 6, type='per'; cbar=true; cmap='default'; alpha=false;
    case 7, type='per'; cbar=true; cmap='default'; 
    case 8, type='per'; cbar=true; 
    case 9, type='per'; 
    case 10
    otherwise, error('Invalid number of input arguments.');
    end
    %__________________________________________________________________________________________________________________
    % DECLARE VARIABLES
    if(isstr(s)), [s, fs] = wavread(s); end; % read audio data from file
    Tw = T(1);                          % frame width (ms)
    Ts = T(2);                          % frame shift (ms)
    Nw = round(fs*Tw*0.001);            % frame width (samples)
    nfft = 1024;
    Ns = round(fs*Ts*0.001);            % frame shift (samples)
    N  = length(s);                     % length of speech signal (samples)
    Smin = Slim(1);                     % lower normalized dynamic range limit 
    Smax = Slim(2);                     % upper normalized dynamic range limit
    if(isstr(w)), w = str2func(w); end; % obtain window function handle from string input
    %__________________________________________________________________________________________________________________
    % SPEECH PREPROCESSING
    if(islogical(alpha) && alpha), s = filter([1 -0.95],1,s);   % apply a typical preemphasis filter
    elseif(~islogical(alpha)) s = filter(alpha,1,s); end;       % apply custom preemphasis filter

    %__________________________________________________________________________________________________________________
    % GET SPECTROGRAM DATA 
    %[S,F,T] = spectrogram(s,w(Nw).',Nw-Ns,nfft,fs);            % MATLAB's new spectrogram function
    %[S,F,T] = specgram(s,nfft,fs,w(Nw).',Nw-Ns);               % MATLAB's depreciated spectrogram function
    [S,F,T] = toframes(s,w,T,fs,nfft,type);                     % Framing function, use this if you do not have the Signal Processing Toolbox

    %__________________________________________________________________________________________________________________
    % SET DYNAMIC RANGE 
    S = abs(S);                         % compute magnitude spectrum 
    S = S/max(max(S));                  % normalize magntide spectrum
%    S(S<eps) = eps;                     % prevent zero spectral magnitude values 
    S = 20*log10(S);                    % compute power spectrum in dB
    
    %__________________________________________________________________________________________________________________
    % PLOT AND LABEL RESULTS 
    handle = imagesc(T, F, S, [Smin Smax]);
    %handle = imagesc(T, F, S, 'CDataMapping', 'direct');
    axis('xy');
    axis([0 N/fs  0 fs/2]);
%    xlabel('time (s)', 'FontSize', 8, 'FontWeight', 'n');
%    ylabel('frequency (Hz)', 'FontSize', 8, 'FontWeight', 'n');
%    set(gca,'YDir','normal', 'FontSize', 6);

    if(cbar), colorbar('FontSize',6); end

    %__________________________________________________________________________________________________________________
    % DEFINE CUSTOM COLOR MAPS
    switch(lower(cmap))
    case {'default'}
        colormap('gray');
%         colormap(jet);
        map=colormap;
        colormap(1-map);    
    otherwise, colormap(cmap);
    end
    function [S,F,T] = toframes(s,w,T,fs,nfft,type)

    %__________________________________________________________________________________________________________________
    % VALIDATE INPUTS 
    switch nargin
    case 1, type='per'; nfft=1024; fs=8000; T=[32 32/8]; w={@hamming, @hamming};
    case 2, type='per'; nfft=1024; fs=8000; T=[32 32/8]; 
    case 3, type='per'; nfft=1024; fs=8000; 
    case 4, type='per'; nfft=1024; 
    case 5, type='per';
    case 6
    otherwise, error('Invalid number of input arguments.');
    end

    %__________________________________________________________________________________________________________________
    % DEFINE VARIABLES
    if(isstr(s)) [s, fs, nbits] = wavread(s); else, nbits=16; end;

    s = s(:).';
%    s = s-mean(s);
    smax = max(abs(s))/0.999;
    s = s/smax; % normalize

    Tw = T(1);                                              % frame length [ms]
    Ts = T(2);                                              % frame frameshift [ms]
    N  = round(fs*Tw*0.001);                                % frame length [samples]
    Z  = round(fs*Ts*0.001);                                % frame shift [samples]
    D  = mod(length(s), Z);                                 % add N-D zeros to the end
    G  = (ceil(N/Z)-1)*Z;                                   % add G zeros to the beginning
    s  = [zeros(1,G) s zeros(1,N-D)];                       % zero pad signal to allow an integer number of segments
    ss = length(s);                                         % length of the signal for processing
    t  = [0:ss-1]/fs;                                       % time vector
    M  = ((ss-N)/Z)+1;                                      % number of overlapping segments
    if(isstr(w)), w=str2func(w); end;
    wa = w(N).';                                            % analysis window A (for magnitude component)
    wsyn = 0.5-0.5*cos((2*pi*((0:N-1)+0.5))/N);             % synthesis window


    %__________________________________________________________________________________________________________________
    % SPLIT SPEECH INTO OVERLAPPED FRAMES (EACH ROW IS A FRAME), AND WINDOW THE FRAMES
    indf = Z*(0:(M-1)).';                                   % indexes for frames
    inds = (1:N);                                           % indexes for samples
    refs = indf(:,ones(1,N)) + inds(ones(M,1),:);           % sample indexes for each frame

    segments_s = s(refs);                                   % split into overlapped frames (using indexing)
    segments_sm = segments_s .* wa(ones(M,1),:);            % apply magnitude spectrum analysis window 


    %__________________________________________________________________________________________________________________
    % PERFORM COMPLEX FOURIER SPECTRUM ANALYSIS 
    F = [0:nfft-1]/(nfft-1)*fs;
    T = [0:M-1]/(M-1)*ss/fs;

    switch(lower(type))
    case {'per','periodogram'}
        S = fft(segments_sm, nfft, 2);                      % short-time Fourier analysis
        S = abs(S).^2/N;                                    % periodogram PSD estimates
        S = sqrt(S);                                        % magnitude spectrum (for consistency)
    case {'lp','lpc','ar'}
        p = 12;                                             % order of AR model
        [A,E] = lpc(segments_sm.',p);
        S = repmat(sqrt(E),1,nfft)./abs(fft(A,nfft,2));     % LP-based (AR) magnitude spectrum
    otherwise
    end
    
    S = S.';
    F = F.';
    T = T;
    