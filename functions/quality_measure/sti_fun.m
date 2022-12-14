function [STI] = sti_fun(clean, noisy,fs,method)
% @inputs       clean  - reference speech vector
%               noisy  - modified speech vector
%               fs     - sampling frequency (Hz)
%               method - sti method (i.e. 'ansi' or 'payton' or 'drullman')
%
% @outputs      STI    - sti score
%
% @usage        score = sti(clean, noisy, 16000, 'drullman');

    wav     = clean(:);
    wav_deg = noisy(:);

    if (fs==16000)
        nbands=7;
        W = [0.13 0.14 0.11 0.12 0.19 0.17 0.14];
    elseif (fs==8000)
        nbands=6;
        W = [0.13 0.14 0.11 0.12 0.19 0.14];
        W = W/sum(W);
    else
        wav = resample(wav, 8000, fs);
        wav_deg = resample(wav_deg, 8000, fs);
        fs = 8000;
        nbands=6;
        W = [0.13 0.14 0.11 0.12 0.19 0.14];
        W = W/sum(W);
    end

    nmods=16;
    mtf = zeros(nmods,nbands);
    N = min(length(wav),length(wav_deg)); % the length of speech
    wav = wav(1:N);
    wav_deg = wav_deg(1:N);
    delta_f = fs/N; % frequency resolustion

    %the cutoff frequencies of octave-bands with last center frequency at fs/2.
    centre_freqs = 0.5*fs./power(2,fliplr([0:1:nbands-1])); % the center frequencies of octave-bands
    bwfact=2^(1/2);
    trans_freqs = centre_freqs/bwfact; %the cutoff frequencies of octave-band filters

    %the cutoff frequencies of 1/3-octave-bands with last center frequency at 16 Hz.
    fc=16./power(2,fliplr([0:1:nmods-1])/3);
    bwfact=2^(1/6);
    fl=fc/bwfact;
    fh=fc*bwfact;
    indl = round(fl/delta_f+1);
    indh = round(fh/delta_f+1);

    for k = 1:nbands
        % 1) The speech was band-pass filtered by seven octave-band filters.
        if (k == nbands)
            [b, a] = butter(6, trans_freqs(k)/(fs/2), 'high');
        else
            [b, a] = butter(3, [trans_freqs(k) trans_freqs(k+1)]/(fs/2)); 
        end
        bp = filter(b, a, wav);
        bp_deg=filter(b,a,wav_deg);
        % 2) Squared
        bp = bp.^2;
        bp_deg=bp_deg.^2;
        % 3) Low-pass filter to derive the energy envelop and normalize to unit mean.
        [b, a] = butter(4, 50/(fs/2), 'low');
        r = filter(b, a, bp); r=r/mean(r);
        r_deg = filter(b, a, bp_deg); r_deg=r_deg/mean(r_deg);
        % 4) get the energy spectrum
        r_fft = fft(r)/length(r);
        r_fft_deg = fft(r_deg)/length(r_deg);
        % 5) 1/3-octave-band analysis: summing the components over
        % 1/3-octave-band frequencies.
        switch(lower(method))
        case {'ansi'}, numer=abs(r_fft_deg).^2;                     % original STI (ANSI Standard) 
        case {'payton'}, numer=abs(r_fft.*conj(r_fft_deg));  % Payton & Braida, JASA, 106, 3637-3648, 1999.
        case {'drullman'}, numer=real(r_fft.*conj(r_fft_deg));      % Drullman et al., JASA, 95, 2670-2680, 1994.
        otherwise, error('STI type not supported.');
        end
        denom=abs(r_fft).^2;
        for jj = 1:nmods
            m_numer = sum(numer(indl(jj):1:indh(jj)));
            m_denom = sum(denom(indl(jj):1:indh(jj)));
            switch (lower(method))
                 case {'ansi'}, mtf(jj,k)=sqrt(m_numer/m_denom);
                 otherwise, mtf(jj,k)=m_numer/m_denom;
            end
        end
    end
%max(mtf)
    mtf(find(mtf >= 1.0)) = 0.99;
    mtf(find(mtf <= 0.0)) = 0.00001;
    SNR = 10 * log10(mtf ./ (1-mtf));
    SNR(find(SNR <= -15)) = -15;
    SNR(find(SNR >= 15)) = 15;
    SNR_mean1 = mean(SNR);
    TI = (SNR_mean1+15)/30;
    STI = sum(TI.*W);
end
