function fbanks = filter_banks(obj)
    mel_edges = linspace(0,obj.f2mel(obj.Fs/2),obj.N_fbanks+2);
    f_edges = obj.mel2f(mel_edges);
    idx_edges = (obj.N_analysis/2-1)*f_edges/(obj.Fs/2) + 1;
    fbanks = zeros(obj.N_fbanks,obj.N_analysis/2+1);

    for k = 1:obj.N_fbanks
        first_bin = ceil(idx_edges(k));
        last_bin  = floor(idx_edges(k+2));
        fbanks(k, first_bin:last_bin) = interp1(idx_edges(k:(k+2)), [0 1 0], first_bin:last_bin);
        
        if (obj.norm_fbanks == 1)
            fbanks(k,:) = fbanks(k,:)/sum(fbanks(k,:)); % normalized.
        end
    end
end