function noise_psd_init = init_noise(noisy_mag,frame_len,shift_len)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%This m-file computes an initial noise PSD estimate by means of a
%%%%Bartlett estimate.
%%%%Input parameters:   noisy_mag:      noisy magnitude spectrum
%%%%                    frame_len:      frame size
%%%%                    shift_len:      shift size of frame
%%%%Output parameters:  noise_psd_init: initial noise PSD estimate
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch frame_len
    case 20
        if shift_len == 10
            noise_psd_init = (noisy_mag(1:6,:)).^2;
        elseif shift_len == 5
            noise_psd_init = (noisy_mag(1:12,:)).^2;
        end
    case 16
        if shift_len == 4
            noise_psd_init = (noisy_mag(1:20,:)).^2;
        elseif shift_len == 8
            noise_psd_init = (noisy_mag(1:10,:)).^2;
        end
    case 32
        if shift_len == 8
            noise_psd_init = (noisy_mag(1:10,:)).^2;
        elseif shift_len == 16
            noise_psd_init = (noisy_mag(1:5,:)).^2;
        end
    case 40
        noise_psd_init = (noisy_mag(1:3,:)).^2;
end

noise_psd_init = mean(noise_psd_init);
end