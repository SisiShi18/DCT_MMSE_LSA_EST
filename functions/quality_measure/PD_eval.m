function phase_dev_score = PD_eval(x_ref, x_est, x_noisy, param)
% ----------------------------------------------------------------------
%          Phase Deviation (PD) Objective Speech Quality Measure
%
%   This function implements the Phase Deviation distance measure.
%
%   Usage:  [PD_score PD_inst] = PD_eval(cleanFile, enhancedFile, noisyFile, param)
%           
%         cleanFile     - clean speech vector
%         enhancedFile  - enhanced speech vector
%         noisyFile     - noisy speech vector
%         param         - additional parameter-struct:
%                         * param.fs            - sampling frequnecy of speech vectors
%                         * param.wl1           - Window lenght of STFT in ms (default: 32ms)
%                         * param.step          - Window hop size of STFT in ms (default: 4ms)
%                         * param.plot_pd       - 1...plot PD, 0...do not plot (default: 0)
%                         * param.plot_pdmetric - 1...plot PD-metric, 0...do not plot (default: 0)
%
%         PD_score      - computed Phase Deviation distance measure
%         
%         Note that the PD measure is limited in the range [0, 4].
%
%  
% Copyright (c) 2014 by Andreas Gaich
% $Revision: 0.0 $  $Date: 03/11/2014 $

%% Init
param.fs_calc = 16000;
param.window_case='chebwin';
param.resolution = 2048;
if length(x_ref) ~= length(x_est) || length(x_ref) ~= length(x_noisy)
     error('Speech vectors do not have the same length')
end
if ~exist('param','var');
    error('Sampling frequency of speech signals missing')
end
if ~isfield(param,'wl1')
    param.wl1 = 0.032*param.fs_calc;
else
    param.wl1 = param.wl1/1000*param.fs_calc;
    param.wl1 = fix(param.wl1);
end
if ~isfield(param,'step')
    param.step = 0.004*param.fs_calc;
else
    param.step = param.step/1000*param.fs_calc;
    param.step = fix(param.step);
end
if ~isfield(param,'plot_pdmetric')
    param.plot_pdmetric = 0;
end
if ~isfield(param,'plot_pd')
    param.plot_pd = 0;
end
param.Frame_no = fix((length(x_ref)-param.wl1)/param.step) + 1;

%%
x_ref = x_ref(:);   x_est = x_est(:);   x_noisy = x_noisy(:);

x_ref = resample(x_ref,param.fs_calc,param.fs);     x_est = resample(x_est,param.fs_calc,param.fs);     x_noisy = resample(x_noisy,param.fs_calc,param.fs);

x_est = gain_normalization(x_est, x_ref);

[Xmat Xmat_disp] = sig2mat(x_ref,param);            
Xmat = Xmat(:,1:param.Frame_no);            
Xmat_disp = Xmat_disp(:,1:param.Frame_no);
[Xmat_est Xmat_est_disp] = sig2mat(x_est,param);    
Xmat_est = Xmat_est(:,1:param.Frame_no);    
Xmat_est_disp = Xmat_est_disp(:,1:param.Frame_no);
[Ymat Ymat_disp] = sig2mat(x_noisy,param);          
Ymat = Ymat(:,1:param.Frame_no);            
Ymat_disp = Ymat_disp(:,1:param.Frame_no); 
        
phi_x=angle(Xmat);              phi_x=phi_x(1:param.wl1/2+1,:);             phi_x=unwrap(phi_x,[],1);
phi_x_est=angle(Xmat_est);      phi_x_est=phi_x_est(1:param.wl1/2+1,:);     phi_x_est=unwrap(phi_x_est,[],1);
phi_y=angle(Ymat);              phi_y=phi_y(1:param.wl1/2+1,:);             phi_y=unwrap(phi_y,[],1);

phase_dev_diff = cos(phi_y-phi_x)-cos(phi_y-phi_x_est); % max = +-2;
phase_dev_score = mean(sum(phase_dev_diff.^2/(param.wl1/2+1)));     % fit PD_metric in the range of [0,4]

%% Plot

    phi_x=angle(Xmat_disp);              phi_x=phi_x(1:param.resolution/2+1,:);             phi_x=unwrap(phi_x,[],1);
    phi_x_est=angle(Xmat_est_disp);      phi_x_est=phi_x_est(1:param.resolution/2+1,:);     phi_x_est=unwrap(phi_x_est,[],1);
    phi_y=angle(Ymat_disp);              phi_y=phi_y(1:param.resolution/2+1,:);             phi_y=unwrap(phi_y,[],1);
    
    t = 0:param.step/param.fs_calc:(param.Frame_no-1)*param.step/param.fs_calc;
    f = 0:param.fs_calc/param.resolution:param.fs_calc/2;

if param.plot_pdmetric
  
    phase_dev_inst_disp = (cos(phi_y-phi_x)-cos(phi_y-phi_x_est)).^2;

    figure
    subplot(121)
    STFT_webpage(x_ref,param.fs_calc,'hamming');
    title('Spectrogram','FontSize',12);

    subplot(122)
    surf(t,f,phase_dev_inst_disp);
    axis([0 t(end) 0 f(end)]);
    view(2);
    shading flat;
    caxis([0 4]);
    xlabel('Time(s)','FontSize',12); ylabel('Frequency (Hz)','FontSize',12); title('Phase Deviation Metric','FontSize',12)
    set(gca,'FontSize',12);
end

if param.plot_pd
    
    phase_dev_inst_disp = 1-(cos(phi_y-phi_x));
    phase_dev_inst_est_disp = 1-(cos(phi_x_est-phi_x));

    figure
    subplot(131)
    STFT_webpage(x_ref,param.fs_calc,'hamming');
    title('Spectrogram','FontSize',12)
    
    subplot(132)
    surf(t,f,phase_dev_inst_disp);
    axis([0 t(end) 0 f(end)]);
    view(2);
    shading flat;
    caxis([0 2]);
%     colormap(hot);
    xlabel('Time(s)','FontSize',12); ylabel('Frequency (Hz)','FontSize',12); title('Phase Deviation (Noisy)','FontSize',12)
    set(gca,'FontSize',12);
    
    subplot(133)
    surf(t,f,phase_dev_inst_est_disp);
    axis([0 t(end) 0 f(end)]);
    view(2);
    shading flat;
    caxis([0 2]);
%     colormap(hot);
    xlabel('Time(s)','FontSize',12); ylabel('Frequency (Hz)','FontSize',12); title('Phase Deviation (Enhanced)','FontSize',12)
    set(gca,'FontSize',12);
end

function [sig_norm] = gain_normalization(sig_processed, sig_ref)

gain_factor = sqrt(var(sig_ref)/var(sig_processed));
sig_norm = gain_factor*sig_processed;

function [Amat Amat_disp] = sig2mat(asig,param)
wl1=param.wl1;
shift=param.step;
if strcmp(param.window_case,'hamming')
    ham1=(hamming(wl1));
elseif strcmp(param.window_case,'chebwin')
    ham1=chebwin(wl1,25);   % 25, 40
elseif strcmp(param.window_case,'rect')
    ham1=ones(wl1,1);
end

for ii = 1:2
    if ii == 1
        nfft2=wl1;
        Lx1=length(asig);
        frame_np=fix((Lx1-wl1)/shift) + 1;
        Amat=zeros(nfft2,frame_np+1);
        for frame_indx=1:frame_np; %100     %   5:round(length(x1)/step-2)
            seg=(frame_indx-1)*shift + 1:(frame_indx-1)*shift + 1 + wl1 - 1;
            ga=asig(seg);
            gaa=ga.*ham1;    ffga=fft(gaa,nfft2);  ffgaa=ffga(1:nfft2);
%             ffgaa=ffga(1:wl1);
            Amat(:,frame_indx)=ffgaa;
        end
    else
        nfft2=param.resolution;
        Lx1=length(asig);
        frame_np=fix((Lx1-wl1)/shift) + 1;
        Amat_disp=zeros(nfft2,frame_np+1);
        for frame_indx=1:frame_np; %100     %   5:round(length(x1)/step-2)
            seg=(frame_indx-1)*shift + 1:(frame_indx-1)*shift + 1 + wl1 - 1;
            ga=asig(seg);
            gaa=ga.*ham1;    ffga=fft(gaa,nfft2);  ffgaa=ffga(1:nfft2);
%             ffgaa=ffga(1:wl1);
            Amat_disp(:,frame_indx)=ffgaa;
        end
    end
end
