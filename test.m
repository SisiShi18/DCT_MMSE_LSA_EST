%% Script Description:
%   Implementaion of the DCT Short-Time Log-Spectral Amplitude MMSE estimators (CLSA);
%   Output enhanced speech files in the audios folder.
%   Output plots of spectrograms in the figs folder.
% %--------------------------------------------------------------------------
clear all; 
close all; 
clc;
%--------------------------------------------------------------------------
%%                           Define Globals
%--------------------------------------------------------------------------
global min_db max_db Cl Cg FLSA_sg CLSA_l_gain_double CLSA_n_gain CLSA_l_gain  
global include_clean acc acc_ideal
%--------------------------------------------------------------------------
min_db =-10; max_db = 35;
include_clean= 0;
method_start = include_clean + 2;
%--------------------------------------------------------------------------
%%                          Progress indication
%--------------------------------------------------------------------------
progress = {'-','\','\','|','/'}; P = length(progress); pp = 0;
%--------------------------------------------------------------------------
%%                 Define working directory and file path
%--------------------------------------------------------------------------
fprintf('Setting working directory and file paths...\n')
currentFolder = pwd;
addpath(genpath(currentFolder));
exp_dir = [currentFolder filesep];
aud_dir = [exp_dir 'audios' filesep];
fig_dir = [exp_dir 'figs' filesep];
LUT_dir = [exp_dir 'gain_LUT' filesep];
%--------------------------------------------------------------------------
%%                    Create clean and noise file paths
%--------------------------------------------------------------------------
clean_files = dir([aud_dir 'clean' filesep '*.wav']); % path to clean speech test files.
% %% ---------------------------Clean file name------------------------------
fidx = 1; % list 7, '{'FB07_09.wav'}','The dune rose from the edge of the water.'
% fidx = 2; % list32, 'FF32_09.wav', 'The purple tie was ten years old.'
%% ----------------------------Noise type----------------------------------
% 2  : 'Pink noise';
% 3  : 'white';
% 6  : 'Speech noise';
% 19 : 'Voice Babble';
% 20 : 'F-16 two-seat';
% 21 : 'Car Factory electrical welding';
% 22 : 'Car Factory production hall';
% 23 : 'Car Volvo-340 asphalt road';
%--------------------------------------------------------------------------
% noise_files = dir([aud_dir 'noise_babble' filesep '*.wav']); % path to noise test files.
% noise_type = [19];
% noise_files = dir([aud_dir 'noise_F_16' filesep '*.wav']); % path to noise test files.
% noise_type = [20];
noise_files = dir([aud_dir 'noise_speech' filesep '*.wav']); % path to noise test files.
noise_type = [6];
% noise_files = dir([aud_dir 'noise_white' filesep '*.wav']); % path to noise test files. 
% noise_type = [3];
% noise_files = dir([aud_dir 'noise_pink' filesep '*.wav']); % path to noise test files. 
% noise_type = [2];
% noise_files = dir([aud_dir 'noise_car_Volvo340_asphalt_road' filesep '*.wav']); % path to noise test files. 
% noise_type = [23];
%--------------------------------------------------------------------------
%%                  Define SNR range 
%--------------------------------------------------------------------------
% plt.SNR_arr = 0:5:15; % SNR range
plt.SNR_arr = 5; % SNR range
x = plt.SNR_arr;
%--------------------------------------------------------------------------
%%                       Load Look Up Tables (LUT)
%--------------------------------------------------------------------------
% don't have to load every time; comment out 'clear all;' and this section
%--------------------------------------------------------------------------
fprintf('Loading Gain Look Up Table...\n');
savefile = strcat(LUT_dir,'CSA_gain_LUT.mat'); load(savefile);
savefile = strcat(LUT_dir,'LSAGain_40.mat'); load(savefile);
%--------------------------------------------------------------------------
%%                      Gain function Interpolation
%--------------------------------------------------------------------------
fprintf('Interpolating gain functions...\n');
nu=0.3;
Rprior= -40:1:40; % in dB
Rpost=-40:1:40; % in dB
[FLSA_sg]=TabulateGainGamma2logmmse(Rprior,Rpost,nu);
%======================== Double precision ================================
[CLSA_l_gain_double] = Tabulate_CLSA_l_gain(Rprior,Rpost);
%--------------------------------------------------------------------------
Cl = griddedInterpolant({xi_CSA,(gammaK_CSA')},(CSA_l_gain),'spline','nearest');
% Cg = griddedInterpolant({xi_CSA,(gammaK_CSA')},(CSA_g_gain),'spline','nearest');
% --------------------------------------------------------------------------
%%          Define MMSE estimator types/parameters & functions
%--------------------------------------------------------------------------
% CLSA/FLSA MMSE:
%   short-time (Cosine/Fourier) Log Spectral Amplitude MMSE estimator
% CSA MMSE:
%   short-time Cosine Spectral Amplitude MMSE estimator
% beta-order MMSE + PC: Peceptually Weighted beta-order spectual amplitude
% estimator + Phase Compensation technique
% PE: DCT Polarity Spectrum Estimator
%--------------------------------------------------------------------------
MMSE_types = {'CLSA MMSE','FLSA MMSE','CLSA MMSE double','CSA MMSE','beta-order MMSE+PC',...
    'CLSA MMSE + PE','CLSA MMSE + PE ideal'};
%--------------------------------------------------------------------------
%           (1)        (2)              (3)            (4)           
myfun = {@CLSA_MMSE, @FLSA_MMSE, @CLSA_MMSE_double, @CSA_MMSE_t, ...
                @PC_Beta_MMSE, @CLSA_MMSE_PE, @CLSA_MMSE_PE_ideal} ;
%                      (5)           (6)                (7)          
%--------------------------------------------------------------------------
%%                plot spectrum / output audios
%--------------------------------------------------------------------------
plot_spect = true;
output_audios=1;
%--------------------------------------------------------------------------
%%                              Methods
%--------------------------------------------------------------------------
if include_clean
    methods = {'clean','Unprocessed',...
    'FLSA normal', 'FLSA super-gaussian','beta-order+PC',...
    'CSA laplace',...
    'CLSA normal','CLSA laplace',...
    'CLSA normal + PE','CLSA laplace + PE',...
    'CLSA normal + PE ideal','CLSA laplace + PE ideal'};
     para_idx = [[0,0,0];[0,0,0];...
        [2,1,2];[2,3,2];[5,1,5];...
        [4,2,4];...
        [1,1,1];[1,2,1];...
        [6,1,6];[6,2,6];...
        [7,1,8];[8,2,7];];
%--------------------------------------------------------------------------
else
    methods = {'Unprocessed',...
               'EM-LSA','SG-LSA','$\beta$-PC',...
               'L-STSA',...
               'N-LSA','L-LSA',...
               'N-PoE','L-PoE',...
               'N-PoE-O','L-PoE-O'};
    para_idx = [[0,0,0];...
         [2,1,2];[2,3,2];[5,1,5];...
         [4,2,4];...
         [1,1,1];[1,2,1];...
         [6,1,6];[6,2,6];...
         [7,1,7];[7,2,7]];
end
%--------------------------------------------------------------------------
mmse_inputs.xi_min = 10^(-25/10);
mmse_inputs.gain_min = 0.1;
%--------------------------------------------------------------------------
%%             Define speech/noise variance estimation method
% --------------------------------------------------------------------------
% xi_type = {'noise_ideal'}; %ideal noise psd estimate
xi_type = {'DD'}; % Decision-direct a priori SNR estimation
%--------------------------------------------------------------------------
%%                 Define speech quality measure types
%--------------------------------------------------------------------------
measure = {'pesq',...% Quality Metric,
    'DP',... % Phase-Aware Metric
    'stoi',... % Intelligibility Metric
    'seg'};% Quality Metrics
m_idx = [1,2,3]; % measure index
%--------------------------------------------------------------------------
%%              Set signal analysis/synthesis parameters
%--------------------------------------------------------------------------
Fs = 16e3;
synWinType = {'hamming'};
speech_priori = {'norm','lap','gamma'};
frame_len = [32]; % duration in time (ms)
shift_len = [8];  % 75% over-lapping
Tw_arr= round(Fs.*frame_len./1000);
%--------------------------------------------------------------------------
%%             Set the parameters for save the results
%--------------------------------------------------------------------------
plt.methods = methods;
plt.pesq_mode ='wb';
plt.save_plots = false; % save the plots
plt.save_data = false; % save the outputs
plt.plotType = 'bar';
plt.save_dir = fig_dir;

% plt.xTicks = {'0dB','5dB','10dB','15dB'};
plt.xTicks = {'5 dB'};

plt.unit = 'inch';
plt.width= 7;
plt.height= plt.width/((1+sqrt(5))/2);
plt.fontType = 'Times';
plt.fontSize = 10;
plt.legendLoc = 'northeastoutside';
plt.STI_method = 'ansi'; % sti method (i.e. 'ansi' or 'payton' or 'drullman')
plt.improv = true;
%--------------------------------------------------------------------------
MM = length(methods);RR = length(plt.SNR_arr);WW = length(synWinType);
NN = length(noise_type);SS = length(shift_len);XX = length(xi_type);
%--------------------------------------------------------------------------
latex.horiz =  plt.xTicks;
latex.vert = methods(2:end);
latex.hline =  [0,1,NaN];
latex.vline =  [1];
latex.format = '%.3f';
%--------------------------------------------------------------------------
%%                      Start to process speech files
%--------------------------------------------------------------------------
for ww = 1:WW %length(synWinType)
    mmse_inputs.synWinType = synWinType{ww};
    plt.synWinType = synWinType{ww};
    %---------------------------------------------------------------
    for ss = 1: SS % length(shift_len)
        
        mmse_inputs.frame_len = frame_len(ss);
        mmse_inputs.shift_len = shift_len(ss);
        plt.frame_len = frame_len(ss);
        plt.shift_len = shift_len(ss);
        Tw = Tw_arr(ss);
        %---------------------------------------------------------------
        for n = 1:XX % length(xi_type)
            
            mmse_inputs.xi_type = xi_type{n};
            
            plt.xi_type = xi_type{n};
            
            plt.save_note = [mmse_inputs.xi_type,'_',mmse_inputs.synWinType,...
                '_Syn_',num2str(mmse_inputs.shift_len),'ms_Shif']; % differentiate the tests
            
            output(ss).xi_type =  plt.xi_type;
            %---------------------------------------------------------------
            for nn = 1:NN % length(noise_type)
                %% Read noise file
                noise_no = noise_type(nn);
                
                [noise_init, ~] = audioread([noise_files.folder,...
                    filesep, noise_files.name]); % noise waveform.
                
                plt.noise_name = get_noise_name(noise_no);
                output(ss).noise(nn).noise_name = plt.noise_name;
                %---------------------------------------------------------
                
                %% Read clean speech signal                
                [clean, Fs] = audioread([clean_files(fidx).folder, filesep , clean_files(fidx).name]); % clean waveform.
                output(ss).noise(nn).audio.clean = clean;
                
                %% Create noisy signal according to global SNR
                for rr = 1:RR
                    [noisy,noise]= get_noisy(clean,noise_init,plt.SNR_arr(rr),Tw,'global');
                    
                    mmse_inputs.noisy= noisy;
                    mmse_inputs.noise = noise;
                    mmse_inputs.clean = clean;
                    mmse_inputs.Fs = Fs;
                    
                    S = length(noisy);
                    
                    output(ss).noise(nn).audio.noisy(rr,:) = noisy;
                    
                    clean = clean./norm(clean);
                    noisy = noisy./norm(noise);
                    
                    if output_audios
                        est_file = ['audios' filesep clean_files(1).name(1:end-4),'_',...
                            'clean.wav'];
                        
                        est_file = regexprep(est_file, ' ', '_');
                        audiowrite([exp_dir est_file],clean./max(abs(clean)),Fs,'BitsPerSample',16);
                        
                        est_file = ['audios' filesep clean_files(1).name(1:end-4),'_',...
                            'noisy','_',plt.noise_name,'_',...
                            num2str(plt.SNR_arr(rr)), 'dB.wav'];
                        
                        est_file = regexprep(est_file, ' ', '_');
                        audiowrite([exp_dir est_file],noisy./max(abs(noisy)),Fs,'BitsPerSample',16);
                    end
                end
                
                %% ========================== Apply MMSE estimator ==========================
                for m = 2:MM
                    %---------------------------------------------------------------
                    output(ss).noise(nn).mmse(m-1).name = MMSE_types{para_idx(m,1)};
                    
                    mmse_inputs.speech_prior = speech_priori{para_idx(m,2)};
                    
                    output(ss).noise(nn).mmse(m-1).prior = mmse_inputs.speech_prior;
                    
                    plt.prior{m} = mmse_inputs.speech_prior;
                    
                    [mod_noisy,mod_clean,mod_noise] = myfun{para_idx(m,3)}(mmse_inputs);
                    
                    output(ss).noise(nn).audio.mod_sig(rr,m-1,:) = mod_noisy;
                    
                    S = length(mod_noisy);
                    
                    if output_audios
                        est_file = ['audios' filesep MMSE_types{para_idx(m,1)},'_',...
                            mmse_inputs.speech_prior,'_', plt.xi_type,'_',plt.noise_name,'_',...
                            num2str(plt.SNR_arr(rr)), 'dB.wav'];
                        
                        est_file = regexprep(est_file, ' ', '_');
                        audiowrite([exp_dir est_file],mod_noisy./max(abs(mod_noisy)),Fs,'BitsPerSample',16);
                    end
                end
            end
        end
        
    end
end
%--------------------------------------------------------------------------
%%                              Evaluation
%--------------------------------------------------------------------------
S = length(output.noise.audio.mod_sig(1,1,:));
param.weightmode=1; param.trend=1; % R_d and with noisy amplitude-weight
param.fs = Fs;pesq_mode ='wb';
for p = m_idx
    switch measure{p}
        case 'pesq'
            fprintf('\n Calculating PESQ ');
            PESQScores.Clean = pesq2(clean,clean,Fs);
            PESQScores.Noisy = pesq2(clean,noisy,Fs);
            PESQScores.EM_LSA = pesq2(clean(1:S),squeeze(output.noise.audio.mod_sig(1,1,:)),Fs);
            PESQScores.SG_LSA = pesq2(clean(1:S),squeeze(output.noise.audio.mod_sig(1,2,:)),Fs);
            PESQScores.beta_PC = pesq2(clean(1:S),squeeze(output.noise.audio.mod_sig(1,3,:)),Fs);
            PESQScores.L_STSA = pesq2(clean(1:S),squeeze(output.noise.audio.mod_sig(1,4,:)),Fs);
            PESQScores.N_LSA = pesq2(clean(1:S),squeeze(output.noise.audio.mod_sig(1,5,:)),Fs);
            PESQScores.L_LSA = pesq2(clean(1:S),squeeze(output.noise.audio.mod_sig(1,6,:)),Fs);
            PESQScores.N_PoE = pesq2(clean(1:S),squeeze(output.noise.audio.mod_sig(1,7,:)),Fs);
            PESQScores.L_PoE = pesq2(clean(1:S),squeeze(output.noise.audio.mod_sig(1,8,:)),Fs);
            PESQScores.N_PoE_O = pesq2(clean(1:S),squeeze(output.noise.audio.mod_sig(1,9,:)),Fs);
            PESQScores.L_PoE_O = pesq2(clean(1:S),squeeze(output.noise.audio.mod_sig(1,10,:)),Fs);
            
            PESQScores
            
        case 'seg'
            fprintf('\n Calculating Segmental SNR ');
            type = 'seg';
            SegScores.Clean = quality_score(clean,clean,S,noise,Tw,type);
            SegScores.Noisy = quality_score(clean,noisy,S,noise,Tw,type);
            SegScores.EM_LSA = quality_score(clean,squeeze(output.noise.audio.mod_sig(1,1,:)),S,noise,Tw,type);
            SegScores.SG_LSA = quality_score(clean,squeeze(output.noise.audio.mod_sig(1,2,:)),S,noise,Tw,type);
            SegScores.beta_PC = quality_score(clean,squeeze(output.noise.audio.mod_sig(1,3,:)),S,noise,Tw,type);
            SegScores.L_STSA = quality_score(clean,squeeze(output.noise.audio.mod_sig(1,4,:)),S,noise,Tw,type);
            SegScores.N_LSA = quality_score(clean,squeeze(output.noise.audio.mod_sig(1,5,:)),S,noise,Tw,type);
            SegScores.L_LSA = quality_score(clean,squeeze(output.noise.audio.mod_sig(1,6,:)),S,noise,Tw,type);
            SegScores.N_PoE = quality_score(clean,squeeze(output.noise.audio.mod_sig(1,7,:)),S,noise,Tw,type);
            SegScores.L_PoE = quality_score(clean,squeeze(output.noise.audio.mod_sig(1,8,:)),S,noise,Tw,type);
            SegScores.N_PoE_O = quality_score(clean,squeeze(output.noise.audio.mod_sig(1,9,:)),S,noise,Tw,type);
            SegScores.L_PoE_O = quality_score(clean,squeeze(output.noise.audio.mod_sig(1,10,:)),S,noise,Tw,type);
            
            SegScores
%         case 'UnRMSE'
%             fprintf('\n Calculating UnRMSE score');
%             [UnRMSEScores.Clean , ~] = Calc_UnwrappedPWSNR(clean,clean,param);
%             [UnRMSEScores.Noisy , ~] = Calc_UnwrappedPWSNR(clean,noisy,param);
% 
%             UnRMSEScores.EM_LSA = Calc_UnwrappedPWSNR(clean(1:S),squeeze(output.noise.audio.mod_sig(1,1,:)),param);
%             UnRMSEScores.SG_LSA = Calc_UnwrappedPWSNR(clean(1:S),squeeze(output.noise.audio.mod_sig(1,2,:)),param);
%             UnRMSEScores.beta_PC = Calc_UnwrappedPWSNR(clean(1:S),squeeze(output.noise.audio.mod_sig(1,3,:)),param);
%             UnRMSEScores.L_STSA = Calc_UnwrappedPWSNR(clean(1:S),squeeze(output.noise.audio.mod_sig(1,4,:)),param);
%             UnRMSEScores.N_LSA = Calc_UnwrappedPWSNR(clean(1:S),squeeze(output.noise.audio.mod_sig(1,5,:)),param);
%             UnRMSEScores.L_LSA = Calc_UnwrappedPWSNR(clean(1:S),squeeze(output.noise.audio.mod_sig(1,6,:)),param);
%             UnRMSEScores.N_PoE = Calc_UnwrappedPWSNR(clean(1:S),squeeze(output.noise.audio.mod_sig(1,7,:)),param);
%             UnRMSEScores.L_PoE = Calc_UnwrappedPWSNR(clean(1:S),squeeze(output.noise.audio.mod_sig(1,8,:)),param);
%             UnRMSEScores.N_PoE_O = Calc_UnwrappedPWSNR(clean(1:S),squeeze(output.noise.audio.mod_sig(1,9,:)),param);
%             UnRMSEScores.L_PoE_O = Calc_UnwrappedPWSNR(clean(1:S),squeeze(output.noise.audio.mod_sig(1,10,:)),param);
%             
%             UnRMSEScores
%         case 'UnHPSNR'
%             fprintf('\n Calculating UnHPSNR score');
%             [~, UnHPSNRScores.Clean] = Calc_UnwrappedPWSNR(clean,clean,param);
%             [~, UnHPSNRScores.Noisy] = Calc_UnwrappedPWSNR(clean,noisy,param);
%             
%             [~,UnHPSNRScores.EM_LSA] = Calc_UnwrappedPWSNR(clean(1:S),squeeze(output.noise.audio.mod_sig(1,1,:)),param);
%             [~,UnHPSNRScores.SG_LSA] = Calc_UnwrappedPWSNR(clean(1:S),squeeze(output.noise.audio.mod_sig(1,2,:)),param);
%             [~,UnHPSNRScores.beta_PC] = Calc_UnwrappedPWSNR(clean(1:S),squeeze(output.noise.audio.mod_sig(1,3,:)),param);
%             [~,UnHPSNRScores.L_STSA] = Calc_UnwrappedPWSNR(clean(1:S),squeeze(output.noise.audio.mod_sig(1,4,:)),param);
%             [~,UnHPSNRScores.N_LSA] = Calc_UnwrappedPWSNR(clean(1:S),squeeze(output.noise.audio.mod_sig(1,5,:)),param);
%             [~,UnHPSNRScores.L_LSA] = Calc_UnwrappedPWSNR(clean(1:S),squeeze(output.noise.audio.mod_sig(1,6,:)),param);
%             [~,UnHPSNRScores.N_PoE] = Calc_UnwrappedPWSNR(clean(1:S),squeeze(output.noise.audio.mod_sig(1,7,:)),param);
%             [~,UnHPSNRScores.L_PoE] = Calc_UnwrappedPWSNR(clean(1:S),squeeze(output.noise.audio.mod_sig(1,8,:)),param);
%             [~,UnHPSNRScores.N_PoE_O] = Calc_UnwrappedPWSNR(clean(1:S),squeeze(output.noise.audio.mod_sig(1,9,:)),param);
%             [~,UnHPSNRScores.L_PoE_O] = Calc_UnwrappedPWSNR(clean(1:S),squeeze(output.noise.audio.mod_sig(1,10,:)),param);
%             
%             UnHPSNRScores
        case 'DP'
            fprintf('\n Calculating PD ');
            PDScores.Clean = PD_eval(clean, clean, noisy, param);
            PDScores.Noisy = PD_eval(clean, noisy, noisy, param);
            PDScores.EM_LSA = PD_eval(clean(1:S)', squeeze(output.noise.audio.mod_sig(1,1,:))', noisy(1:S)', param);
            PDScores.SG_LSA = PD_eval(clean(1:S)', squeeze(output.noise.audio.mod_sig(1,2,:))', noisy(1:S)', param);
            PDScores.beta_PC = PD_eval(clean(1:S)', squeeze(output.noise.audio.mod_sig(1,3,:))', noisy(1:S)', param);
            PDScores.L_STSA = PD_eval(clean(1:S)',squeeze(output.noise.audio.mod_sig(1,4,:))', noisy(1:S)', param);
            PDScores.N_LSA = PD_eval(clean(1:S)',squeeze(output.noise.audio.mod_sig(1,5,:))', noisy(1:S)', param);
            PDScores.L_LSA = PD_eval(clean(1:S)',squeeze(output.noise.audio.mod_sig(1,6,:))', noisy(1:S)', param);
            PDScores.N_PoE = PD_eval(clean(1:S)',squeeze(output.noise.audio.mod_sig(1,7,:))', noisy(1:S)', param);
            PDScores.L_PoE = PD_eval(clean(1:S)',squeeze(output.noise.audio.mod_sig(1,8,:))', noisy(1:S)', param);
            PDScores.N_PoE_O = PD_eval(clean(1:S)',squeeze(output.noise.audio.mod_sig(1,9,:))', noisy(1:S)', param);
            PDScores.L_PoE_O = PD_eval(clean(1:S)',squeeze(output.noise.audio.mod_sig(1,10,:))', noisy(1:S)', param);
            
            PDScores
        case 'stoi'
            fprintf('\n Calculating STOI \n');
            % STOI Evaluation
            STOIScores.Clean = stoi(clean, clean, Fs);
            STOIScores.Noisy = stoi(clean, noisy, Fs);
            STOIScores.EM_LSA = stoi(clean(1:S)',squeeze(output.noise.audio.mod_sig(1,1,:))',Fs);
            STOIScores.SG_LSA = stoi(clean(1:S)',squeeze(output.noise.audio.mod_sig(1,2,:))',Fs);
            STOIScores.beta_PC = stoi(clean(1:S)',squeeze(output.noise.audio.mod_sig(1,3,:))',Fs);
            STOIScores.L_STSA = stoi(clean(1:S)',squeeze(output.noise.audio.mod_sig(1,4,:))',Fs);
            STOIScores.N_LSA = stoi(clean(1:S)',squeeze(output.noise.audio.mod_sig(1,5,:))',Fs);
            STOIScores.L_LSA = stoi(clean(1:S)',squeeze(output.noise.audio.mod_sig(1,6,:))',Fs);
            STOIScores.N_PoE = stoi(clean(1:S)',squeeze(output.noise.audio.mod_sig(1,7,:))',Fs);
            STOIScores.L_PoE = stoi(clean(1:S)',squeeze(output.noise.audio.mod_sig(1,8,:))',Fs);
            STOIScores.N_PoE_O = stoi(clean(1:S)',squeeze(output.noise.audio.mod_sig(1,9,:))',Fs);
            STOIScores.L_PoE_O = stoi(clean(1:S)',squeeze(output.noise.audio.mod_sig(1,10,:))',Fs);
            
            STOIScores
%         case 'CSII'
%             fprintf('\n Calculating CSII ');
%             % CSII
%             [iq.Clean,sdr.Clean] = iq3(clean,clean,'ic',Fs);
%             [iq.Noisy,sdr.Noisy] = iq3(noisy,clean,'ic',Fs);
%             [iq.CSA_N,sdr.CSA_N] = iq3(squeeze(output.noise.audio.mod_sig(1,4,:))',clean(1:S)','ic',Fs);
%             [iq.CSA_L,sdr.CSA_L] = iq3(squeeze(output.noise.audio.mod_sig(1,5,:))',clean(1:S)','ic',Fs);
%             [iq.CSA_G,sdr.CSA_G] = iq3(squeeze(output.noise.audio.mod_sig(1,6,:))',clean(1:S)','ic',Fs);
%             [iq.FSA_ADP,sdr.FSA_ADP] = iq3(squeeze(output.noise.audio.mod_sig(1,1,:))',clean(1:S)','ic',Fs);
%             [iq.FSA_AUP,sdr.FSA_AUP] = iq3(squeeze(output.noise.audio.mod_sig(1,2,:))',clean(1:S)','ic',Fs);
%             [iq.FSA_AUPM,sdr.FSA_AUPM] = iq3(squeeze(output.noise.audio.mod_sig(1,3,:))',clean(1:S)','ic',Fs);
%             iq
    end
end

%--------------------------------------------------------------------------
%                     Plot spectrogram of outputs
%--------------------------------------------------------------------------
if plot_spect == 1
    a = {'(a) ','(b) ','(c) ','(d) ','(e) ','(f) ','(g) ','(h) ','(i) ','(j) ','(k) ','(l) ','(m) ','(n) ','(o) ','(p) '};
    
    sl = length(output.noise.mmse)+2;
    si = ceil(sl/2);
    
    f = figure; 
    p = tiledlayout(si,2);
    
    p.XLabel.FontName = 'Times Roman';
    p.XLabel.FontSize = 14;
    p.XLabel.FontWeight = 'bold'; 
    p.YLabel.FontName = 'Times Roman';
    p.TileSpacing = 'compact';
    p.Padding = 'compact';
     
    nexttile;
    myspectrogram(output.noise.audio(1).clean,Fs);
    
    title('(a) Clean speech','Interpreter','latex','fontsize',12,'FontWeight','bold');
    set(gca,'XTick',[0 0.5 1 1.5 2 2.5], 'YTick',[0 4000 8000],'ylim',[0 8000]);
    ylabel('Freq.[kHz]','fontsize',12,'FontWeight','bold');
    
    set(gca,'xticklabel',{[]},'yticklabel',{'0' '4' '8'}) ;
    
    nexttile;
  
    myspectrogram(output.noise.audio(1).noisy(1,:),Fs);
    
    title(['(b) Noisy speech, ' plt.noise_name ' SNR = ' num2str(plt.SNR_arr(1)) ' dB'],'Interpreter','latex','fontsize',12,'FontWeight','bold');
    set(gca,'XTick',[0 0.5 1 1.5 2 2.5], 'YTick',[0 4000 8000],'ylim',[0 8000]);
    ylabel('Freq.[kHz]','fontsize',12,'FontWeight','bold');
    set(gca,'xticklabel',{[]},'yticklabel',{'0' '4' '8'});
    
    k = 3;
    for r = 2:si
        for c = 1:2
            nexttile;
            myspectrogram(squeeze(output.noise.audio(1).mod_sig(1,k-2,:)),Fs);
            
            title([a{k} methods{k-1}],'Interpreter','latex','fontsize',12,'FontWeight','bold');
            set(gca,'XTick',[0 0.5 1 1.5 2 2.5], 'YTick',[0 4000 8000],'ylim',[0 8000]);
            
            ylabel('Freq.[kHz]','fontsize',12,'FontWeight','bold');
            set(gca,'xticklabel',{[]},'yticklabel',{'0' '4' '8'}) ;
            
            if k == sl-1|| k== sl
                set(gca,'xticklabel',{[0 0.5 1 1.5 2 2.5]}) ;
                xlabel('Time[s]','fontsize',12,'FontWeight','bold');
%                 c.TickDirection = 'out';

                preset_name='red-white-blue'; 
                J = customcolormap_preset(preset_name);
               
                colormap(J);
                c = colorbar; 
                c.Location = 'southoutside';
            end
            
            k=k+1;
            
            if k == sl+1
                break
            end
        end
    end
    
    save_name = strcat('spect_MMSE_compare_', plt.noise_name,'_',num2str(plt.SNR_arr(1)),'_dB','.eps');
    
    print('-depsc', '-r600',sprintf(strcat(plt.save_dir,save_name)));
    
end
%% ============================== Functions ===============================
function noise_t = get_noise_name(noise_no)
switch noise_no
    case 1
        noise_t = 'Sinusoid';
    case 2
        noise_t = 'Pink noise';
    case 3
        noise_t = 'White noise';
    case 4
        noise_t = 'White noise_n6dB';
    case 5
        noise_t = 'White noise_n12dB';
    case 6
        noise_t = 'Speech noise';
    case 7
        noise_t = 'M 109';
    case 8
        noise_t = 'Buccaneer 190 Knots 1000 Feet';
    case 9
        noise_t = 'Leopard 2';
    case 10
        noise_t = 'Wheel carrier';
    case 11
        noise_t = 'Buccaneer 450 Knots 300 Feet';
    case 12
        noise_t = 'Lynx';
    case 13
        noise_t = 'Leopard 1';
    case 14
        noise_t = 'Operation room of destroyer';
    case 15
        noise_t = 'Engine room of destroyer';
    case 16
        noise_t = 'Machine gun repeated';
    case 17
        noise_t = 'HF radio channel';
    case 18
        noise_t =  'STITEL test signal';
    case 19
        noise_t = 'Voice babble';
    case 20
        noise_t = 'F-16 two-seat';
    case 21
        noise_t = 'Car Factory electrical welding';
    case 22
        noise_t = 'Car Factory production hall';
    case 23
        noise_t = 'Car Volvo-340 asphalt road';
    case 24
        noise_t = 'Car Volvo-340 brick road';
    otherwise
        noise_t = 'unknown';
end
end
function score = quality_score(clean, mod, S, noise, Tw,type)

global min_db max_db
sd_fun = @(ref,est)min(max(10*log10(ref), min_db), max_db)-min(max(10*log10(est), min_db), max_db);
% mod = mod./norm(mod);
% clean = clean./norm(clean);
% noise = noise./norm(noise);
%--------------------------------------------------------------------------
switch type
    case{'seg','seg_S','seg_N','Seg','SegSNR'}
        %% no-overlapping
        [frame_clean] = framing(clean,Tw/2,Tw/2,'window_type','rect');
        [frame_noise] = framing(noise,Tw/2,Tw/2,'window_type','rect');
        [frame_mod] = framing(mod,Tw/2,Tw/2,'window_type','rect');
        %--------------------------------------------------------------------------
        % function snr = get_SNR(sig,noise,mod,varargin)
        score = get_SNR(frame_clean,frame_noise,frame_mod,'type',type);
    case{'SD','sd'}
        % no-overlapping
        [frame_clean] = framing(clean,Tw,Tw,'window_type','rect');
        [frame_mod] = framing(mod,Tw,Tw,'window_type','rect');
        
        score = mean(mean(bsxfun(sd_fun,frame_clean,frame_mod),2));
    otherwise
        score = 0;
end
end
%--------------------------------------------------------------------------
%==========================================================================
%--------------------------------E.O.F.------------------------------------
%==========================================================================
