function [noisy,noise_new,segSNR_new]= get_noisy(obj)
            % Be careful this noise is depend on the original noise signal
            % type, as well as the type of SNR. The original noise is
            % modifed according to the requirement. see addnoise function.
            
            noisePara.noise_init = obj.noise_init;
            framePara.window_type = 'rect';
            
            switch obj.SNR_type
                case 'seg'
                    segSNR = obj.SNR;    
                    SNR_Para.type = 'seg';
                    % get two initial values, to start the converge process
                    SNR_ratio = 2; % approximated ratio between the seg and global SNR
                    globSNR_old = segSNR*SNR_ratio; % first initial value
                    globSNR_new = globSNR_old + 1; % second initial value
                    
                    % create noisy signal using global SNR
                    noisy = obj.addnoise(obj.orig_Sig,globSNR_old,noisePara);
                    % no-overlapping
                    [frame_clean] = framing(obj.orig_Sig,obj.Tw,obj.Tw,framePara);
                    [frame_noisy] = framing(noisy,obj.Tw,obj.Tw,framePara);
                    
                    % find seg SNR
                    SNR_seg_old = obj.getSNR(frame_clean,frame_noisy-frame_clean,SNR_Para);
                    err_old = segSNR-SNR_seg_old;
%% check convergence                   
%                     globSNR_arr = [];
%                     globSNR_arr = [globSNR_arr globSNR_old globSNR_new];
%                     segSNR_arr = [];
%                     segSNR_arr = [segSNR_arr SNR_seg_old];
%                     err_arr = [];      
%                     err_arr = [err_arr err_old];
%% find the SNR that when the target segSNR within the threshold                   
                    while (1)
                        [noisy,noise_new] = obj.addnoise(obj.orig_Sig,globSNR_new,noisePara);
                        [frame_noisy] = framing(noisy,obj.Tw,obj.Tw,framePara);
                        
                        segSNR_new = obj.getSNR(frame_clean,frame_noisy-frame_clean,SNR_Para);
                        err_new = segSNR-segSNR_new;

                        if abs(err_new) < obj.SNR_thresh  
%                            fprintf('segSNR : %6.4f . done . \n',segSNR_new); 
                           break; 
                        end
                        
                        err_diff = err_new-err_old;
                        SNR_temp = globSNR_new - err_new/(err_diff);
                        err_old = err_new;
                        globSNR_new = SNR_temp;
                                                               
                    end               
                case {'global','overall'}
                    segSNR_new = [];
                    [noisy,noise_new] = obj.addnoise(obj.orig_Sig,obj.SNR,noisePara);
%                     fprintf('global SNR : %6.4f . done . \n',obj.SNR);
            end
end
        
