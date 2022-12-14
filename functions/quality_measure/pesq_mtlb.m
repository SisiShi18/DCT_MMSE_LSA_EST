function [ scores ] = pesq_mtlb( reference, degraded, fs, mode )
%   PESQ2_MTLB is MATLAB driver for speech quality calculation
% in accordance with PESQ version 2.0 [1] on Windows operating system.
%   PESQ2_MTLB uses executable file PESQ2.EXE, which should be placed
% at the same directory where PESQ2_MTLB is present.
%   PESQ2.EXE can be compiled from C sources, which are freely available
% for downloading from ITU-T website [1].
%   PESQ2_MTLB uses STDOUT2SCORES function developed by Kamil Wojcicki [2].
%
%   PESQ version 2.0 supports both NB-PESQ (narrowband PESQ measure)
%   as well as WB-PESQ (wideband PESQ measure).
%   This driver supports both modes through the MODE parameter.
%
%   SYNTAX:
%            [ SCORES ]=PESQ2_MTLB( REF, DEG, FS, MODE, BNR, PATH )
%
%   INPUTS:
%       REF is  file name of reference signal (wav or raw)
%           (file extension must be explicitly stated)
%
%       DEG is file name of degraded signal (wav or raw)
%           (file extension must be explicitly stated)
%
%       FS is the sampling frequency in Hz for REF and DEG
%           (FS values may be only 8000 or 16000)
%
%       MODE specifies narrowband: 'nb' or wideband: 'wb'
%
%       BNR is name (string variable) of binary file PESQ2.EXE
%
%       PATH is path to folder which contains reference and degraded audio
%       files. Full path to the folder may be pointed here if it is necessary.
%
%   OUTPUTS:
%       SCORES is a two element array: [ PESQ_MOS, MOS_LQO ]
%       for narrowband mode, or a scalar for the
%       wideband mode: MOS_LQO
%
%       PESQ_RESULTS.TXT is text file which contains results history
%       and is saved and renewed at the same foler, were PESQ2.EXE is placed
%
%   EXAMPLES of computing NB-PESQ and WB-PESQ scores for wav-files
%       and raw-files see at TEST_PESQ2_MTLB m-file
%
% [1] ITU-T (2005), "P.862: Revised Annex A - Reference implementations
%     and conformance testing for ITU-T Recs P.862, P.862.1 and P.862.2"
%     URL: http://www.itu.int/rec/T-REC-P.862-200511-I!Amd2/en
% [2] Wojcicki K. PESQ MATLAB Wrapper.
%     URL: http://www.mathworks.com/matlabcentral/fileexchange/33820-pesq-matlab-wrapper
%
%   See also TEST_PESQ2_MTLB.
%
%   Author: Arkadiy Prodeus, email: aprodeus@gmail.com, Jule 2014
%==========================================================================
% Select platform dependent execution path
if isunix()
    % use Linux-tailored routines
    scores = pesqbin_linux( reference, degraded, fs, mode, 'pesq' );
else
    % use Window-tailored routines
    scores = pesqbin_windows( reference, degraded, fs, mode, 'pesq.exe' );
end
%==========================================================================
% MATLAB wrapper for PESQ binary on the Linux operating system
    function [ scores ] = pesqbin_linux( reference, degraded, fs, mode, binary )
        % identify current folder
        expdir = pwd();
        % select conditional mode of processing:
        % narrowband mode or wideband mode
        switch lower( mode )
            
            % computed prediction for narrowband speech
            case { [], '', 'nb', '+nb', 'narrowband', '+narrowband' }
                
                command = sprintf( 'cd %s && %s +%i %s %s && cd -', expdir, binary, fs, reference, degraded );
                
                [ status, stdout ] = system( command );
                
                if status~=0
                    warning( 'The %s binary exited with error code %i:\n%s\n', binary, status, stdout );
                end
                
                scores = stdout2scores( stdout, mode );
                
                % computed prediction for wideband speech
            case { 'wb', '+wb', 'wideband', '+wideband' }
                
                command = sprintf( 'cd %s && %s +%i +wb %s %s && cd -', expdir, binary, fs, reference, degraded );
                
                [ status, stdout ] = system( command );
                
                if status~=0
                    warning( 'The %s binary exited with error code %i.\n%s\n', binary, status, stdout );
                end
                
                scores = stdout2scores( stdout, mode );
                
                % otherwise declare an error
            otherwise
                error( sprintf('Mode: %s is unsupported!\n',mode) );
        end
    end
%==========================================================================
    function [ scores ] = pesqbin_windows( reference, degraded, fs, mode, binary )
        % identify current folder
        expdir = pwd();
        % select conditional mode of processing:
        % narrowband mode or wideband mode
        switch lower( mode )
            
            % computed prediction for narrowband speech
            case { [], '', 'nb', '+nb', 'narrowband', '+narrowband' }
                
                command = sprintf( 'pushd %%CD%% && cd %s && %s +%i %s %s && popd', expdir, binary, fs, reference, degraded );
                
                [ status, stdout ] = system( command );
                
                if status~=0
                    warning( 'The %s binary exited with error code %i:\n%s\n', binary, status, stdout );
                end
                
                scores = stdout2scores( stdout, mode );
                
                % computed prediction for wideband speech
            case { 'wb', '+wb', 'wideband', '+wideband' }
                
                command = sprintf( 'pushd %%CD%% && cd %s && %s +%i +wb %s %s && popd', expdir, binary, fs, reference, degraded );
                
                [ status, stdout ] = system( command );
                
                if status~=0
                    warning( 'The %s binary exited with error code %i:\n%s\n', binary, status, stdout );
                end
                
                scores = stdout2scores( stdout, mode );
                
                % otherwise declare an error
            otherwise
                error( sprintf('Mode: %s is unsupported!\n',mode) );
                
        end % switch lower( mode )
    end
%==========================================================================
    function [ scores ] = stdout2scores( stdout, mode )
        % The PESQ binary outputs results, along with some other
        % information to STDOUT. This function is used to extract
        % the actual scores from the STDOUT output of the PESQ binary.
        
        % Author: Kamil Wojcicki, UTD, November 2011
        % URL: http://www.mathworks.com/matlabcentral/fileexchange/33820-pesq-matlab-wrapper
        
        % select conditional mode of processing:
        % narrowband mode or wideband mode
        switch lower( mode )
            
            % computed prediction for narrowband speech
            case { [], '', 'nb', '+nb', 'narrowband', '+narrowband' }
                tag = 'P.862 Prediction (Raw MOS, MOS-LQO):  = ';
                defaults = [ NaN, NaN ];
                
                % computed prediction for wideband speech
            case { 'wb', '+wb', 'wideband', '+wideband' }
                tag = 'P.862.2 Prediction (MOS-LQO):  = ';
                defaults = NaN;
                % otherwise declare an error
            otherwise
                error( sprintf('Mode: %s is unsupported!\n',mode) );
                
        end % switch lower( mode )
        
        % length of standard output (in characters)
        S = length( stdout );
        
        % location of MOS score predictions
        idx = strfind( stdout, tag );
        
        % sanity check...
        if isempty(idx) || length(idx)~=1 || idx>S
            scores = defaults;
            return;
        end
        
        % truncate to keep MOS info at the start
        stdout = stdout(idx+length(tag):end);
        
        % scan for at most two floats
        scores = sscanf( stdout, '%f', [1,2] );
        
        % sanity check...
        if isempty( scores )
            scores = defaults;
        end
    end
end
%==========================================================================
% EOF