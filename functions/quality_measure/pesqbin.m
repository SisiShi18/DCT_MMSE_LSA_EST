function [ scores ] = pesqbin( reference, degraded, fs, mode )
% PESQBIN MATLAB wrapper for the PESQ binary version 2.0.
%
%   [SCORES]=PESQBIN(REF,DEG,FS,MODE) calls the PESQ binary for 
%   reference signal REF and degraded signal DEG, each sampled 
%   at FS Hz and specified as either array of samples or audio 
%   filename. Note that this wrapper routine expects that PESQ 
%   version 2.0 [1] is installed on the system and added to
%   system path. PESQ 2.0 supports both NB-PESQ (narrowband 
%   PESQ measure) as well as WB-PESQ (wideband PESQ measure).
%   This wrapper supports both modes through the MODE parameter.
%   For narrowband mode this function returns two element
%   array: [ PESQ_MOS, MOS_LQO ], while for the wideband mode
%   single MOS_LQO score is returned. Works on both Linux and
%   Windows operating systems, again, as long as PESQ binary
%   is installed, see [1].
%
%   Inputs
%           REF is the reference signal as vector, 
%               or as filename of an audio file (wav or raw)
%
%           DEG is the degraded signal as vector, 
%               or as filename of an audio file (wav or raw)
%
%           FS is the sampling frequency in Hz for REF and DEG
%              (optional if REF and/or DEG are passed as filenames)
%
%           MODE specifies narrowband: 'nb' or wideband: 'wb'
%                mode of operation
%
%   Outputs 
%           SCORES is a two element array: [ PESQ_MOS, MOS_LQO ] 
%                  for narrowband mode, or a scalar for the 
%                  wideband mode: MOS_LQO 
%
%   Example
%           % read reference and degraded samples from audio
%           [ reference, fs ] = wavread( 'speech.wav' );
%           [ degraded, fs ] = wavread( 'speech_bab_0dB.wav' );
% 
%           % compute NB-PESQ and WB-PESQ scores 
%           nb = pesqbin( reference, degraded, fs, 'nb' );
%           wb = pesqbin( reference, degraded, fs, 'wb' );
%
%           % display results to screen
%           fprintf( 'NB PESQ MOS = %5.3f\n', nb(1) );
%           fprintf( 'NB MOS LQO  = %5.3f\n', nb(2) );
%           fprintf( 'WB MOS LQO  = %5.3f\n', wb );
%
%           % example output
%           %    NB PESQ MOS = 1.969
%           %    NB MOS LQO  = 1.607
%           %    WB MOS LQO  = 1.083
%
%   PESQ MOS mappings to LQO:
%
%       % P.862.1->P.800.1 (PESQ_MOS->MOS_LQO)
%           MOS_LQO = 0.999 + ( 4.999-0.999 ) ./ ( 1+exp(-1.4945*PESQ_MOS+4.6607) ) 
%
%       % P.862.2->P.800.1 (PESQ_MOS->MOS_LQO)
%           MOS_LQO = 0.999 + ( 4.999-0.999 ) ./ ( 1+exp(-1.3669*PESQ_MOS+3.8224) ) 
%
%   [1] ITU-T (2005), "P.862: Revised Annex A - Reference implementations 
%       and conformance testing for ITU-T Recs P.862, P.862.1 and P.862.2"
%       url: http://www.itu.int/rec/T-REC-P.862-200511-I!Amd2/en
%
%   See also TEST_PESQBIN.

%   Author: Kamil Wojcicki, UTD, November 2011


    % usage information
    usage = 'usage: [ pesq_mos ] = pesqbin( reference, degraded, fs, mode );';

    % default settings and 'lite' input validation
    switch( nargin )
    case { 0, 1 }, error( usage );
    case 2, mode='nb'; if ~isstr(reference) || ~isstr(degraded), error( usage ); end;
    case 3, mode='nb';
    case 4, 
    otherwise, error( usage );
    end 

    % Select platform dependent execution path
    if isunix()
        % use Linux-tailored routines
        scores = pesqbin_linux( reference, degraded, fs, mode, 'pesq' );
    else 
        % use Window-tailored routines
        scores = pesqbin_windows( reference, degraded, fs, mode, 'pesq.exe' );
    end


% MATLAB wrapper for PESQ binary on the Linux operating system
function [ scores ] = pesqbin_linux( reference, degraded, fs, mode, binary )


    % temporary directory in which wav files will be stored (if needed)
    % and in which the PESQ binary will be executed, hence also where
    % the pesq_results.txt logfile will also be written
    tmpdir = sprintf( '%s%s', filesep, 'tmp' );
    tmpdir = regexprep( tmpdir, '\n.*', '' );


    % default precision for audio samples (in number of bits)
    nbits = 16;

    
    % determine if the user provided an audio signal as vector 
    % as the reference signal, or an audio filename instead
    if isstr(reference)
        file.reference = reference; 
        [ ~, fs] = audioread( file.reference );
    else
        file.reference = sprintf( '%s%s%s', tmpdir, filesep, '~80b4eb734d.wav' );
        audiowrite( file.reference, 0.999*reference./max(abs(reference)), fs);
    end % isstr(reference)


    % determine if the user provided an audio signal as vector 
    % as the degraded signal, or an audio filename instead
    if isstr(degraded)
        file.degraded = degraded; 
        [ ~, fs] = audioread( file.degraded );
    else
        file.degraded = sprintf( '%s%s%s', tmpdir, filesep, '~be4dfad7fba.wav' );
        audiowrite( file.degraded, 0.999*degraded./max(abs(degraded)), fs );
    end % isstr(degraded)


    % select conditional mode of processing: 
    % narrowband mode or wideband mode
    switch lower( mode )
    
        % computed prediction for narrowband speech
        case { [], '', 'nb', '+nb', 'narrowband', '+narrowband' }

            command = sprintf( 'cd %s && %s +%i %s %s && cd -', tmpdir, binary, fs, file.reference, file.degraded );
    
            [ status, stdout ] = system( command );

            if status~=0 
                warning( 'The %s binary exited with error code %i:\n%s\n', binary, status, stdout );
            end

            scores = stdout2scores( stdout, mode );
    
        % computed prediction for wideband speech
        case { 'wb', '+wb', 'wideband', '+wideband' }

            command = sprintf( 'cd %s && %s +%i +wb %s %s && cd -', tmpdir, binary, fs, file.reference, file.degraded );
    
            [ status, stdout ] = system( command );

            if status~=0 
                warning( 'The %s binary exited with error code %i.\n%s\n', binary, status, stdout );
            end

            scores = stdout2scores( stdout, mode );
    
        % otherwise declare an error
        otherwise
            error( sprintf('Mode: %s is unsupported!\n',mode) );
    
    end % switch lower( mode )


% MATLAB wrapper for PESQ binary on the Windows operating system
function [ scores ] = pesqbin_windows( reference, degraded, fs, mode, binary )


    % temporary directory in which wav files will be stored (if needed)
    % and in which the PESQ binary will be executed, hence also where
    % the pesq_results.txt logfile will also be written
    %tmpdir = sprintf( '%s%s%s%s%s', 'C:', filesep, 'WINDOWS', filesep, 'Temp' );
    [ ~, tmpdir ] = system( 'echo %TEMP%' ); % e.g., C:\DOCUME~1\USERNAME\LOCALS~1\Temp
    tmpdir = regexprep( tmpdir, '\n.*', '' );
    expdir = pwd();


    % default precision for audio samples (in number of bits)
    nbits = 16;

    
    % determine if the user provided an audio signal as vector 
    % as the reference signal, or an audio filename instead
    if isstr(reference)
        file.reference = reference; 
        [ ~, fs] = audioread( file.reference );
    else
        file.reference = sprintf( '%s%s%s', tmpdir, filesep, '~80b4eb734d.wav' );
        audiowrite( file.reference, 0.999*reference./max(abs(reference)), fs );
    end % isstr(reference)


    % determine if the user provided an audio signal as vector 
    % as the degraded signal, or an audio filename instead
    if isstr(degraded)
        file.degraded = degraded; 
        [ ~, fs ] = audioread( file.degraded );
    else
        file.degraded = sprintf( '%s%s%s', tmpdir, filesep, '~be4dfad7fba.wav' );
        audiowrite( file.degraded, 0.999*degraded./max(abs(degraded)), fs );
    end % isstr(degraded)


    % select conditional mode of processing: 
    % narrowband mode or wideband mode
    switch lower( mode )
    
        % computed prediction for narrowband speech
        case { [], '', 'nb', '+nb', 'narrowband', '+narrowband' }

            command = sprintf( 'pushd %%CD%% && cd %s && %s +%i %s %s && popd', tmpdir, binary, fs, file.reference, file.degraded );
    
            cd( tempdir );
            [ status, stdout ] = system( command );
            cd( expdir );

            if status~=0 
                warning( 'The %s binary exited with error code %i:\n%s\n', binary, status, stdout );
            end

            scores = stdout2scores( stdout, mode );
    
        % computed prediction for wideband speech
        case { 'wb', '+wb', 'wideband', '+wideband' }

            command = sprintf( 'pushd %%CD%% && cd %s && %s +%i +wb %s %s && popd', tmpdir, binary, fs, file.reference, file.degraded );
    
            cd( tempdir );
            [ status, stdout ] = system( command );
            cd( expdir );

            if status~=0 
                warning( 'The %s binary exited with error code %i:\n%s\n', binary, status, stdout );
            end

            scores = stdout2scores( stdout, mode );
    
        % otherwise declare an error
        otherwise
            error( sprintf('Mode: %s is unsupported!\n',mode) );
    
    end % switch lower( mode )


% The PESQ binary outputs results, along with some other
% information to STDOUT. This function is used to extract
% the actual scores from the STDOUT output of the PESQ binary.
function [ scores ] = stdout2scores( stdout, mode )


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


% EOF
