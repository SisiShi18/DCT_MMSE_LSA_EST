function [frame_data,idx_mat]= framing (signal,Tw,Ts,varargin)
    in = inputParser;
    % Only use modified version when window shift divides the
    % window length evenly
    default_opt = 'hamming';
    valid_opt = {'hamming','hanning','rect','mod_hamming','mod_hanning','mod_rect'};
    check_opt = @(x) any(validatestring(x,valid_opt));
    addParameter(in,'window_type',default_opt,check_opt);

    in.parse(varargin{:});
    window_type = in.Results.window_type;

    % framing data set and windowing
    N = length (signal); %apply on raw data, for computing the raw energy of each frame
    end_idx_row = N-(Tw);
    idx_mat = bsxfun(@plus,(0:Ts:end_idx_row)',1:Tw);
    frame_data = signal(idx_mat);

    switch window_type
        % Only use modified version when window shift divides the
        % window length evenly
        case 'hamming'
            frame_data = bsxfun(@times,frame_data,hamming(Tw)'); % apply window
        case 'hanning'
            frame_data = bsxfun(@times,frame_data,hanning(Tw)');
        case 'rect'
            %frame_data = bsxfun(@times,frame_data,rectwin(obj.Tw)');
            frame_data = bsxfun(@times,frame_data,ones(1,Tw)); % apply window
        case 'mod_hamming'
            frame_data = bsxfun(@times,frame_data,obj.modified_window(Tw,Ts,'type','hamming')'); % apply window
        case 'mod_hanning'
            frame_data = bsxfun(@times,frame_data,obj.modified_window(Tw,Ts,'type','hanning')');
        case 'mod_rect'
            frame_data = bsxfun(@times,frame_data,obj.modified_window(Tw,Ts,'type','rect')');
    end
end