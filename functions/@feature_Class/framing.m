function [frame_data,idx_mat]= framing (obj,varargin)
    in = inputParser;
    % Only use modified version when window shift divides the
    % window length evenly
    default_opt = 'hamming';
    valid_opt = {'hamming','hanning','rect','mod_hamming','mod_hanning','mod_rect'};
    check_opt = @(x) any(validatestring(x,valid_opt));
    addParameter(in,'winType',default_opt,check_opt);
    
    addParameter(in,'signal',obj.signal,@(x) isnumeric(x));
    
    in.parse(varargin{:});
    signal = in.Results.signal;
    winType = in.Results.winType;

    % framing data set and windowing
    if obj.preE
        N = length (obj.preEmphasis); %apply pre-emphesis
        end_idx_row = N-(obj.Tw);
        idx_mat = bsxfun(@plus,(0:obj.Ts:end_idx_row)',1:obj.Tw);
        frame_data = obj.preEmphasis(idx_mat);

        frame_data = bsxfun(@times,frame_data,hamming(obj.Tw)'); % apply window
    else
        N = length (signal); %apply on raw data, for computing the raw energy of each frame
        end_idx_row = N-(obj.Tw);
        idx_mat = bsxfun(@plus,(0:obj.Ts:end_idx_row)',1:obj.Tw);
        frame_data = signal(idx_mat);
        
        switch winType
            % Only use modified version when window shift divides the
            % window length evenly
            case 'hamming'
               frame_data = bsxfun(@times,frame_data,hamming(obj.Tw)'); % apply window
            case 'hanning'
               frame_data = bsxfun(@times,frame_data,hanning(obj.Tw)'); 
            case 'rect'
               frame_data = bsxfun(@times,frame_data,ones(1,obj.Tw)); % apply window
            case 'mod_hamming'
            frame_data = bsxfun(@times,frame_data,obj.modified_window(obj.Tw,obj.Ts,'type','hamming')'); % apply window            
            case 'mod_hanning'
            frame_data = bsxfun(@times,frame_data,obj.modified_window(obj.Tw,obj.Ts,'type','hanning')');
            case 'mod_rect'
            frame_data = bsxfun(@times,frame_data,obj.modified_window(obj.Tw,obj.Ts,'type','rect')');
        end
    end
end