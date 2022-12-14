function h = imagesc(obj,varargin)
    % Color coded image plot of traces in object.
    % USAGE:
    %   h = obj.plot('attr')
    %   or s.imagesc('mag')
    % where:
    %   h    = plot handle as returned by buildin plot function
    %          (optional)
    %   attr = One of the following strings: mag,phase
    %          default: mag

    if isempty(varargin)
        attr = 'mag';
    else
        attr = varargin{1};
    end

    switch attr
        case 'mag'
            plotHandle = imagesc (rot90(obj.mag));
        case 'phase'
            plotHandle = imagesc (rot90(obj.pha));       
        case 'abs'
            plotHandle = imagesc (rot90(obj.mag));
        case 'polar'
            plotHandle = imagesc (rot90(obj.pha));
        case 'power'
            plotHandle = imagesc (rot90(obj.PSD));
        case 'sse'
            plotHandle = imagesc (rot90(10*log10(obj.SSE))); colorbar ; colormap(1-gray(256));
    end
    if nargout ~= 0;
        h = plotHandle;
    end
end