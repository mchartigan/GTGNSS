classdef (Abstract) Measurement < handle
    %MEASUREMENT Abstract class for defining various types of measurements.
    %Newer versions of navigation filters require this.

    properties (Abstract)
        dim (1,1)   {mustBeInteger,mustBePositive}  % dimension of state
    end

    methods (Abstract)
        y = computemeas(obj,t,x)    % get computed measurements
        H = measpartials(obj,t,x)   % get measurement partial matrix
    end
end