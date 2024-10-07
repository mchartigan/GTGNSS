classdef (Abstract) Propagator < handle
    %PROPAGATOR Abstract class for defining various types of propagators.
    %This applies most to orbital propagators and clock propagators.
    
    properties (Abstract)
        t0                          % starting time of propagation
        x0                          % starting state of propagation
        ts                          % times of resultant propagation
        xs                          % states of resultant propagation
    end
    
    methods (Abstract)
        out = run(obj,tf,n)         % run to final time w/ equispaced steps
        out = runat(obj,ts)         % run to specific times
        out = modelfit(obj)         % fit an approximation model to the result
    end
end

