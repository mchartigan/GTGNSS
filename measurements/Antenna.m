classdef Antenna < handle
    %ANTENNA Stores properties of communications antenna.
    
    properties
        % dBi, gain above isotropic antenna (default 0)
        gain    (1,1)   double = 0
        % rad, antenna mask angle (pi/2 - off-boresight angle)
        mask    (1,1)   double = 5 * pi/180
    end
end

