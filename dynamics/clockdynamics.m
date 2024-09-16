function dxdt = clockdynamics(x)
%CLOCKDYNAMICS Continuous-time dynamics of oscillator as a random walk- and
%run-type process.
%   Inputs:
%    - t; simulation time
%    - x; clock state [phase offset (s); frequency offset (s/s); frequency
%         drift (1/s)]
arguments
    x (3,1) double
end

dxdt = [0 1 0; 0 0 1; 0 0 0] * x;
end

