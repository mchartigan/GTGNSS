# Notes on formatting of clock JSON files

Data is stored as standard JSON files; minimum required data is listed in the example below:

    {
        "frequency" : ...,          % output frequency of oscillator, Hz
        "stability" : {             % Allan deviation characteristics
            "int" : [...],              % observation intervals for frequency stability, s
            "dev" : [...]               % deviation for each interval, in s/s
        },
        "aging" : ...,              % frequency drift, or aging, in s/s/day
        "phase-noise" : {           % oscillator phase noise characteristics
            "freq"  : [...],            % frequency offset from nominal, Hz
            "noise" : [...]             % noise, in dBc/Hz
        }
    }