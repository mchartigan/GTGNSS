% setup_GTGNSS.m
% Author: Mark Hartigan
% Date  : December 12, 2024
% Description:
%    Configure the GTGNSS toolbox for first-time setup or redeployment.

fprintf("Beginning toolbox initialization...\n");

% get full path to SPKPATH file
here = split(mfilename('fullpath'), 'GTGNSS_setup');
spkpath = strcat(here{1}, 'res/SPKPATH');

% check if SPKPATH exists -- this file stores the path to the directory
% where you store SPICE kernels on your local machine. This is so that
% paths can be referenced locally in scripts with relative ease while
% maintaining modularity across different machines. If it doesn't exist,
% you'll be prompted to enter a path, which will be saved in ./res/SPKPATH.
if ~isfile(spkpath)
    candidate = uigetdir('', "Select SPICE folder " + ...
        "(choose parent folder of 'generic/mk/generic_*.tm' to run" + ...
        " included examples)");

    fid = fopen(spkpath, "w");
    fprintf(fid, "%s", candidate);
    fclose(fid);

    fprintf("SPICE directory set to: %s\n", candidate);
end

fprintf("Toolbox initialization complete.\n");