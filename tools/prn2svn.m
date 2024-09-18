function svns = prn2svn(epoch,prns)
%PRN2SVN Given #s of psuedorandom noise codes for GPS satellites and the epoch 
%of use, returns the corresponding space vehicle # associated at that time.
%   Inputs:
%    - epoch; seconds past J2000
%    - prns; list of PRN numbers
arguments
    epoch   (1,1)   double {mustBePositive}
    prns    (1,:)   {mustBePositive,mustBeInteger}
end

% convert epoch to seconds since the Unix epoch, 1/1/1970 00:00:00
epoch = datetime(epoch, 'ConvertFrom', 'epochtime', 'Epoch', '2000-01-01 12:00:00');

% Open data file and throw out the header
fid = fopen('GPS_PRN_Usage_Data.txt','r');
for ii = 1:8
    % throw out the headers
    hdr = fgetl(fid);
    if ii < 4
        % Print the info?
        % fprintf('%s\n',hdr);
    end
end

% For now we discard the clock information
% SVN #/T PRN Start_Time - Stop_Time   => Clk_Turn-on
[A, COUNT] = fscanf(fid,'%d %*s %d %d %*s %d %*s %d');
Data = reshape(A,5,[])';

svns = zeros(size(prns));

for i = 1:length(prns)

    % Find data by prn
    prn = prns(i);
    prn_rows = find(Data(:,2) == prn);
    unixTimeStart = Data(prn_rows,3);
    unixTimeEnd = Data(prn_rows,4);
    
    % Time of applicability
    % WARNING: this code (maybe?) does not handle leap seconds properly
    % times are in seconds since Unix Epoch (1/1/1970 00:00:00)
    % Get the datetime representing the UNIX epoch
    unixepoch = datetime(1970,1,1,0,0,0);
    
    % Add seconds to start/end to the datetime representing the epoch
    tstart = unixepoch + seconds(unixTimeStart);
    tend = unixepoch + seconds(unixTimeEnd);
    
    doa_row = prn_rows(tstart < epoch & epoch < tend);
   
    if ~isempty(doa_row)
        % this prn was transmitting at this epoch, set svn to not zero
        svns(i) = Data(doa_row,1);
    end
end

