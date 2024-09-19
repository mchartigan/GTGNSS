function gcfs = svn2gcf(svns)
%SVN2GCF Given GPS Space Vehicle #s, return the associated antenna Gain
%Correction Factors (GCFs) for the L1 and L2 antennae.
%   Input:
%    - svns, list of space vehicle #s
arguments
    svns    (1,:)   {mustBePositive,mustBeInteger}
end

file = importdata("GCF_from_SVN.txt");
data = file.data;

means = mean(data,1);

gcfs = diag(means(2:3)) * ones(2, length(svns));
for i=1:length(svns)
    idx = find(data(:,1) == svns(i));
    if size(idx,1) == 0
        warning("svn2gcf:invalidSVN", ...
            "Provided SVN %d is not in the database.", svns(i));
    else
        gcfs(:,i) = data(idx,2:3)';
    end
end
end

