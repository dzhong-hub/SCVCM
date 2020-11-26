function [mtws,mstd,dtws] = averageFromInputFiles(dayfiles, row, col, data_type, nodata)
%
% function for averaging the daily EALCO TWS to monthly TWS
%
% Input:
%   dayfiles -- a file list that contains all daily EALCO TWS for the
%               period (month)
%   row, col -- demension of the EALCO TWS grid saved in the daily file
%   data_type -- data type saved in the files
%
% Output:
%   mtws -- mean EALCO TWS for the period (month)
%   mstd -- standard deviation of mtws
%   dtws -- daily TWS time series
%
%   Detang Zhong, CCRS, NRCan 2020
%
nf = length(dayfiles);
dtws = zeros(row, col, nf);
mtws = zeros(row, col)*nodata;
mstd = zeros(row, col)*nodata;

%% read the daily EALCO TWS data
for i = 1:nf
    dfile = char(dayfiles(i));
    fid = fopen(dfile);
        dtws(:,:,i) = (fread(fid,[col, row], data_type))';
    fclose(fid);
end

%% calculate the mean and std of each grid cell
for i = 1:row
    for j = 1:col
        %% remove possible invalid values
        dij = dtws(i,j,:);
        %idx = find(dij == -32760 | dij == 9999);
        idx = find(abs(dij) >= 9999);
        if~isempty(idx)
            dij(idx) = [];
        end
        if ~isempty(dij)
            mtws(i,j) = mean(dij);
            mstd(i,j) = std(dij);
        end
    end
end

end

