function [scal,sstd] = scaleFactorFromEALCO(ts, idx)
%
% function for determining priori scale factors from daily EALCO TWS
%
% Input:
%   ts -- EALCO TWS time series for a JPL mascon
%   idx -- idx with valid data
%
% Output:
%   scal -- scal factors of each grid cell of EALCO TWS
%   sstd -- standard deviation of the scal factors
%
%   Detang Zhong, CCRS, NRCan 2020
%
[r,c,b] = size(ts);
n = length(idx);
dij = zeros(n,b);
scal = zeros(n,1);
sstd = zeros(n,1);
nout = zeros(n,1);

%% calculate the scale factor and its std of each grid cell
%% calculate mascon mean for each day
mascon_mean = zeros(1,b);     % mascon average
for i = 1:b
    day = ts(:,:,i);
    %% retreive valid data indexed by idx for the mascon only
    %% and then calculate the mean for the mascon average
    dij(:,i) = day(idx);
    %% check for invalid values and exclude them for the mean
    tmp = dij(:,i);
    tmp(abs(tmp)>=9999)=[]; % cover invalid values <-32760 and 9999
    mascon_mean(i) = mean(tmp);
end
%% least squares fit the scale factor of each grid cell 
for i = 1:n
    y = dij(i,:)';
    x = mascon_mean';
    %% remove possible invalid data
    %inval = find(y == -32760 | y == 9999);
    inval = find(abs(y)>= 9999); % cover invalid values <-32760 and 9999
    if~isempty(inval)
        y(inval) = [];
        x(inval) = [];
    end
    m = length(y);
    if m > 1
        k = x\y;
        v = y - x*k;            
        sigma = sqrt(v'*v/(m - 1));
        
        %% remove outliers and estimate the scale factor again
        inval = find(abs(v)>3*sigma);
        if~isempty(inval)
            y(inval) = [];
            x(inval) = []; 
            k = x\y;
            v = y - x*k;            
            sigma = sqrt(v'*v/(m - 1));
        end              
        
        scal(i) = k;
        sstd(i) = sigma;
        nout(i) = length(inval);
    else
        scal(i) = NaN;
        sstd(i) = NaN;
        nout(i) = NaN;
    end
end

end

