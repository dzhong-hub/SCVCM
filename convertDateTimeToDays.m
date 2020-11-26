function [days] = convertDateTimeToDays(dt_str)
%
% Function to convert a date time string list to a days list
%
% Input:
%   dt_str -- a date time string list in the format: *_yyyy_ddd or 
%            input monthly file name in the format *_yyyy-mm-dd_yyyy-mm-dd,
%            example SWC+SWE_2002-04-05_2002-04-30
%
% Output:
%   days  -- total number of days from 2002_01_01
%

%% retreive the year and days from the file name
n = length(dt_str);
days = zeros(n,1);

for i = 1:n
    % DateString = dt_str{i};   % for Ontario
    % formatIn = 'yyyy/mm/dd';  % for Ontario
    % [yyyy,mm,dd,H,M,S] = datevec(DateString,formatIn);
    
    DateString = dt_str(i);     % fo BC
    %formatIn = 'dd-mm-yyyy';    % fo BC
    [yyyy,mm,dd,H,M,S] = datevec(DateString);
    t2 = datetime(yyyy,mm,dd);
    % t2 = datetime(yyyy,1,1);
    t1 = datetime(2002,1,1);
    dt = between(t1,t2,'Days');
    days(i) = caldays(dt);
end
end

