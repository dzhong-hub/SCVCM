function [dt] = convertDaysToDateTime(NoOfDays)
%
% Function to convert a day number from 2002-01-01 to date time string
%
% Input:
%   NoOfDays  -- total number of days from 2002_01_01
%
% Output:
%   dt -- a date time string list in the format: *_yyyy_ddd or 
%            input monthly file name in the format *_yyyy-mm-dd_yyyy-mm-dd,
%            example SWC+SWE_2002-04-05_2002-04-30
%

dt = datetime(2002,1,1)+days(NoOfDays);
dt = datetime(dt,'Format','yyyy-MM-dd');

