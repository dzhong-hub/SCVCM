function [NumDays, Time_Str] = daysFromFileName(fname, daily)
%
% Function to determine the time stamp from the input file name
%
% Input:
%   fname -- input daily file name in the format: *_yyyy_ddd or 
%            input monthly file name in the format *_yyyy-mm-dd_yyyy-mm-dd,
%            example SWC+SWE_2002-04-05_2002-04-30
%   daily -- indicator for daily files 1: daily 0: monthly
%
% Output:
%   NumDays  -- total number of days from 2002_01_01
%   Time_Str -- time string
%

%% retreive the year and days from the file name
[path, name,ext] = fileparts(fname);

idx = length(name);

switch daily
    case 1
        yyyy = str2double(name((idx-7):(idx-4)));
        ddd = str2double(name((idx-2):idx));
        t2 = datetime(yyyy,1,1);
        t1 = datetime(2002,1,1);
        dt = between(t1,t2,'Days');
        NumDays = caldays(dt) + ddd;
        Time_Str = name((idx-7):idx);
    case 2
        yyyy1 = str2double(name((idx-20):(idx-17)));
        mm1 = str2double(name((idx-15):(idx-14)));
        dd1 = str2double(name((idx-12):(idx-11)));
        yyyy2 = str2double(name((idx-9):(idx-6)));
        mm2 = str2double(name((idx-4):(idx-3)));
        dd2 = str2double(name((idx-1):(idx)));
        t2 = datetime(yyyy2,mm2,dd2);
        t1 = datetime(yyyy1,mm1,dd1);
        t0 = datetime(2002,1,1);
        dt1 = between(t0,t1,'Days');
        dt2 = between(t0,t2,'Days');
        NumDays = (caldays(dt1) + caldays(dt2))/2;
        Time_Str = name((idx-20):(idx));
    case 3       
        yyyy1 = str2double(name((idx-16):(idx-13)));
        dd1 = str2double(name((idx-11):(idx-9)));       
        yyyy2 = str2double(name((idx-7):(idx-4)));
        dd2 = str2double(name((idx-2):(idx)));
        t2 = datetime(yyyy2,1,1);
        t1 = datetime(yyyy1,1,1);
        t0 = datetime(2002,1,1);
        dt1 = between(t0,t1,'Days');
        dt2 = between(t0,t2,'Days');
        NumDays = (caldays(dt1) + dd1 + caldays(dt2) + dd2)/2;
        Time_Str = name((idx-16):(idx));
end
end

