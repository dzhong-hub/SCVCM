function logmsg(flog, handles, msg)
%
% log computation status message into a log file
%
% Input:
%   flog - file id of the message log file
%   msg - msg to be logged
%

%set(handles.texProcessStatus,'String',msg);
msg = [datestr(now) '   ' msg];
disp(msg);
fprintf(flog,'%s\n', msg);