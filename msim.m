function ms = msim(x,y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
n = length(x);
d = abs(x-y);
s = abs(x)+abs(y);
sim = d./s;
ms = mean(sim);
end

