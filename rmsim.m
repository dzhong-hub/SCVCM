function rms = rmsim(x,y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
n = length(x);
d = abs(x-y);
s = abs(x)+abs(y);
sim = d./s;
sim2 = sim.^2;
rms = sqrt(mean(sim2));
end

