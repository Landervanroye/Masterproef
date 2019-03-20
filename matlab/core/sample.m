function [ samp ] = sample( nodes, cum )
%SAMPLE Summary of this function goes here
%   Detailed explanation goes here
randi = rand(1);
if randi == 1
    samp = nodes(end);
else
i = 1;
while randi>=cum(i)
    i = i+1;
end
if i == 1
    cumbeg = 0;
else 
    cumbeg = cum(i-1);
end
cumeind = cum(i);
samp = nodes(i) + (randi - cumbeg)/(cumeind-cumbeg)*(nodes(i+1)- nodes(i));
if samp < 0
    disp('problem')
end
    
end
end



