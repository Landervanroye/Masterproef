function [ r ] = trapezium(y,x,i)
%PARALLELLOGRAM Summary of this function goes here
%   Detailed explanation goes here
r=0;
deltax = x(2)-x(1);
for it = 1:i
    r = r + deltax*y(it);
end
end


