function [ distr ] = pthistogram( parlist, deltax, weights,Nxint)
%PARTTODISTR Summary of this function goes here
%   Detailed explanation goes here
%N = length(x)-1;
distr = zeros(1,Nxint);
for i = 1:length(parlist)
    distri = floor(parlist(i)/deltax)+1;
    distr(distri) = distr(distri)+ weights(i);
end
end

