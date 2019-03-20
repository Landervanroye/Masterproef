function [ distr ] = pthistogram( parlist, deltax, weights,Nxint)
%PARTTODISTR Summary of this function goes here
%   Detailed explanation goes here
%%%%%%%%%%%%%% slower than the other one!

distr = zeros(1,Nxint);
distrindices = floor(parlist/deltax)+1;
for i = 1:length(distr)
    %distri = floor(parlist(i)/deltax)+1;
    distr(i) = sum(weights(distrindices==i));
end
end

