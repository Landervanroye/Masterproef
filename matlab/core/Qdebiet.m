function [ r ] = Qdebiet( debdiscr, xdiscr)
%FDEBIET Summary of this function goes here
%   Detailed explanation goes here
x = xdiscr.x;
deltaxdeb = debdiscr.deltax;
Nxd = debdiscr.N;
r = zeros(length(x), length(debdiscr.x));
for i = 1:(xdiscr.N)
%% uitzoeken in welk deel we ons bevinden
    xi = x(i);
    index = floor(xi/deltaxdeb)+1;
    if index == Nxd + 1
        index = Nxd;
    end
    r(i, index) = 1;
end 
end

