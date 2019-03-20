function [ r ] = Uf( t,Uout,deltat)
%UF Summary of this function goes here
%   Detailed explanation goes here
index = floor(t/deltat)+1;
r=Uout(index,:)';
end

