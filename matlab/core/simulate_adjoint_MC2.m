function [ DJ ] = simulate_adjoint_MC2(T, U,W,X,D,deltax,deltaxp,nu)
%SIMULATE_ADJOINT_MC Summary of this function goes here
%   Detailed explanation goes here
deltat = T(2)-T(1);
NT = length(T);
nbp = size(X,2);
lambda = zeros(1,nbp); %%% 
DJ = zeros(size(D));

basin_n = floor(X(NT, :)/deltax)+1;
dbasin = floor(X(NT, :)/deltaxp)+1;
dbasin_nmin1 = floor(X(NT-1, :)/deltaxp)+1;
lambda(:) = deltat*deltax*U(NT,basin_n); % lambda N-1
for i = 1:nbp
    DJ(dbasin_nmin1(i)) = DJ(dbasin_nmin1(i)) - deltat*lambda(i)*W(NT-1,i)*exp(-D(dbasin_nmin1(i))*deltat);
end

for k = NT-1:-1:2
    basin_n = floor(X(k, :)/deltax)+1; % kan basin_n-1 in volgende iteratie worden!!
    %basin_nmin1 = floor(X(k-1, :)/deltax)+1; %
    dbasin = floor(X(k-1, :)/deltaxp)+1;
    dbasin_nmin1 = floor(X(k-1, :)/deltaxp)+1;
    lambda= deltat*deltax*U(k,basin_n)+ exp(-D(dbasin)*deltat).*lambda;
    for i = 1:nbp
        DJ(dbasin_nmin1(i)) = DJ(dbasin_nmin1(i)) - deltat*lambda(i)*W(k-1,i)*exp(-D(dbasin_nmin1(i))*deltat);
    end
end
DJ =DJ+ nu*deltaxp*D;
end

