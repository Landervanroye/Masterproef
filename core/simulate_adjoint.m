function [ Lout ] = simulate_adjoint( Uout, T,problem,debdiscr, xdiscr,d)
%ADJUNCT Summary of this function goes here
%   Detailed explanation goes here
A = xdiscr.A;
alpha = problem.alpha;
deltax = xdiscr.deltax;
deltat = T(2)-T(1);
N = xdiscr.N;
Nd = debdiscr.N;
Q = debdiscr.Q;
Lout = zeros(length(T),size(A,1));
Lout(length(T)-1,:) = deltax*deltat*Uout(end,:);
for i = length(T)-2:-1:1
    Lout(i,:) = Lout(i+1,:)+deltat*deltax*Uout(i+1,:) + deltat*(Lout(i+1,:)*alpha*A - (d*Q').*Lout(i+1,:));
end
end

