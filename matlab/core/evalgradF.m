function [ dp ] = evalgradF( Uout, Lout, d, xdiscr, debdiscr, problem,deltat, T)
%EVALGRADF Summary of this function goes here
%   Detailed explanation goes here
A = xdiscr.A;
alpha = problem.alpha;
deltax = xdiscr.deltax;
beginvwdn = problem.beginvwdn;
nu = problem.nu;
N = xdiscr.N;
Nd = debdiscr.N;
deltaxp = debdiscr.deltax;
Q = debdiscr.Q;

dp = nu*deltaxp*d;
%%for i = 2:length(T)
for i = 1:(length(T)-1)
    dp = dp -deltat*(Lout(i,:)*diag(Uout(i,:))*Q);
end

% T2 = T(1:2:length(T));
% Lambda= @(t) Lout(find(T2==t,1),:)';
% integrand = @(t) Lambda(t)'*diag(Uf(t,Uout,deltat))*Q + deltaxp*nu*d;
% %% midpoint integration
% integral = zeros(1,Nd);
% for i = 1:(length(T2)-1)
%     integral = integral + (T2(i+1)-T2(i))*integrand(T2(i));
% end
% dp = integral';
end

