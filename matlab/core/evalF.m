function [ r ] = evalF(Uout, d, xdiscr, debdiscr, problem, T)
%EVALF Summary of this function goes here
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
deltat = T(2)-T(1);
r= 1/2*deltaxp*nu*d*d';
for i2 = 1:(length(T))
    r = r + deltat*1/2*(Uout(i2,:)*Uout(i2,:)'*deltax);
end
end

