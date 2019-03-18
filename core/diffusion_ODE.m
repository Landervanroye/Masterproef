function [ u ] = diffusion_ODE(t,u,problem,debdiscr, xdiscr,d)
%COOL_DOWN_DAE Summary of this function goes here
%   Detailed explanation goes here
A = xdiscr.A;
alpha = problem.alpha;
Q = debdiscr.Q;
u =  (alpha*A - diag(Q*d'))*u;
end

