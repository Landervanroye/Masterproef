function [ Uout ] = simulate( T,problem,debdiscr, xdiscr,d)
%EVALF Summary of this function goes here
%   Detailed explanation goes here
beginvwdn= problem.beginvwdn;
dudt = @(t,u) diffusion_ODE(t,u,problem,debdiscr, xdiscr,d);
Uout = ode1(dudt, T, beginvwdn);
end

