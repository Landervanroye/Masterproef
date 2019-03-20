function [ f,g ] = optimfuncMC( T, deltat,problem,MC_discr, debdiscr, xdiscr,d, xim, samples_beg, weights_beg)
%OPTIMFUNC Summary of this function goes here
%   Detailed explanation g
[Uout_MC, Xout_MC, Weights] = simulate_MC(T,deltat,problem,debdiscr, MC_discr,d, xim,samples_beg, weights_beg);
f = evalF(Uout_MC, d, xdiscr, debdiscr, problem, T);
g= simulate_adjoint_MC(T, Uout_MC,Weights,Xout_MC,d,MC_discr.deltax,debdiscr.deltax, problem.nu);
end

