function [ f,g ] = optimfunc( T, deltat,problem,debdiscr, xdiscr,d)
Uout =simulate( T,problem,debdiscr, xdiscr,d);
f = evalF(Uout, d, xdiscr, debdiscr, problem, T);
Lout = simulate_adjoint( Uout, T,problem,debdiscr, xdiscr,d);
g = evalgradF( Uout, Lout, d, xdiscr, debdiscr, problem,deltat, T);
end

