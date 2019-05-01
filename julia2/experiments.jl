
using Random
include("core.jl")

#wat correctheidstesten en snelheidstesten

L = 1.0;
N = 50;
deltax = L/N;
alpha = 0.01;
nu = Float64(0.5);
deltat= 0.01;
T = Array(0:0.01:1)
Nd = 10;
deltaxd = L/Nd;
db = Array(range(0.0, stop=0.1, length=Nd));
xdiscrx = Array(range(deltax/2,stop = (L-deltax/2),length=N));
BEGINVWDN = cos.((xdiscrx*2*pi/L)).+1.1;
problem = problem_obj(nu, alpha, BEGINVWDN);
debdiscr = debdiscr_obj(deltaxd, Nd, Array(range(deltaxd/2,stop = (L-deltaxd/2),length=Nd)));
xdiscr = xdiscr_obj(deltax, N, Array(range(deltax/2, stop = L-deltax/2, length = N)));
MC_discr = MC_discr_obj(10^5, deltax, xdiscr.x, Array(range(0,stop=L, length = N+1)), xdiscr.N)
rng = MersenneTwister(1234);

## wat gepruts om random getallen beide implementaties overéén te laten komen
xim2 = randn(rng,length(T)-1, MC_discr.Np);
xim = zeros(length(T), MC_discr.Np)
xim[2:end,:]=xim2;


samples_beg, weights_beg =init_MC(problem,MC_discr);


@time Uout_MC, Xout_MC, Weights = simulate_MC(T,deltat,problem,debdiscr, MC_discr,db, xim,samples_beg, weights_beg);
@time grad =  simulate_adjoint_MC(T, Uout_MC,Weights,Xout_MC,db,MC_discr.deltax,debdiscr.deltax,problem.nu);

rng = MersenneTwister(1234);
@time Uout_MC2, samples_end, weights_end= simulate_MC_rng(T,deltat,problem,debdiscr, MC_discr,db, rng,samples_beg, weights_beg);


rng = MersenneTwister(1234);
@time grad2= simulate_adjoint_MC_rng(T,Uout_MC2,db,samples_beg, weights_beg, rng,MC_discr.deltax,debdiscr.deltax,problem.nu, problem, debdiscr, MC_discr);

rng = MersenneTwister(1234);
@time grad2= simulate_adjoint_MC_rng_alt(T,Uout_MC2,db,samples_beg, weights_beg, rng,MC_discr.deltax,debdiscr.deltax,problem.nu, problem, debdiscr, MC_discr);
@time evalF(Uout_MC2, db, xdiscr, debdiscr, problem, T)
#
#MC_discr = MC_discr_obj(10^7, deltax, xdiscr.x, Array(range(0,stop=L, length = N+1)), xdiscr.N)
#samples_beg, weights_beg =init_MC(problem,MC_discr);
#rng = MersenneTwister(1234);
#@time Uout_MC2, samples_end, weights_end= simulate_MC_rng(T,deltat,problem,debdiscr, MC_discr,db, rng,samples_beg, weights_beg);
#rng = MersenneTwister(1234);
#@time grad3= simulate_adjoint_MC_rng(T,Uout_MC2,db,samples_beg, weights_beg, rng,MC_discr.deltax,debdiscr.deltax,problem.nu, problem, debdiscr, MC_discr);
#

#rng = MersenneTwister(1234);
#@time grad2= simulate_adjoint_MC_rng2(T,Uout_MC2,db,samples_end, weights_end, rng,MC_discr.deltax,debdiscr.deltax,problem.nu, problem, debdiscr, MC_discr);


### SNELHEID: simulate_MC_rng 3x zo snel als simulate_MC (omdat gewichten en pos niet moeten worden opgeslagen) -> gewichten en posities opslaan even duur als ze gewoon te berekenen?
### matlab experiment: opslaan op harde schijf ongeveer 7 keer trager dan alles in geheugen houden
### nu 1 miljoen deeltjes -> 3 seconden
