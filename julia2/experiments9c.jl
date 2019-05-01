### experiments on correlated sampling maar met meer deeltjes (zie exp9b)
### nu geen 50 keer opnieuw!!! maar 1 keer opnieuw
using MAT
using Optim
using Random
include("core.jl")

#wat correctheidstesten en snelheidstesten
buckets = 100;
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



function f(d, seed, deltat, problem, debdiscr, MC_discr, samples_beg, weights_beg)
    rng = MersenneTwister(seed);
    Uout, samples_end, weights_end= simulate_MC_rng(T,deltat,problem,debdiscr, MC_discr,d, rng,samples_beg, weights_beg);
    result = evalF(Uout, d, xdiscr, debdiscr, problem, T);
    return result
end

function g!(G, d, seed, deltat, problem, debdiscr, MC_discr, samples_beg, weights_beg)
    rng = MersenneTwister(seed);
    Uout_MC, samples_end, weights_end= simulate_MC_rng(T,deltat,problem,debdiscr, MC_discr,d, rng,samples_beg, weights_beg);
    rng = MersenneTwister(seed);
    G[:]= simulate_adjoint_MC_rng_alt(T,Uout_MC,d,samples_beg, weights_beg, rng,MC_discr.deltax,debdiscr.deltax,problem.nu, problem, debdiscr, MC_discr);
end


v = range(2,stop=7,length=50);
plist = [10^(i) for i in v];
optimsave = zeros(1,10);


for i in 1:length(plist)
    print(i,"   p = ", plist[i], "\n")
    plisti = plist[i]
    deltat= 0.01;
    T = Array(0:0.01:1)
    MC_discr = MC_discr_obj(Int64(floor(plisti)), deltax, xdiscr.x, Array(range(0,stop=L, length = N+1)), xdiscr.N)
    for j = 1:1
        samples_beg, weights_beg =init_MC(problem,MC_discr);
        seed = 1234+100*i+j;
        optimres = optimize(d-> f(d, seed, deltat, problem, debdiscr, MC_discr, samples_beg, weights_beg), (G,d) -> g!(G, d, seed, deltat, problem, debdiscr, MC_discr, samples_beg, weights_beg), db, ConjugateGradient(), Optim.Options(show_trace =false));
        optimsave[j,:] = optimres.minimizer;
    end
    file = matopen(string("exp9cres/b_", buckets, "optim", i, ".mat"), "w")
    write(file, "optimsave", optimsave)
    close(file)
end
