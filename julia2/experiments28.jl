#v = range(2,stop=7,length=50);
#plist = [10^(i) for i in v];
gradsave = zeros(50,10);

N0 = 5;
Nb0 = 10;
Np0 = 100;

a = Nb0/Np0^(1/3);
b = N0/Np0^(2/3);


plist = zeros(1,50);
Nlist = zeros(size(plist));
lrlist = zeros(size(plist));
Nblist = zeros(size(plist));
lr = 0.4;
Np = Np0;

for i = 1:50
    global Np;
    global lr;
    plist[i] = Int64(Np);
    Nlist[i] = Int64(ceil(b*Np^(2.0/3)));
    Nblist[i] = Int64(ceil(a*Np^(1.0/3)));
    lrlist[i] = lr;

    lr = lr*1.00^(1.0/2.0);
    #Np = ceil(Np*(1.0)^(3.0/2.0));
    Np = ceil(Np*(1.3));
end


using MAT
using Optim
using Random
include("core.jl")

#wat correctheidstesten en snelheidstesten


function f(T,xdiscr, d, seed, deltat, problem, debdiscr, MC_discr, samples_beg, weights_beg)
    rng = MersenneTwister(seed);
    Uout, samples_end, weights_end= simulate_MC_rng(T,deltat,problem,debdiscr, MC_discr,d, rng,samples_beg, weights_beg);
    result = evalF(Uout, d, xdiscr, debdiscr, problem, T);
    return result
end

function g!(T,xdiscr,G, d, seed, deltat, problem, debdiscr, MC_discr, samples_beg, weights_beg)
    rng = MersenneTwister(seed);
    Uout_MC, samples_end, weights_end= simulate_MC_rng(T,deltat,problem,debdiscr, MC_discr,d, rng,samples_beg, weights_beg);
    rng = MersenneTwister(seed);
    G[:]= simulate_adjoint_MC_rng_alt(T,Uout_MC,d,samples_beg, weights_beg, rng,MC_discr.deltax,debdiscr.deltax,problem.nu, problem, debdiscr, MC_discr);
end

optimsave = zeros(50,10);


for i in 1:length(plist)
    print(i,"   p = ", plist[i], "\n")
    plisti = plist[i];
    Nt= Nlist[i];
    N = Int64(Nblist[i]);
    Nd = 10;

    L = 1.0;
    deltax = L/N;
    alpha = 0.01;
    nu = Float64(0.5);
    deltat= 1.0/Nt;
    T = Array(0:deltat:1);
    deltaxd = L/Nd;
    db = Array(range(0.0, stop=0.1, length=Nd));
    xdiscrx = Array(range(deltax/2,stop = (L-deltax/2),length=N));
    BEGINVWDN = cos.((xdiscrx*2*pi/L)).+1.1;
    problem = problem_obj(nu, alpha, BEGINVWDN);
    debdiscr = debdiscr_obj(deltaxd, 10, Array(range(deltaxd/2,stop = (L-deltaxd/2),length=Nd)));
    xdiscr = xdiscr_obj(deltax, N, Array(range(deltax/2, stop = L-deltax/2, length = N)));
    MC_discr = MC_discr_obj(Int64(floor(plisti)), deltax, xdiscr.x, Array(range(0,stop=L, length = N+1)), xdiscr.N)



    for j = 1:50
        samples_beg, weights_beg =init_MC(problem,MC_discr);
        seed = 1234+100*i+j;
        optimres = optimize(d-> f(T,xdiscr,d, seed, deltat, problem, debdiscr, MC_discr, samples_beg, weights_beg), (G,d) -> g!(T,xdiscr,G, d, seed, deltat, problem, debdiscr, MC_discr, samples_beg, weights_beg), db, ConjugateGradient(), Optim.Options(show_trace =false));
        optimsave[j,:] = optimres.minimizer;
    end
    file = matopen(string("exp28res/b_", "optim", i, ".mat"), "w")
    write(file, "optimsave", optimsave)
    close(file)
end
