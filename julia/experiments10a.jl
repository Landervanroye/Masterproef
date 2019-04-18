using MAT
# in dit experiment wordt de invloed van het aantal deeltjes bekeken (deltat = 0.01)
buckets = parse(Int64,ARGS[1]);
print(buckets)


include("core.jl")
L = 1.0;
N = buckets;
deltax = L/N;
alpha = 0.01;
nu = 0.5;
Nd = 10;
deltaxd = L/Nd;
db = Array(range(0.0, stop=0.1, length=Nd));
xdiscrx = Array(range(deltax/2,stop = (L-deltax/2),length=N));
BEGINVWDN = cos.((xdiscrx*2*pi/L)).+1.1;
problem = problem_obj(nu, alpha, BEGINVWDN);
debdiscr = debdiscr_obj(deltaxd, Nd, Array(range(deltaxd/2,stop = (L-deltaxd/2),length=Nd)));
xdiscr = xdiscr_obj(deltax, N, Array(range(deltax/2, stop = L-deltax/2, length = N)));


v = range(2,stop=7,length=50);
plist = [10^(i) for i in v];
solsave = zeros(50,10);


for i in 1:length(plist)
    print(i,"   p = ", plist[i], "\n")
    plisti = plist[i]
    deltat= 0.01;
    T = Array(0:0.01:1)
    MC_discr = MC_discr_obj(Int64(floor(plisti)), deltax, xdiscr.x, Array(range(0,stop=L, length = N+1)), xdiscr.N)
    for j = 1:50
        samples_beg, weights_beg =init_MC(problem,MC_discr);
        rng = MersenneTwister(1234+100*i+j);
        Uout_MC2, samples_end, weights_end= simulate_MC_rng(T,deltat,problem,debdiscr, MC_discr,db, rng,samples_beg, weights_beg);
        sol = Uout_MC2[end,:];
        solsave[j,:] = sol;
    end
    file = matopen(string("exp10res/b_", buckets, "solp", i, ".mat"), "w")
    write(file, "solsave", solsave)
    close(file)
end
