#### stel fout ~c1*sqrt(1/N) + c2*(Δt) -> Δt =cte* N^(-1/2) (0.01/2000^(-0.5)) = cte bij 100 buckets zie experiments4)
### zie ook experiments8
cte = (0.01/1000^(-0.5));
using MAT
buckets = 100;
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
gradsave = zeros(50,10);


for i in 46:length(plist)
    print(i,"   p = ", plist[i], "\n")
    plisti = plist[i]
    deltat = (cte*plisti^(-0.5));
    timesteps = floor(1/deltat);
    deltat = 1/timesteps;
    T = Array(0:deltat:1);
    MC_discr = MC_discr_obj(Int64(floor(plisti)), deltax, xdiscr.x, Array(range(0,stop=L, length = N+1)), xdiscr.N)
    for j = 1:50
        samples_beg, weights_beg =init_MC(problem,MC_discr);

        rng = MersenneTwister(1234+100*i+j);
        Uout_MC2, samples_end, weights_end= simulate_MC_rng(T,deltat,problem,debdiscr, MC_discr,db, rng,samples_beg, weights_beg);

        rng = MersenneTwister(1234+100*i+j);
        grad= simulate_adjoint_MC_rng_alt(T,Uout_MC2,db,samples_beg, weights_beg, rng,MC_discr.deltax,debdiscr.deltax,problem.nu, problem, debdiscr, MC_discr);
        gradsave[j,:] = grad;
    end
    file = matopen(string("exp8res/b_", buckets, "gradpt", i, ".mat"), "w")
    write(file, "gradsave", gradsave)
    close(file)
end
