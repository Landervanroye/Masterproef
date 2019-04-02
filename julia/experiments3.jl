using MAT
include("core.jl")
L = 1.0;
N = 50;
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


v = range(1,stop=3,length=20);
deltatlist = [10^(i) for i in v];
gradsave = zeros(50,10);

for i in 1:length(deltatlist)
    print(i,"   t = ",  deltatlist[i] , "\n");
    deltat= 1/deltatlist[i];
    T = Array(0:deltat:1);
    MC_discr = MC_discr_obj(10^6, deltax, xdiscr.x, Array(range(0,stop=L, length = N+1)), xdiscr.N)
    for j = 1:50
        samples_beg, weights_beg =init_MC(problem,MC_discr);

        rng = MersenneTwister(1234+100*i+j);
        Uout_MC2, samples_end, weights_end= simulate_MC_rng(T,deltat,problem,debdiscr, MC_discr,db, rng,samples_beg, weights_beg);

        rng = MersenneTwister(1234+100*i+j);
        grad= simulate_adjoint_MC_rng_alt(T,Uout_MC2,db,samples_beg, weights_beg, rng,MC_discr.deltax,debdiscr.deltax,problem.nu, problem, debdiscr, MC_discr);
        gradsave[i,:] = grad;
    end
    file = matopen(string("exp3res/gradt", i, ".mat"), "w")
    write(file, "gradsave", gradsave)
    close(file)
end
