# in dit experiment wordt de SGD met een vaste leanring rate, voor verschillende learning rates. We starten van ongeveer het optimum

using MAT
include("core.jl")
L = 1.0;
N = 100;
deltax = L/N;
alpha = 0.01;
nu = 0.5;
Nd = 10;
deltaxd = L/Nd;
db = Array(range(0.0, stop=0.1, length=Nd));
#db = [1.2397    1.0841    0.8199    0.5471    0.3815    0.3815    0.5471    0.8199    1.0841    1.2397];
db = [0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8];
xdiscrx = Array(range(deltax/2,stop = (L-deltax/2),length=N));
BEGINVWDN = cos.((xdiscrx*2*pi/L)).+1.1;
problem = problem_obj(nu, alpha, BEGINVWDN);
debdiscr = debdiscr_obj(deltaxd, Nd, Array(range(deltaxd/2,stop = (L-deltaxd/2),length=Nd)));
xdiscr = xdiscr_obj(deltax, N, Array(range(deltax/2, stop = L-deltax/2, length = N)));
deltat= 0.01;
T = Array(0:0.01:1)
MC_discr = MC_discr_obj(10^4, deltax, xdiscr.x, Array(range(0,stop=L, length = N+1)), xdiscr.N)

v = range(0,stop=3,length=20);
lrlist = [10^(-i) for i in v];

for i in 1:length(lrlist)
    lr = lrlist[i];
    poskeep = zeros(5000,10);
    poskeep[1,:] = db';
    print(i, "\n");
    for j = 2:5000
        samples_beg, weights_beg =init_MC(problem,MC_discr);

        rng = MersenneTwister(1234+100000*i+j);
        Uout_MC2, samples_end, weights_end= simulate_MC_rng(T,deltat,problem,debdiscr, MC_discr,db, rng,samples_beg, weights_beg);

        rng = MersenneTwister(1234+100000*i+j);
        grad= simulate_adjoint_MC_rng_alt(T,Uout_MC2,db,samples_beg, weights_beg, rng,MC_discr.deltax,debdiscr.deltax,problem.nu, problem, debdiscr, MC_discr);

        poskeep[j,:] = poskeep[j-1] .- lr*grad';
    end

    file = matopen(string("exp5res/poskeep", i, ".mat"), "w")
    write(file, "poskeep", poskeep)
    close(file)
end
