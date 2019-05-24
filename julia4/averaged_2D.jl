using MAT
include("core.jl")
L = 1.0;
N = 20;
deltax = L/N;
alpha = 0.5;
nu = 0.05;
Nd = 10;
deltaxd = L/Nd;
#db = Array(range(0.0, stop=0.1, length=Nd));
#db = [1.2397    1.0841    0.8199    0.5471    0.3815    0.3815    0.5471    0.8199    1.0841    1.2397];
db = Float64(0.06)*rand(10,10);
xdiscrx = [0,0];
BEGINVWDN = Float64(1.0)*ones(N,N);
problem = problem_obj(nu, alpha, BEGINVWDN);
debdiscr = debdiscr_obj(deltaxd, Nd, [0,0]);
xdiscr = xdiscr_obj(deltax, N, [0,0]);
deltat= 0.01;
T = Array(0:0.02:1);
deeltjes = 10^5;
MC_discr = MC_discr_obj(deeltjes, deltax, xdiscr.x, Array(range(0,stop=L, length = N+1)), xdiscr.N)
lrlist = [10];
for i in 1:length(lrlist)
    lr = lrlist[i];
    poskeep = zeros(5000,10,10);
    dbtemp = copy(db);
    poskeep[1,:,:] = dbtemp';
    print(i, "\n");
    wbeg = 1.0/(deeltjes/(N*N));
    possafe = zeros(10,10);
    for j = 2:5000
        print(j, "\n");
        samples_beg = rand(2,deeltjes);
        weights_beg = wbeg*ones(1,deeltjes);
        rng = MersenneTwister(1234+100000*i+j);
        Uout_MC2, samples_end, weights_end= simulate_MC_rng_2D(T,deltat,problem,debdiscr, MC_discr,dbtemp, rng,samples_beg, weights_beg);
        rng = MersenneTwister(1234+100000*i+j);
        #print(Uout_MC2[end,:,:], "\n")
        grad= simulate_adjoint_MC_rng_alt_2D(T,Uout_MC2,dbtemp,samples_beg, weights_beg, rng,MC_discr.deltax,debdiscr.deltax,problem.nu, problem, debdiscr, MC_discr);
        poskeep[j,:,:] = poskeep[j-1,:,:] .- lr*grad;
        dbtemp[:,:] = poskeep[j,:,:];
        possafe[:,:] = poskeep[j,:,:];
        print(possafe, '\n')
        file = matopen(string("2dres/poskeep", j, ".mat"), "w")
        write(file, "poskeep", possafe)
        close(file)
        end
    end
end
