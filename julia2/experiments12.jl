using MAT


include("core.jl")


#v = range(2,stop=7,length=50);
#plist = [10^(i) for i in v];
gradsave = zeros(50,10);

N0 = 5;
Nb0 = 10;
Np0 = 100;
a = Nb0/Np0^(1/3);
b = N0/Np0^(2/3);
plist = exp10.(range(2, stop=6, length=50));
Nlist = zeros(size(plist));
Nblist = zeros(size(plist));
for i = 1:length(plist)
    Np = plist[i];
    Nlist[i] = Int64(ceil(b*Np^(2.0/3)));
    Nblist[i] = Int64(ceil(a*Np^(1.0/3)));
end
for z in 1:length(plist)
    i =1;
    L = 1.0;
    buckets = Int64(Nblist[z]);
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



    print(i,"   p = ", plist[z], "\n")
    plisti = plist[z]
    deltat= 1.0/Nlist[z];
    T = Array(0.0:deltat:1.0)
    MC_discr = MC_discr_obj(Int64(floor(plisti)), deltax, xdiscr.x, Array(range(0,stop=L, length = N+1)), xdiscr.N)
    for j = 1:50
        samples_beg, weights_beg =init_MC(problem,MC_discr);

        rng = MersenneTwister(1234+100*z+j);
        Uout_MC2, samples_end, weights_end= simulate_MC_rng(T,deltat,problem,debdiscr, MC_discr,db, rng,samples_beg, weights_beg);

        rng = MersenneTwister(1234+100*z+j);
        grad= simulate_adjoint_MC_rng_alt(T,Uout_MC2, db ,samples_beg, weights_beg, rng,MC_discr.deltax,debdiscr.deltax,problem.nu, problem, debdiscr, MC_discr);
        gradsave[j,:] = grad;
    end
    file = matopen(string("exp12res/gradp", z, ".mat"), "w")
    write(file, "gradsave", gradsave)
    close(file)
end
