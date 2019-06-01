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
plist = zeros(1,100);
Nlist = zeros(size(plist));
lrlist = zeros(size(plist));
Nblist = zeros(size(plist));
lr = 1.0;
Np = Np0;

for i = 1:20
    global Np;
    global lr;
    plist[i] = Int64(Np);
    Nlist[i] = Int64(ceil(b*Np^(2.0/3)));
    Nblist[i] = Int64(ceil(a*Np^(1.0/3)));
    lrlist[i] = lr;

    lr = lr/1.5;
    Np = ceil((Np0/lr)^(3.0/2));
end

db = [0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8];
poskeep = db;
gradv = zeros(size(db));
j = 1;
stationary = false;
gempos = zeros(size(db));
seed = 1;
z=1;
gradv = zeros(size(db));
while z < length(plist)
    global z;
    global seed;
    global j;
    global stationary;
    global db;
    global gradv;
    global gempos;
    file = matopen(string("exp13res/poskeep", seed, ".mat"), "w")
    write(file, "poskeep", poskeep)
    close(file)

    info = [z,j,lr[z],stationary];
    file = matopen(string("exp13res/info", seed, ".mat"), "w")
    write(file, "info", info)
    close(file)

    i =1;
    L = 1.0;
    buckets = Int64(Nblist[z]);
    N = buckets;
    deltax = L/N;
    alpha = 0.01;
    nu = 0.5;
    Nd = 10;
    deltaxd = L/Nd;
    #db = Array(range(0.0, stop=0.1, length=Nd));
    xdiscrx = Array(range(deltax/2,stop = (L-deltax/2),length=N));
    BEGINVWDN = cos.((xdiscrx*2*pi/L)).+1.1;
    problem = problem_obj(nu, alpha, BEGINVWDN);
    debdiscr = debdiscr_obj(deltaxd, Nd, Array(range(deltaxd/2,stop = (L-deltaxd/2),length=Nd)));
    xdiscr = xdiscr_obj(deltax, N, Array(range(deltax/2, stop = L-deltax/2, length = N)));

    plisti = plist[z]
    deltat= 1.0/Nlist[z];

    T = Array(0.0:deltat:1.0)
    MC_discr = MC_discr_obj(Int64(floor(plisti)), deltax, xdiscr.x, Array(range(0,stop=L, length = N+1)), xdiscr.N)
    samples_beg, weights_beg =init_MC(problem,MC_discr);

    rng = MersenneTwister(1234+seed);
    Uout_MC2, samples_end, weights_end= simulate_MC_rng(T,deltat,problem,debdiscr, MC_discr,db, rng,samples_beg, weights_beg);

    rng = MersenneTwister(1234+seed);
    grad= simulate_adjoint_MC_rng_alt(T,Uout_MC2, db ,samples_beg, weights_beg, rng,MC_discr.deltax,debdiscr.deltax,problem.nu, problem, debdiscr, MC_discr);
    file = matopen(string("exp13res/grad", seed, ".mat"), "w")
    write(file, "grad", grad)
    close(file)
    print(size(db))
    print(size(grad))
    db = db .- lrlist[z]*grad;

    if stationary
        gempos = gempos + db/500;
        j = j+1;
    end

    if !stationary && dotp(gradv, grad) < 0
        j = 0;
        gempos = zeros(size(db));
        stationary = true;
    end
    gradv = grad;
    if stationary && j == 500
        db = gempos;
        file = matopen(string("exp13res/gempos", seed, ".mat"), "w")
        write(file, "gempos", gempos)
        close(file)
        gempos = zeros(size(db));
        stationary = false;
        j = 1;
        z = z+1;
    end

    seed = seed + 1;
    # z en j nog updaten!!
end
