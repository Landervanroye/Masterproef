using MAT


include("core.jl")


#v = range(2,stop=7,length=50);
#plist = [10^(i) for i in v];
gradsave = zeros(50,10);

N0 = 1;
Nb0 = 10;
Np0 = 1000;

a = Nb0/Np0^(1.0/2);
b = N0/Np0^(1.0/2);


plist = zeros(1,500);
Nlist = zeros(size(plist));
lrlist = zeros(size(plist));
Nblist = zeros(size(plist));
lr = 0.9;
Np = Np0;

for i = 1:500
    global Np;
    global lr;
    plist[i] = Int64(Np);
    Nlist[i] = Int64(ceil(b*Np^(1.0/2)));
    Nblist[i] = Int64(ceil(a*Np^(1.0/2)));
    lrlist[i] = lr;

    lr = lr;
    Np = ceil(Np*(1.05));
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
steps = 1;
# 1.5 100 1.05 werkt redelijk snel maar wel plateaus
while z < length(plist)
    global z;
    global seed;
    global j;
    global stationary;
    global db;
    global gradv;
    global gempos;
    global steps;
    file = matopen(string("exp29res/poskeep", seed, ".mat"), "w")
    write(file, "poskeep", db)
    close(file)

    info = [z,j,lrlist[z],stationary, steps];
    file = matopen(string("exp29res/info", seed, ".mat"), "w")
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
    file = matopen(string("exp29res/grad", seed, ".mat"), "w")
    write(file, "grad", grad)
    close(file)
    #print(size(db))
    #print(size(grad))
    db = db .- lrlist[z]*grad;


    if j == steps
        #db = gempos;
        file = matopen(string("exp24res/gempos", seed, ".mat"), "w")
        write(file, "gempos", gempos)
        close(file)
        gradv = zeros(size(db));
        gempos = zeros(size(db));
        j = 0;
        z = z+1;
        steps = ceil(steps/1.0^(1.0/2));
    end
    j = j +1;
    gempos = gempos + db/steps;
    seed = seed + 1;
    doptim = [1.05297419495814	0.865937785927921	0.550081595087884	0.246567832768274	0.0913055272180670	0.0913084900456580	0.246568810439528	0.550080718778183	0.865940131211980	1.05297700126303]
    if seed%1 == 0
        err = sqrt(sum((doptim - db).^2));
        print("j  = ", seed, "  err = ", err, "\n");
    end
    # z en j nog updaten!!
end
