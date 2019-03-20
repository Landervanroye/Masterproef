function trapezium(y,x,i)
    r = 0;
    deltax = x[2] - x[1];
    for it in 1:i
        r = r + deltax * y[it];
    end
    r
end

function simulationstep( particles,weights,dt,deltaQx, alpha,lambda, xi)
    # weights
    for i in 1:length(particles)
        particle = particles[i];
        xpos = Int64(floor(particle/deltaQx)+1);
        weights[i] = weights[i]*exp(-lambda[xpos]*dt);
    end
    # posities
    particles = particles .+ sqrt(2*dt*alpha).*xi';
    # randvoorwaarden
    for i in 1:length(particles)
        if particles[i] < 0
            particles[i] = -particles[i];
        end
        if particles[i] > 1
            particles[i] = 1-(particles[i]-1)
        end
    end
    particles, weights
end

function simulate_MC(T,deltat,problem,debdiscr, MCdiscr,db, xim, samples_beg, weights_beg)
    Xbeg = MCdiscr.x;
    deltax = MCdiscr.deltax;
    Nxint = MCdiscr.Nxint;

    Uout = zeros(length(T), length(Xbeg));
    Uout[1,:] = pthistogram(samples_beg, deltax, weights_beg, Nxint);

    Xout =zeros(length(T), length(samples_beg));
    Xout[1,:] = samples_beg;

    Weights =zeros(length(T), length(samples_beg));
    Weights[1,:] = weights_beg;

    samples_prev = samples_beg;
    weights_prev = weights_beg;

    for i in 2:length(T)
        xi = xim[i,:];
        samples_prev, weights_prev = simulationstep(samples_prev, weights_prev, deltat,debdiscr.deltax, problem.alpha, db, xi);
        Uout[i,:]= pthistogram(samples_prev, deltax, weights_prev, Nxint);
        Xout[i,:] = samples_prev;
        Weights[i,:] = weights_prev;
    end
    Uout, Xout, Weights
end

function pthistogram( parlist, deltax, weights,Nxint)
    distr = zeros(1,Nxint);
    for i = 1:length(parlist)
        distri = Int64(floor(parlist[i]/deltax)+1);
        distr[distri] = distr[distri]+ weights[i];
    end
    distr
end

function init_MC(problem,MCdiscr)
    beginvwdn = problem.beginvwdn;
    Xbeg = MCdiscr.x;
    Np = MCdiscr.Np;
    deltax = MCdiscr.deltax;
    Nxint = MCdiscr.Nxint;

    norm = trapezium(beginvwdn, Xbeg,length(Xbeg));
    gewicht = norm/Np/deltax;

    cum = zeros(size(Xbeg));
    for i = 1:length(Xbeg)
        cum[i] = trapezium(beginvwdn,Xbeg, i)/norm;
    end

    samples_beg = zeros(1,Np);
    weights_beg = gewicht.*ones(1,Np);

    for i = 1:Np
        samples_beg[i] = sample(MCdiscr.nodes, cum);
    end
    samples_beg, weights_beg
end

function sample(nodes, cum)
    randi = rand();
    if randi == 1
        samp = nodes[end];
    else
        i = 1;
        while randi>=cum[i]
            i = i+1;
        end
        if i == 1
            cumbeg = 0;
        else
            cumbeg = cum[i-1];
        end
        cumeind = cum[i];
        samp = nodes[i] + (randi - cumbeg)/(cumeind-cumbeg)*(nodes[i+1]- nodes[i]);
        if samp < 0
            print("problem")
        end
    end
    samp
end

function simulate_adjoint_MC(T, U,W,X,D,deltax,deltaxp,nu)
    deltat = T[2]-T[1];
    NT = length(T);
    nbp = size(X,2);
    lambda = zeros(NT, nbp);
    DJ = zeros(size(D));

    for s = 1:nbp
        bNT = Int64(floor(X[NT, s]/deltax)+1);
        bdkm = Int64(floor(X[NT-1, s]/deltaxp)+1);
        lambda[NT, s]= deltat*deltax*U[NT,bNT];
        DJ[bdkm] = DJ[bdkm] - lambda[NT, s]*deltat*W[NT-1,s]*exp(-D[bdkm]*deltat);
        for k = NT-1:-1:2
            bk = Int64(floor(X[k, s]/deltax)+1);
            bdkm = Int64(floor(X[k-1, s]/deltaxp)+1);
            bdkmpp = Int64(floor(X[k, s]/deltaxp)+1);
            lambda[k,s] = deltat*deltax*U[k,bk] + exp(-D[bdkmpp]*deltat)*lambda[k+1,s];
            DJ[bdkm] = DJ[bdkm] - lambda[k,s]*deltat*W[k-1,s]*exp(-D[bdkm]*deltat);
        end
    end
    DJ =DJ.+ nu*deltaxp.*D;
    DJ
end



struct xdiscr_obj
    deltax::Float64
    N::Int64
    x::Array{Float64}
end

struct debdiscr_obj
    deltax::Float64
    N::Int64
    x::Array{Float64}
end

struct problem_obj
    nu::Float64
    alpha::Float64
    beginvwdn::Array{Float64}
end

struct MC_discr_obj
    Np::Int64
    deltax::Float64
    x::Array{Float64}
    nodes::Array{Float64}
    Nxint::Int64
end


L = 1.0;
N = 50;
deltax = L/N;
alpha = 0.01;
nu = 0.5;
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

MC_discr = MC_discr_obj(1000000, deltax, xdiscr.x, Array(range(0,stop=L, length = N+1)), xdiscr.N)
xim = randn(length(T), MC_discr.Np);
samples_beg, weights_beg =init_MC(problem,MC_discr);
Uout_MC, Xout_MC, Weights = simulate_MC(T,deltat,problem,debdiscr, MC_discr,db, xim,samples_beg, weights_beg);
grad =  simulate_adjoint_MC(T, Uout_MC,Weights,Xout_MC,db,MC_discr.deltax,debdiscr.deltax,problem.nu);
