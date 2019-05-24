using Random

function simulate_MC_rng_2D(T,deltat,problem,debdiscr, MCdiscr,db, randgen, samples_beg, weights_beg)
    Xbeg = MCdiscr.x;
    deltax = MCdiscr.deltax;
    Nxint = MCdiscr.Nxint;

    Uout = zeros(length(T), Nxint, Nxint);
    for s in 1:length(weights_beg)
        xpos = Int64(floor(samples_beg[1,s]/deltax)+1);
        ypos = Int64(floor(samples_beg[2,s]/deltax)+1);
        Uout[1,xpos,ypos] = Uout[1,xpos,ypos] + weights_beg[s];
    end
    samples_prev = copy(samples_beg);
    weights_prev = copy(weights_beg);

    deltaQx = debdiscr.deltax;
    alpha = problem.alpha;
    constant = sqrt(2*deltat*alpha);
    for s in 1:length(weights_prev)
        weight = weights_prev[s];
        particle = samples_prev[:,s];


        for i in 2:length(T)
            xpos = Int64(floor(particle[1]/deltaQx)+1);
            ypos = Int64(floor(particle[2]/deltaQx)+1);

            weight = weight*exp(-db[xpos,ypos]*deltat);
            # posities
            particle = particle + constant*randn(randgen,2);
            # randvwdn
            if particle[1] < 0
                particle[1]= -particle[1];
            end
            if particle[1]> 1
                particle[1] = 1-(particle[1]-1)
            end

            if particle[2] < 0
                particle[2]= -particle[2];
            end
            if particle[2]> 1
                particle[2] = 1-(particle[2]-1)
            end
            #histogram
            xpos = Int64(floor(particle[1]/deltax)+1);
            ypos = Int64(floor(particle[2]/deltax)+1);

            Uout[i,xpos,ypos] = Uout[i,xpos,ypos] + weight;
        end
        weights_prev[s]= weight;
        samples_prev[:,s]= particle;
    end
    Uout, samples_prev, weights_prev
end


function simulate_adjoint_MC_rng_alt_2D(T,U,D,samples_beg::Array{Float64}, weights_beg::Array{Float64}, randgen,deltax::Float64,deltaxp::Float64,nu::Float64, problem, debdiscr, MC_discr)
    deltat = (T[2]-T[1])::Float64;
    NT = length(T);
    nbp = length(weights_beg);
    DJ = zeros(size(D))::Array{Float64};
    Np = MC_discr.Np;
    X = zeros(2,length(T))::Array{Float64};
    W = zeros(length(T))::Array{Float64};

    deltaQx = debdiscr.deltax::Float64;
    alpha = problem.alpha::Float64;
    constant = sqrt(2*deltat*alpha)::Float64;
    for s = 1:nbp
        particle = samples_beg[:,s];
        weight = weights_beg[s]::Float64;


        X[:,1] = particle;
        W[1] = weight;
        for i in 2:length(T)
            xpos = Int64(floor(particle[1]::Float64/deltaQx::Float64)+1);
            ypos = Int64(floor(particle[2]::Float64/deltaQx::Float64)+1);

            weight::Float64 = (weight::Float64)*exp((-D[xpos::Int64,ypos::Int64]::Float64*deltat::Float64)::Float64)::Float64;

            # posities
            particle = particle + constant::Float64*randn(randgen,2);
            # randvwdn
            if particle[1] < 0
                particle[1]= -particle[1];
            end
            if particle[1]> 1
                particle[1] = 1-(particle[1]-1)
            end

            if particle[2] < 0
                particle[2]= -particle[2];
            end
            if particle[2]> 1
                particle[2] = 1-(particle[2]-1)
            end
        end
        # adjoint

        bNTx = Int64(floor(X[1,NT]/deltax)+1);
        bNTy = Int64(floor(X[2,NT]/deltax)+1);

        bdkmx = Int64(floor(X[1,NT-1]/deltaxp)+1);
        bdkmy = Int64(floor(X[2,NT-1]/deltaxp)+1);

        lambda= 0.5*deltat*deltax*deltax*U[NT,bNTx,bNTy];
        DJ[bdkmx, bdkmy] = DJ[bdkmx, bdkmy] - lambda*deltat*W[NT-1]*exp(-D[bdkmx, bdkmy]*deltat);
        for k = NT-1:-1:2
            xkx =X[1,k];
            xky =X[2,k];

            xkmin1x=X[1,k-1];
            xkmin1y=X[2,k-1];

            wkmin1 = W[k-1];
            bkx = Int64(floor(xkx/deltax)+1);
            bky = Int64(floor(xky/deltax)+1);

            bdkmx = Int64(floor(xkmin1x/deltaxp)+1);
            bdkmy = Int64(floor(xkmin1y/deltaxp)+1);

            bdkmppx = Int64(floor(xkx/deltaxp)+1);
            bdkmppy = Int64(floor(xky/deltaxp)+1);

            lambda = deltat::Float64*(deltax^2)::Float64*U[k,bkx,bky]::Float64 + exp(-D[bdkmppx, bdkmppy]::Float64*deltat::Float64)::Float64*lambda::Float64;
            DJ[bdkmx, bdkmy] = DJ[bdkmx, bdkmy]::Float64 - lambda::Float64*deltat::Float64*wkmin1::Float64*exp((-D[bdkmx,bdkmy]::Float64*deltat::Float64)::Float64)::Float64;
        end
    end
    DJ =DJ.+ nu::Float64*(deltaxp^2)::Float64.*D;
    DJ
end

function dotp(a,b)
    sum = 0;
    for i = 1:length(a)
        sum = sum + a[i]*b[i]
    end
    sum
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
