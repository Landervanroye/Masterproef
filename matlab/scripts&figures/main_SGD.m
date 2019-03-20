clear all
close all
%%%% problem parameters

xdiscr=struct();
debdiscr = struct();
problem = struct();


L = 1;
N = 100; % N is het aantal gridcells
% neumann geïmplementeerd volgens eindige volume methode,http://www.csc.kth.se/utbildning/kth/kurser/DN2255/ndiff13/Lecture3.pdf,  https://www.math.ubc.ca/~gustaf/M31611/fd.pdf
deltax = L/N;
alpha = 0.01; %alpha^1 in diffusievergelijking!
nu = 0.5;
deltat= 0.0025;
T = [0:deltat:10];

Nd = 10; % Nd is het aantal EQUIDISTANTE INTERVALLEN waarover lucht geblazen wordt
deltaxd = L/Nd;
db = linspace(0.0,0.1,Nd); % eerste debiet voor optimalisatie
db = zeros(1,Nd);

problem.nu = nu;
problem.alpha = alpha;

debdiscr.deltax = deltaxd;
debdiscr.N = Nd;
debdiscr.x = linspace(deltaxd/2,L-deltaxd/2,Nd);


xdiscr.deltax= deltax;
xdiscr.N = N;
xdiscr.x = linspace(deltax/2,L-deltax/2,N);

debdiscr.Q = Qdebiet( debdiscr, xdiscr);

%% eindige differenties om spatiele dimensie te discretiseren
A = zeros(N);
for i = 1:N
    A(i,i) = -2/deltax^2;
end
for i = 2:N
    A(i,i-1) = 1/deltax^2;
    A(i-1,i) = 1/deltax^2;
end

%% Neumann BC geïmplementeerd in ODE
A(1,1)=-1/deltax^2;
A(end, end-1) = 1/deltax^2;
A(end, end) = -1/deltax^2;

xdiscr.A = A;
%% beginvoorwaarden
%problem.beginvwdn = ((0.25*L<=xdiscr.x) .* (xdiscr.x<=0.75*L))';
problem.beginvwdn = cos(xdiscr.x*2*pi/L)'+1.1;

%% MC
MC_discr = struct();
MC_discr.deltax = deltax;
MC_discr.x = xdiscr.x;
MC_discr.nodes = linspace(0,L,N+1);
MC_discr.Nxint = xdiscr.N;

MC_discr.Np = 100;
xim = randn(length(T), MC_discr.Np);
[samples_beg, weights_beg] =init_MC(problem,MC_discr);

% deterministisch optimum
fund = @(d) optimfunc( T, deltat,problem,debdiscr, xdiscr,d);
options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true);
doptim = fminunc(fund,db,options);

funs = @(d, xi) optimfuncMC( T, deltat,problem,MC_discr, debdiscr, xdiscr,d,xi, samples_beg, weights_beg);

[result, gradlist, flist] = SGD( @(x) funs2(x,T,MC_discr,funs), 200, db );
figure
plot(result(4,:), result(5,:), 'k.')
hold on
plot(doptim(4), doptim(5), 'r.')
title('4de en 5de element')
figure
semilogy(sqrt(sum(gradlist.^2,1)))
title('norm gradient')
figure
semilogy(flist)
title('f')
norm(mean(result(:,end-50:end),2)-doptim)



%%% variërend aantal deeltjes
Nlist = [100,500,1000,2000,3000,4000,5000,6000,7000,8000];
biaslist = zeros(size(Nlist));
for i = 1:length(Nlist)
    fprintf('N %4.0f, deltat, %4.2f\n', Nlist(i), deltat)
    MC_discr.Np = Nlist(i);
    xim = randn(length(T), MC_discr.Np);
    [samples_beg, weights_beg] =init_MC(problem,MC_discr);
    funs = @(d, xi) optimfuncMC( T, deltat,problem,MC_discr, debdiscr, xdiscr,d,xi, samples_beg, weights_beg);
    [result, gradlist, flist] = SGD( @(x) funs2(x,T,MC_discr,funs), 200, db );
    biaslist(i) = norm(mean(result(:,end-50:end),2)-doptim);
end
figure
semilogy(Nlist,biaslist)
figure
loglog(Nlist,biaslist)

title('bias verschillend aantal deeltjes')


%%% variërende tijdsstap
deltatlist = [0.1,0.05,0.01,0.005,0.001];
biaslisttijd = zeros(size(deltatlist));
MC_discr.Np = 1000;
[samples_beg, weights_beg] =init_MC(problem,MC_discr);

MC_discr.Np = 1000;
xim = randn(length(T), MC_discr.Np);

for i = 1:length(deltatlist)
    fprintf('N %4.0f, deltat, %4.2f\n', N, deltatlist(i))
    deltat= deltatlist(i);
    T = [0:deltat:10];
    
    xim = randn(length(T), MC_discr.Np);
    funs = @(d, xi) optimfuncMC( T, deltat,problem,MC_discr, debdiscr, xdiscr,d,xi, samples_beg, weights_beg);
    [result, gradlist, flist] = SGD( @(x) funs2(x,T,MC_discr,funs), 200, db );
    biaslisttijd(i) = norm(mean(result(:,end-50:end),2)-doptim);
end
figure
semilogy(deltatlist, biaslisttijd)
title('bias verschillende tijdstap')


%veel deeltjes, kleine tijdstap

MC_discr.Np = 10000;
xim = randn(length(T), MC_discr.Np);

deltat= 0.001;
T = [0:deltat:10];

xim = randn(length(T), MC_discr.Np);
funs = @(d, xi) optimfuncMC( T, deltat,problem,MC_discr, debdiscr, xdiscr,d,xi, samples_beg, weights_beg);
[result, gradlist, flist] = SGD( @(x) funs2(x,T,MC_discr,funs), 200, db );
bias = norm(mean(result(:,end-50:end),2)-doptim);


function stop = myoutput(x,optimvalues,state);
stop = false;
state
if isequal(state,'iter')
    hold on
    plot(x(4), x(5), 'r.')
end
end
function [res, grad]=  funs2(d,T,MC_discr, funs)
    xim = randn(length(T), MC_discr.Np);
    [res, grad] = funs(d, xim);
end
