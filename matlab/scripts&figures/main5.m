tic
clear all
close all
%%%% problem parameters

xdiscr=struct();
debdiscr = struct();
problem = struct();


L = 1;
N = 50; % N is het aantal gridcells
% neumann geïmplementeerd volgens eindige volume methode,http://www.csc.kth.se/utbildning/kth/kurser/DN2255/ndiff13/Lecture3.pdf,  https://www.math.ubc.ca/~gustaf/M31611/fd.pdf
deltax = L/N;
alpha = 0.01; %alpha^1 in diffusievergelijking!
nu = 0.5;
deltat= 0.01;
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
MC_discr.Np = 1000;
MC_discr.deltax = deltax;
MC_discr.x = xdiscr.x;
MC_discr.nodes = linspace(0,L,N+1);
MC_discr.Nxint = xdiscr.N;
xim = randn(length(T), MC_discr.Np);
[samples_beg, weights_beg] =init_MC(problem,MC_discr);


fund = @(d) optimfunc( T, deltat,problem,debdiscr, xdiscr,d);
funs = @(d, xi) optimfuncMC( T, deltat,problem,MC_discr, debdiscr, xdiscr,d,xi, samples_beg, weights_beg);

% deterministisch optimum
options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true);
doptim = fminunc(fund,db,options);

% stochastisch optimum zelfde getallen
Nstoch = 15;
results = zeros(Nstoch, length(db));
for i = 1:Nstoch
xim = randn(length(T), MC_discr.Np);
options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true);
results(i,:) = fminunc(@(x) funs(x,xim),db,options);
end

% stochastisch optimum andere 
figure
options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true, 'OutputFcn', @myoutput,'MaxIterations',15);
result2 = fminunc(@(x) funs2(x,T,MC_discr,funs),db,options);
hold on
plot(doptim(4), doptim(5), 'k.')

toc
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
