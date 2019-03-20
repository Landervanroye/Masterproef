clear all
close all

%% MC
%initialisatie
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






fund = @(d) optimfunc( T, deltat,problem,debdiscr, xdiscr,d);
[det_sol, det_grad] = fund(db);




solcel = cell(1,2);
ximfull = randn(length(T), 20000);

deltat = 0.01;
T = [0:deltat:10];
MC_discr = struct();
MC_discr.Np = 20000;
MC_discr.deltax = deltax;
MC_discr.x = xdiscr.x;
MC_discr.nodes = linspace(0,L,N+1);
MC_discr.Nxint = xdiscr.N;
xim = ximfull;
[samples_begf, weights_begf] =init_MC(problem,MC_discr);


[Uout_MCfull, Xout_MCfull, Weightsfull] = simulate_MC(T,deltat,problem,debdiscr, MC_discr,db, xim,samples_begf, weights_begf);



xim = ximfull(:, 1:10000);
MC_discr.Np = 10000;
samples_beg = samples_begf(1:10000);
weights_beg = weights_begf(1:10000);
[Uout_MC1, Xout_MC1, Weights1] = simulate_MC(T,deltat,problem,debdiscr, MC_discr,db, xim,samples_beg, weights_beg);

xim = ximfull(:, 10001:20000);
MC_discr.Np = 10000;
samples_beg = samples_begf(10001:20000);
weights_beg = weights_begf(10001:20000);
[Uout_MC2, Xout_MC2, Weights2] = simulate_MC(T,deltat,problem,debdiscr, MC_discr,db, xim,samples_beg, weights_beg);


