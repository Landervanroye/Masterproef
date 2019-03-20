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


%%% convergeren gradient en oplossing?
plijst = [100,200,400,800,1600,3200,6400,12800, 25600,51200];
deltatlijst = 0.01*[32,16, 8, 4, 2, 1, 0.5, 0.25,0.125];

plijst = [51200];
deltatlijst = 0.01*[128,64,32,16,8,2,1];


solcel = cell(length(plijst), length(deltatlijst));

for pi = 1:length(plijst)
for ti = 1:length(deltatlijst)
fprintf('N %4.0f, deltat, %4.10f\n', plijst(pi), deltatlijst(ti))
deltat = deltatlijst(ti);
T = [0:deltatlijst(ti):10];
MC_discr = struct();
MC_discr.Np = plijst(pi);
MC_discr.deltax = deltax;
MC_discr.x = xdiscr.x;
MC_discr.nodes = linspace(0,L,N+1);
MC_discr.Nxint = xdiscr.N;
xim = randn(length(T), MC_discr.Np);
[samples_beg, weights_beg] =init_MC(problem,MC_discr);
[Uout_MC, Xout_MC, Weights] = simulate_MC(T,deltat,problem,debdiscr, MC_discr,db, xim,samples_beg, weights_beg);
solcel{pi, ti} = Uout_MC(end,:);
end 
end

error = zeros(length(plijst),length(deltatlijst));
for pi = 1:length(plijst)
for ti = 1:length(deltatlijst)
error(pi, ti) = norm(solcel{pi,ti}(end,:)- Uout(end,:));
end 
end



