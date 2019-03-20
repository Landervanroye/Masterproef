clear all
close all
%%%% problem parameters

xdiscr=struct();
debdiscr = struct();
problem = struct();


L = 1;
N = 4*100; % N is het aantal gridcells
% neumann geïmplementeerd volgens eindige volume methode,http://www.csc.kth.se/utbildning/kth/kurser/DN2255/ndiff13/Lecture3.pdf,  https://www.math.ubc.ca/~gustaf/M31611/fd.pdf
deltax = L/N;
alpha = 0.01; %alpha^1 in diffusievergelijking!
nu = 0.5;
deltat= 0.0025/16;
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
deltat= 0.001;
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
deltatlijst = 0.01*[16, 8, 4, 2, 1];
gradcell = cell(length(plijst), length(deltatlijst));
solcell = cell(length(plijst), length(deltatlijst));




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
funs = @(d, xi) optimfuncMC( T, deltat,problem,MC_discr, debdiscr, xdiscr,d,xi, samples_beg, weights_beg);
[solcell{pi,ti}, gradcell{pi,ti}] =funs(db, xim);
end 
end
%%%%%% tot hier convergence_test.mat





gradfout = zeros(length(plijst), length(deltatlijst));
solfout = zeros(length(plijst), length(deltatlijst));

for pi = 1:length(plijst)
for ti = 1:length(deltatlijst)
    solfout(pi, ti) = abs(solcell{pi, ti}-det_sol);
    gradfout(pi, ti) = norm(gradcell{pi, ti}- det_grad);
end
end

[Pgrid,Tgrid]= meshgrid(deltatlijst,plijst);
surf(Pgrid, Tgrid, log(gradfout))
legend('P','T', 'error')
% constant Delta t * N
len = min(size(gradfout));
gradl = zeros(1,len);
for i = 1:len
    gradl(i) = gradfout(i,i);
end
    

deltatlijst2 = 0.01*[32,16, 8, 4, 2, 1, 0.5]; % 0.25 -> out of memory

% 16*100 = deltat*sqrt(N) => N = (cte/deltat)
plijst = floor((0.01*0.5*226./deltatlijst2).^2);

% 1*sqr(1000) = 1*30 = deltat*sqrt(N) => N = (cte/deltat)^2
solcel2 = cell(size(deltatlijst2));
gradcel2 = cell(size(deltatlijst2));

for ti = 1:length(deltatlijst2)
fprintf('N %4.0f, deltat, %4.10f\n', plijst(ti), deltatlijst2(ti))
deltat = deltatlijst2(ti);
T = [0:deltatlijst2(ti):10];
MC_discr = struct();
MC_discr.Np = plijst(ti);
MC_discr.deltax = deltax;
MC_discr.x = xdiscr.x;
MC_discr.nodes = linspace(0,L,N+1);
MC_discr.Nxint = xdiscr.N;
xim = randn(length(T), MC_discr.Np);
[samples_beg, weights_beg] =init_MC(problem,MC_discr);
funs = @(d, xi) optimfuncMC( T, deltat,problem,MC_discr, debdiscr, xdiscr,d,xi, samples_beg, weights_beg);
[solcel2{ti}, gradcel2{ti}] =funs(db, xim);
end

gradfout2 = zeros(size(deltatlijst2));
for ti = 1:length(deltatlijst2)
    gradfout2(ti) = norm(gradcel2{ti}- det_grad);
end

deltatlijst3 = 0.01*[32,16, 8, 4, 2, 1, 0.5]; % 0.25 -> out of memory

% 16*100 = deltat*sqrt(N) => N = (cte/deltat)
plijst = (0.01*0.5*51200./deltatlijst3);
solcel3 = cell(size(deltatlijst3));
gradcel3 = cell(size(deltatlijst3));

for ti = 1:length(deltatlijst3)
fprintf('N %4.0f, deltat, %4.10f\n', plijst(ti), deltatlijst3(ti))
T = [0:deltatlijst3(ti):10];
deltat = deltatlijst3(ti);
MC_discr = struct();
MC_discr.Np = plijst(ti);
MC_discr.deltax = deltax;
MC_discr.x = xdiscr.x;
MC_discr.nodes = linspace(0,L,N+1);
MC_discr.Nxint = xdiscr.N;
xim = randn(length(T), MC_discr.Np);
[samples_beg, weights_beg] =init_MC(problem,MC_discr);
funs = @(d, xi) optimfuncMC( T, deltat,problem,MC_discr, debdiscr, xdiscr,d,xi, samples_beg, weights_beg);
[solcel3{ti}, gradcel3{ti}] =funs(db, xim);
end 
gradfout3 = zeros(size(deltatlijst3));
for ti = 1:length(deltatlijst3)
    gradfout3(ti) = norm(gradcel3{ti}- det_grad);
end

figure 
hold on
loglog(deltatlijst2, gradfout2)
loglog(deltatlijst3, gradfout3)
legend('fout grad sqrt(N)*deltat=cte', 'fout N*deltat = cte')

%%%% hoe gedraagt de ruisterm zich ten opzichte van de bias

deltatlijst4 = 0.01*[32,16, 8, 4, 2, 1, 0.5]; % 0.25 -> out of memory

% 16*100 = deltat*sqrt(N) => N = (cte/deltat)
plijst = (0.01*0.5*51200./deltatlijst4);

iterations = 10;
bias = zeros(1,length(deltatlijst));
variance = zeros(1,length(deltatlijst));
for ti = 1:length(deltatlijst4)
fprintf('N %4.0f, deltat, %4.10f\n', plijst(ti), deltatlijst4(ti))
deltat = deltatlijst4(ti);
T = [0:deltatlijst4(ti):10];
MC_discr = struct();
MC_discr.Np = plijst(ti);
MC_discr.deltax = deltax;
MC_discr.x = xdiscr.x;
MC_discr.nodes = linspace(0,L,N+1);
MC_discr.Nxint = xdiscr.N;
sol = zeros(iterations,1);
grad = zeros(iterations,length(db));
for it = 1:iterations
xim = randn(length(T), MC_discr.Np);
[samples_beg, weights_beg] =init_MC(problem,MC_discr);
funs = @(d, xi) optimfuncMC( T, deltat,problem,MC_discr, debdiscr, xdiscr,d,xi, samples_beg, weights_beg);
[sol(it), grad(it,:)] =funs(db, xim);
end
variance(ti) =norm(sum((grad - det_grad).^2,2));
bias(ti) =norm(mean(grad,2) - det_grad);
end 

figure
hold on
plot(deltatlijst4, sqrt(variance))
plot(deltatlijst4, bias)
legend('std', 'bias')