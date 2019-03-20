clear all
close all
finitediff = 0;
optimization = 1;
optimizationMC = 1;
laatstegrafieken =1;


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

Nd = 10; % Nd is het aantal EQUIDISTANTE INTERVALLEN waarover lucht geblazen wordt
deltaxd = L/Nd;
db = linspace(0.0,0.1,Nd); % eerste debiet voor optimalisatie
%db = zeros(1,Nd)
global optimskeep;
optimskeep =zeros(size(db));
global optimdkeep;
optimdkeep =zeros(size(db));


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

%% simulation
T = [0:deltat:10];
%%%%%%%%%%%
T = [0:deltat:10];
Uout = simulate( T,problem,debdiscr, xdiscr,db);

%% resultaat plotten
figure
hold on
for i = 1:100:length(T)
    figure(1)
    plot(linspace(0,1,N)', Uout(i,:))
end

%% adjoint simulatie
Lout = simulate_adjoint( Uout, T,problem,debdiscr, xdiscr,db);

%% integratie ter berekening gradiënten
dp = evalgradF( Uout, Lout, db, xdiscr, debdiscr, problem, deltat, T);
figure
plot(dp, '.')

%% finite differences
if finitediff ~=0
deltap = 0.0000001;
dpfin = zeros(Nd,1);
integralorr= evalF(Uout, db, xdiscr, debdiscr, problem, T);

for i =1:Nd
    daccent = db;
    daccent(i) = daccent(i)+deltap;
    integral = evalF(simulate( T,problem,debdiscr, xdiscr,daccent), daccent, xdiscr, debdiscr, problem, T);
    dpfin(i)= (integral-integralorr)/deltap;
end
figure
hold on
plot(dp)
plot(dpfin)
norm(dp-dpfin')/norm(dp)
end


%% optimalisatieroutine
if optimization ~=0
options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true);
fund = @(d) optimfunc( T, deltat,problem,debdiscr, xdiscr,d);
doptim = fminunc(fund,db,options);
figure
hold on
plot(db)
plot(doptim,'.', 'Markersize', 20)
figure
hold on
Uoutop = simulate( T,problem,debdiscr, xdiscr,doptim);
% plot vaste tijdsintervals
for i = 1:100:length(T)
    plot(linspace(0,L,N)', Uoutop(i,:))
end
end

%% MC
MC_discr = struct();
MC_discr.Np = 10000;
MC_discr.deltax = deltax;
MC_discr.x = xdiscr.x;
MC_discr.nodes = linspace(0,L,N+1);
MC_discr.Nxint = xdiscr.N;
xim = randn(length(T), MC_discr.Np);
[samples_beg, weights_beg] =init_MC(problem,MC_discr);
[Uout_MC, Xout_MC, Weights] = simulate_MC(T,deltat,problem,debdiscr, MC_discr,db, xim, samples_beg, weights_beg);
%% simulate adjoint
DJ  = simulate_adjoint_MC(T, Uout_MC,Weights,Xout_MC,db,MC_discr.deltax,debdiscr.deltax, problem.nu);
figure
plot(DJ,'.')
title('DJ')


%% finite differences
if finitediff ~=0
deltap = 0.00000001;
dpfin = zeros(Nd,1);
integralorr= evalF(Uout_MC, db, xdiscr, debdiscr, problem, T);

for i =1:Nd
    daccent = db;
    daccent(i) = daccent(i)+deltap;
    integral = evalF(simulate_MC(T,deltat,problem,debdiscr, MC_discr,daccent, xim, samples_beg, weights_beg), daccent, xdiscr, debdiscr, problem, T);
    dpfin(i)= (integral-integralorr)/deltap;
end
figure
hold on
plot(DJ)
plot(dpfin)
norm(DJ-dpfin')/norm(DJ);
DJ-dpfin'
end



%% plot MC op vaste tijdsintervals
figure
hold on
for i = 1:100:length(T)
    plot(MC_discr.x, Uout_MC(i,:))
end

%% comparison MC vs method of lines
figure
hold on
for i = 1:200:length(T)
    %plot(MC_discr.x, Uout_MC(i,:), 'LineWidth',2)
    plot(xdiscr.x, Uout(i,:),'-.', 'LineWidth',1)
end
set(gca,'ColorOrderIndex',1)
for i = 1:200:length(T)
    plot(MC_discr.x, Uout_MC(i,:), 'LineWidth',2)
    %plot(xdiscr.x, Uout(i,:),'LineWidth',2)
end
xlabel('x')
ylabel('u')
legend('t=0', 't=2', 't=4','t=6', 't=8', 't=10')

%% calculation of gradients
%% adjoint simulatie det - MC-simulatie
Lout_MC = simulate_adjoint( Uout, T,problem,debdiscr, xdiscr,db);

%% integratie ter berekening gradiënten met MC-simulatie, deterministische achterwaarts
dp_MC = evalgradF( Uout_MC, Lout_MC, db, xdiscr, debdiscr, problem, deltat, T);
figure
plot(dp_MC, '.')
%% comparison
figure
hold on
plot(dp_MC, '.', 'MarkerSize', 20)
plot(dp, '.', 'MarkerSize', 10)
legend('Monte Carlo', 'FWD EU')
title('gradient')

%% optimization routine

if optimizationMC ~=0
options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true);
funs = @(d) optimfuncMC( T, deltat,problem,MC_discr, debdiscr, xdiscr,d,xim, samples_beg, weights_beg);
doptimMC = fminunc(funs,db,options);
figure
hold on
%plot(db)
plot(doptim,'.', 'MarkerSize',20)
plot(doptimMC, '.', 'MarkerSize',20)
legend('det', 'stoch')
figure
hold on
Uoutopmc = simulate( T,problem,debdiscr, xdiscr,doptimMC);
% plot vaste tijdsintervals
for i = 1:100:length(T)
    plot(linspace(0,L,N)', Uoutopmc(i,:))
end
fund(doptim)
funs(doptimMC)

figure
hold on
plot(doptim)
plot(doptimMC)
end

%% grafieken plotten

if laatstegrafieken ~=0
funs = @(d) optimfuncMC( T, deltat,problem,MC_discr, debdiscr, xdiscr,d, xim, samples_beg, weights_beg);
fund = @(d) optimfunc( T, deltat,problem,debdiscr, xdiscr,d);


% deterministische optimalisatieroutine - deterministische berekening f,
% grad
gradd = zeros(size(optimdkeep,1), size(db,2));
funcd = zeros(size(optimdkeep,1), 1);
normd = zeros(size(optimdkeep,1), 1);
for i = 1:size(optimdkeep,1)
    [funcd(i), temp] =  fund(optimdkeep(i,:));
    gradd(i,:) = temp';
    normd(i) = norm(temp);
end

% % stochastische optimalisatieroutine - deterministische berekening f,
% % grad
% gradsd = zeros(size(optimskeep,1), size(db,2));
% funcsd = zeros(size(optimskeep,1), 1);
% normsd = zeros(size(optimskeep,1), 1);
% for i = 1:size(optimskeep,1)
%     [funcsd(i), temp] =  fund(optimskeep(i,:));
%     gradsd(i,:) = temp';
%     normsd(i) = norm(temp);
% end

% stochastische optimalisatieroutine - stochastische berekening f,
% grad
gradss = zeros(size(optimskeep,1), size(db,2));
funcss = zeros(size(optimskeep,1), 1);
normss = zeros(size(optimskeep,1), 1);
for i = 1:size(optimskeep,1)
    [funcss(i), temp] =  funs(optimskeep(i,:));
    gradss(i,:) = temp';
    normss(i) = norm(temp);
end
%stdss = zeros(size(optimskeep,1), 1);
%normss = zeros(size(optimskeep,1), 1);
% for i = 1:size(optimskeep,1)
%     tempfv = zeros(kmax,1);
%     tempgrad = zeros(kmax, size(db,2));
%     for k = 1:kmax
%     [tempfv(k), tempgrad(k,:)] =  funs(optimskeep(i,:));
%     end
%     funcss(i) = mean(tempfv);
%     gradss(i,:) = mean(tempgrad,1);
%     stdss(i,:) = norm(std(tempgrad,1));
%     normss(i) = norm(gradss(i,:));
% end
% figure
%plot(stdss)
figure
hold on
semilogy(funcd)
%semilogy(funcsd)
semilogy(funcss)
legend('det', 'stoch')
title('functiewaarden')


figure
semilogy(normd)
hold on
%semilogy(normsd)
semilogy(normss)
%semilogy(stdss)
legend('det','stoch')
%title('Norm Gradiënt')
xlabel('Optimisation step')
ylabel('Norm gradiënt')

end