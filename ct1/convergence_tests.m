clear all
close all
%% voorwaartse methode
%%%%% deltat = 0.01
% MC
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
T = [0:deltat:1];


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

% eindige differenties om spatiele dimensie te discretiseren
A = zeros(N);
for i = 1:N
    A(i,i) = -2/deltax^2;
end
for i = 2:N
    A(i,i-1) = 1/deltax^2;
    A(i-1,i) = 1/deltax^2;
end

% Neumann BC geïmplementeerd in ODE
A(1,1)=-1/deltax^2;
A(end, end-1) = 1/deltax^2;
A(end, end) = -1/deltax^2;

xdiscr.A = A;
% beginvoorwaarden
problem.beginvwdn = cos(xdiscr.x*2*pi/L)'+1.1;


fund = @(d) optimfunc( T, deltat,problem,debdiscr, xdiscr,d);
[det_sol, det_grad] = fund(db);

iterations = 1000;

for i = 1:iterations
disp(sprintf('index %d', i))
deltat = 0.01;
T = [0:deltat:1];
MC_discr = struct();
MC_discr.Np = 10000;
MC_discr.deltax = deltax;
MC_discr.x = xdiscr.x;
MC_discr.nodes = linspace(0,L,N+1);
MC_discr.Nxint = xdiscr.N;
xim = randn(length(T), MC_discr.Np);
[samples_beg, weights_beg] =init_MC(problem,MC_discr);
[Uout_MC, Xout_MC, Weights] = simulate_MC(T,deltat,problem,debdiscr, MC_discr,db, xim,samples_beg, weights_beg);
save(sprintf('results/%d_Uout.mat',i), 'Uout_MC')
save(sprintf('results/%d_Xout.mat',i), 'Xout_MC')
save(sprintf('results/%d_Wout.mat',i), 'Weights')
save(sprintf('results/%d_sb.mat',i), 'samples_beg')
save(sprintf('results/%d_wb.mat',i), 'weights_beg')
end
save('results/variables.mat', 'T', 'deltat', 'problem', 'debdiscr', 'MC_discr', 'db', 'det_grad')
Uout_det = simulate( T,problem,debdiscr, xdiscr,db);
save('results/det', 'Uout_det')
save('results/it_data', 'iterations')


%% combineren tot oplossing
clear all
% inladen
det = load('results/det');
var = load('results/variables.mat');
it_data = load('results/it_data');
det_grad = var.det_grad;

% combineren oplossingen
U_combined = zeros(length(var.T),var.MC_discr.Nxint);
Ucell = cell(1,it_data.iterations);
error_sol = zeros(1,it_data.iterations);
for i = 1:it_data.iterations
    Uload = load(sprintf('results/%d_Uout.mat',i));
    U_combined = U_combined + Uload.Uout_MC;
    Ucell{i} = U_combined./i;
    error_sol(i) = norm(Ucell{i}(end,:)-det.Uout_det(end,:));
end



% interpreteren resultaten gecombineerde oplossingen
error_sol = zeros(1,it_data.iterations);
for i = 1:it_data.iterations
    error_sol(i) = norm(Ucell{i}(end,:)-det.Uout_det(end,:));
end
figure
loglog(var.MC_discr.Np*[1:it_data.iterations], error_sol)
hold on
loglog(var.MC_discr.Np*[1:it_data.iterations], 1./sqrt(1:it_data.iterations))
legend('error on solution','1/sqrt(N)')
xlabel('N')
ylabel('absolute error')
title('error on solution vs N, \Delta t = 0.01')
%% achterwaartse methode 
Nd = var.debdiscr.N;
grad = zeros(length(1:5:it_data.iterations), Nd);
for i = 1:5:it_data.iterations
    fprintf('index %d\n', i)
    gradh = zeros(1,Nd);
    for j = 1:i
        X = load(sprintf('results/%d_Xout.mat',j), 'Xout_MC');
        W = load(sprintf('results/%d_Wout.mat',j), 'Weights');
        gradh = gradh + simulate_adjoint_MC(var.T, Ucell{i},W.Weights,X.Xout_MC,var.db,var.MC_discr.deltax,var.debdiscr.deltax,var.problem.nu);
    end
    grad((i-1)/5+1, :) = gradh./i;    
end
error_grad = sqrt(sum((grad - det_grad).^2,2));

figure
loglog(10000*(1:5:it_data.iterations), error_grad)
xlabel('N')
ylabel('absolute error')
title('error on gradient vs N, \Delta t = 0.01')
%zelfde figuur met minder deeltjes
Nd = var.debdiscr.N;
grad2 = zeros(length(1:1:100), Nd);
for i = 1:1:200
    fprintf('index %d\n', i)
    gradh = zeros(1,Nd);
    for j = 1:i
        X = load(sprintf('results/%d_Xout.mat',j), 'Xout_MC');
        W = load(sprintf('results/%d_Wout.mat',j), 'Weights');
        gradh = gradh + simulate_adjoint_MC(var.T, Ucell{i},W.Weights,X.Xout_MC,var.db,var.MC_discr.deltax,var.debdiscr.deltax,var.problem.nu);
    end
    grad2((i-1)/1+1, :) = gradh./i;    
end
error_grad2 = sqrt(sum((grad2 - det_grad).^2,2));



%% achterwaartse methode bias-stochastische fout
plist = 1:2:50;
gradcel = cell(1,length(plist));
Nd = var.debdiscr.N;
grad2 = zeros(length(plist), Nd);
for i = 1:length(plist)
    Ncomb = plist(i);
    gradcel{i} = zeros(10, Nd);
    fprintf('index %d\n', Ncomb)
    for k = 1:20
        gradh = zeros(1,Nd);
        U_combined = zeros(length(var.T),var.MC_discr.Nxint);
        for j = 1:Ncomb
            % combineren oplossingen
            index = 50*(k-1)+j;
            Uload = load(sprintf('results/%d_Uout.mat',index));
            U_combined = U_combined + Uload.Uout_MC;
        end
        U_combined = U_combined./Ncomb;
        for j = 1:Ncomb
            index = 50*(k-1)+j;
            X = load(sprintf('results/%d_Xout.mat',index), 'Xout_MC');
            W = load(sprintf('results/%d_Wout.mat',index), 'Weights');
            gradh = gradh + simulate_adjoint_MC(var.T, U_combined,W.Weights,X.Xout_MC,var.db,var.MC_discr.deltax,var.debdiscr.deltax,var.problem.nu);
        end
        gradcel{i}(k,:) = gradh./Ncomb;
    end
end
stdm = zeros(length(plist),Nd);
meanm = zeros(length(plist),Nd);

for i = 1:length(plist)
    huidig = gradcel{i};
    meanm(i,:) = mean(huidig,1);
    stdm(i,:) = sqrt(sum((huidig-meanm(i,:)).^2,1));
end
meanerror = (meanm-det_grad);
meanerrornorm = sqrt(sum(meanerror.^2,2));
stdnorm = sqrt(sum(stdm.^2,2));

figure
loglog(1/100./sqrt(1:100))
hold on
loglog(meanerrornorm)
loglog(stdnorm)


