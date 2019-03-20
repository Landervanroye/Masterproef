clear all
close all
addpath('/home/lander/Masterproef/core')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% time estimation 
%%% T = cte * Npart * Ntijd = 9*10^(-7) * Npart*Ntijd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% EERSTE EXPERIMENT -> VERSCHILLEND AANTAL DEELTJES
plist = floor(logspace(2, 6, 50));
grad_cel = cell(3,length(plist));

for i = 1:length(plist)
    fprintf('index %d\n', i)
    [ps] = create_problem_struct(0.01, 50);
    fund = @(d) optimfunc( ps.T, ps.deltat,ps.problem,ps.debdiscr, ps.xdiscr,ps.db);
    [det_sol, det_grad] = fund(ps.db);
    grad_cel{1,i} = zeros(50,ps.debdiscr.N);
    grad_cel{2,i} = ps;
    grad_cel{3,i} = det_grad;
    for j = 1:50
        grad_cel{1,i}(j,:)  = batched_fwd_bwd(ps.T,ps.deltat,ps.problem,ps.debdiscr, ps.MC_discr,ps.db, 10^4, plist(i));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%load('grad_cel_0_01.mat')
ps = grad_cel{2,1};
stdm = zeros(length(plist),ps.debdiscr.N);
meanm = zeros(length(plist),ps.debdiscr.N);
det_grad_matrix = zeros(length(plist),ps.debdiscr.N); 
for i = 1:length(plist)
    det_grad_matrix(i,:) = grad_cel{3,i};
    %det_grad_matrix(i,:) = det_grad;
    huidig = grad_cel{1,i};
    meanm(i,:) = mean(huidig,1);
    stdm(i,:) = sqrt(var(huidig,1));
end
stdnorm = sqrt(sum(stdm.^2,2));
meanerror = sqrt(sum((meanm-det_grad_matrix).^2,2));
figure
loglog(plist, stdnorm, 'lineWidth', 2)
hold on
loglog(plist, meanerror, 'lineWidth', 2)
loglog(logspace(2,6,10),1./sqrt(logspace(2,6,10)))
loglog(logspace(2,6,10),10./(logspace(2,6,10)))
xlabel('amount of particles')
legend('std', 'bias', '1/sqrt(p)', '1/p')

%%%% EXPERIMENT 2: error als ~ N, deltat, sqrt(buckets) zelfde verhouding
plist = floor(logspace(4, 5, 10));
grad_cel = cell(3,length(plist));

for i = 1:length(plist)
    fprintf('index %d\n', i)
    [ps] = create_problem_struct(0.01*10^4/plist(i), 10*floor(50*sqrt(plist(i)/10^4)/10));
    fund = @(d) optimfunc( ps.T, ps.deltat,ps.problem,ps.debdiscr, ps.xdiscr,ps.db);
    [det_sol, det_grad] = fund(ps.db);
    grad_cel{1,i} = zeros(50,ps.debdiscr.N);
    grad_cel{2,i} = ps;
    grad_cel{3,i} = det_grad;
    for j = 1:50
        grad_cel{1,i}(j,:)  = batched_fwd_bwd(ps.T,ps.deltat,ps.problem,ps.debdiscr, ps.MC_discr,ps.db, 10^5, plist(i));
    end
end
%%%


%%% invloed aantal buckets? -> bias vergelijken voor genoeg deeltjes
%%% invloed deltat t? 0.01->0.001 10 ptn voor genoeg deeltjes

