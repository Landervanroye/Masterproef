clear all
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% time estimation 
%%% T = cte * Npart * Ntijd = 9*10^(-7) * Npart*Ntijd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        grad_cel{1,i}(j,:)  = batched_fwd_bwd(ps.T,ps.deltat,ps.problem,ps.debdiscr, ps.MC_discr,ps.db, 10^5, plist(i));
    end
end

%%% invloed aantal buckets? -> bias vergelijken voor genoeg deeltjes
%%% invloed deltat t? 0.01->0.001 10 ptn voor genoeg deeltjes

