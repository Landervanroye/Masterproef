%zelfde figuur met minder deeltjes
Nd = var.debdiscr.N;
grad2 = zeros(length(1:2:30), Nd);
for i = 1:10:30
    disp(sprintf('index %d', i))
    gradh = zeros(1,Nd);
    for j = 1:i
        X = load(sprintf('results/%d_Xout.mat',j), 'Xout_MC');
        W = load(sprintf('results/%d_Wout.mat',j), 'Weights');
        gradh = gradh + simulate_adjoint_MC(var.T, Ucell{i},W.Weights,X.Xout_MC,var.db,var.MC_discr.deltax,var.debdiscr.deltax,var.problem.nu);
    end
    grad2((i-1)/10+1, :) = gradh./i;    
end
error_grad2 = sqrt(sum((grad - det_grad).^2,2));