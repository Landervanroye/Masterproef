function [ grad ] = batched_fwd_bwd(T,deltat,problem,debdiscr, MC_discr,db, batch_size, particles)
%BATCHED_FWD_BWD Summary of this function goes here
%   Detailed explanation goes here
MC_discr.Np = batch_size;
U_combined = zeros(length(T),MC_discr.Nxint);

% forward
for i = 1:floor(particles/batch_size)
    xim = randn(length(T), MC_discr.Np);
    [samples_beg, weights_beg] =init_MC(problem,MC_discr);
    [Uout_MC, Xout_MC, Weights] = simulate_MC(T,deltat,problem,debdiscr, MC_discr,db, xim,samples_beg, weights_beg);
    U_combined = U_combined + Uout_MC;
    save(sprintf('results/%d_Uout.mat',i), 'Uout_MC')
    save(sprintf('results/%d_Xout.mat',i), 'Xout_MC')
    save(sprintf('results/%d_Wout.mat',i), 'Weights')
    save(sprintf('results/%d_sb.mat',i), 'samples_beg')
    save(sprintf('results/%d_wb.mat',i), 'weights_beg')
end
rest_N = particles - floor(particles/batch_size)*batch_size;
MC_discr.Np= rest_N;
xim = randn(length(T), MC_discr.Np);
[samples_beg, weights_beg] =init_MC(problem,MC_discr);
[Uout_MC, Xout_MC, Weights] = simulate_MC(T,deltat,problem,debdiscr, MC_discr,db, xim,samples_beg, weights_beg);
U_combined = U_combined + (rest_N/batch_size)*Uout_MC;
save(sprintf('results/rest_Uout.mat'), 'Uout_MC')
save(sprintf('results/rest_Xout.mat'), 'Xout_MC')
save(sprintf('results/rest_Wout.mat'), 'Weights')
save(sprintf('results/rest_sb.mat'), 'samples_beg')
save(sprintf('results/rest_wb.mat'), 'weights_beg')






U_combined = U_combined ./ (particles/batch_size);

%backward
grad = zeros(1,debdiscr.N);
for i = 1:floor(particles/batch_size)
    X = load(sprintf('results/%d_Xout.mat',i), 'Xout_MC');
    W = load(sprintf('results/%d_Wout.mat',i), 'Weights');
    grad = grad + simulate_adjoint_MC(T, U_combined,W.Weights,X.Xout_MC,db,MC_discr.deltax,debdiscr.deltax,problem.nu);
end
MC_discr.Np= rest_N;
X = load(sprintf('results/rest_Xout.mat'), 'Xout_MC');
W = load(sprintf('results/rest_Wout.mat'), 'Weights');
grad = grad + (rest_N/batch_size)*simulate_adjoint_MC(T, U_combined,W.Weights,X.Xout_MC,db,MC_discr.deltax,debdiscr.deltax,problem.nu);



grad = grad ./ (particles/batch_size);
end

