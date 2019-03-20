function [ Uout, Xout, Weights ] = simulate_MC(T,deltat,problem,debdiscr, MCdiscr,db, xim, samples_beg, weights_beg)
%SIMULATE_MC Summary of this function goes here
%   Detailed explanation goes here
Xbeg = MCdiscr.x;
deltax = MCdiscr.deltax;
Nxint = MCdiscr.Nxint;
% %% gewicht per deeltje
% norm = trapezium(beginvwdn, Xbeg,length(Xbeg));
% gewicht = norm/Np/deltax;
% %% cumulatieve distributieberekenen
% cum = zeros(size(Xbeg));
% for i = 1:length(Xbeg)
%     cum(i) = trapezium(beginvwdn,Xbeg, i)/norm;
% end
% %% sampelen volgens cum distr
% samples_beg = zeros(1,Np);
% weights_beg = gewicht*ones(1,Np);
% for i = 1:Np
%     samples_beg(i) = sample(MCdiscr.nodes, cum);
% end

%% simulation 
Uout = zeros(length(T), length(Xbeg));
%Uout(1,:) = pthistogram(samples_beg, deltax, weights_beg, Nxint);
Uout(1,:) = pthistogram(samples_beg, deltax, weights_beg, Nxint);


Xout =zeros(length(T), length(samples_beg));
Xout(1,:) = samples_beg;

Weights =zeros(length(T), length(samples_beg));
Weights(1,:) = weights_beg;

samples_prev = samples_beg;
weights_prev = weights_beg;
for i = 2:length(T)
    %[samples_prev, weights_prev] = simulationstep(samples_prev, T(i)-T(i-1), alpha, d.Q,d.Xl);
    xi = xim(i,:);
    [samples_prev, weights_prev] = simulationstep( samples_prev,weights_prev,  deltat,debdiscr.deltax, problem.alpha, db, xi);
    Uout(i,:) = pthistogram(samples_prev, deltax, weights_prev,Nxint);
    Xout(i,:) = samples_prev;
    Weights(i,:) = weights_prev;
end
end

