function [ samples_beg, weights_beg ] = init_MC(problem,MCdiscr)
%INIT_MC Summary of this function goes here
%   Detailed explanation goes here
%SIMULATE_MC Summary of this function goes here
%   Detailed explanation goes here
beginvwdn = problem.beginvwdn;
Xbeg = MCdiscr.x;
Np = MCdiscr.Np;
deltax = MCdiscr.deltax;
Nxint = MCdiscr.Nxint;
%% gewicht per deeltje
norm = trapezium(beginvwdn, Xbeg,length(Xbeg));
gewicht = norm/Np/deltax;
%% cumulatieve distributieberekenen
cum = zeros(size(Xbeg));
for i = 1:length(Xbeg)
    cum(i) = trapezium(beginvwdn,Xbeg, i)/norm;
end
%% sampelen volgens cum distr
samples_beg = zeros(1,Np);
weights_beg = gewicht*ones(1,Np);
for i = 1:Np
    samples_beg(i) = sample(MCdiscr.nodes, cum);
end

end

