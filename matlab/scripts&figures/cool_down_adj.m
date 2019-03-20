function [ lambda ] = cool_down_adj(t,u,lambda,A,Q,d,alpha, deltax)
%COOL_DOWN_ADJ Summary of this function goes here
%   Detailed explanation goes here
lambda = deltax*u -(alpha*(A')-diag(Q*d'))*lambda;
end

