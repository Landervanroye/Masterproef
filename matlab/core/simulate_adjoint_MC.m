function [ DJ ] = simulate_adjoint_MC(T, U,W,X,D,deltax,deltaxp,nu)
%SIMULATE_ADJOINT_MC Summary of this function goes here
%   Detailed explanation goes here
deltat = T(2)-T(1);
NT = length(T);
nbp = size(X,2);
lambda = zeros(NT-1, nbp); %%% 
DJ = zeros(size(D));

for s = 1:nbp
    bNT = floor(X(NT, s)/deltax)+1;
    bdkm = floor(X(NT-1, s)/deltaxp)+1;
    lambda(NT, s)= deltat*deltax*U(NT,bNT); %%% negeer de eerste rij (makkelijker met indices, minder vergissingen)
    DJ(bdkm) = DJ(bdkm) - lambda(NT, s)*deltat*W(NT-1,s)*exp(-D(bdkm)*deltat);
    for k = NT-1:-1:2
        bk = floor(X(k, s)/deltax)+1;
        bdkm = floor(X(k-1, s)/deltaxp)+1;
        bdkmpp = floor(X(k, s)/deltaxp)+1;
        lambda(k,s) = deltat*deltax*U(k,bk) + exp(-D(bdkmpp)*deltat)*lambda(k+1,s);
        DJ(bdkm) = DJ(bdkm) - lambda(k,s)*deltat*W(k-1,s)*exp(-D(bdkm)*deltat);    
    end
end
DJ =DJ+ nu*deltaxp*D;

end

