function [ particles, weights ] = simulationstep( particles,weights,dt,deltaQx, alpha,lambda, xi)
%SIMULATIONSTEP Summary of this function goes here
%   Detailed explanation goes here

%% aanpassen weights
for i = 1:length(particles)
    particle = particles(i);
    xpos = floor(particle/deltaQx)+1;
    weights(i) = weights(i)*exp(-lambda(xpos)*dt);
end
% dit stukje kan efficienter
%particles = particles + sqrt(2*dt*alpha)*randn(size(particles));
particles = particles + sqrt(2*dt*alpha)*xi;
% randvwdn
for i =1:length(particles)
    if particles(i) < 0
        particles(i) = -particles(i);
    end
    if particles(i) > 1
        particles(i) = 1-(particles(i)-1);
    end
end

end

