close all
clear all

for i = 1:5
    poskeep = load(sprintf('exp5res/poskeep%i.mat', i), 'poskeep');
    poskeep = poskeep.poskeep;
end
