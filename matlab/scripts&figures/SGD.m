function [xlist, gradlist, flist] = SGD( optimfunc, steps, init )
%SGD Summary of this function goes here
%   Detailed explanation goes here
xlist=zeros(length(init), steps);
flist=zeros(1,steps);
gradlist =zeros(length(init), steps-1);
xlist(:,1) = init;
learning_rate = 1;
for i = 2:steps
    [f,grad] = optimfunc(xlist(:,i-1)');
    gradlist(:,i-1) = grad;
    flist(i-1) = f;
    if rem(i,10)==0
    fprintf('step %4.0f, f, %4.2f\n', i, f)
    end
    xlist(:,i) = xlist(:,i-1) - learning_rate*grad';
end
end

