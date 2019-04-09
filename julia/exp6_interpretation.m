close all
clear all
det_grad = load('grad_nauwkeurig.mat', 'det_grad'); % deltat = 0.0001, deltax = 1/500
%det_grad = load('grad_100.mat', 'det_grad');

det_grad = det_grad.det_grad;


%% toenemend aantal deeltjes T = 0:01:1, buckets = 100,200,400,800

bucketlist = [100];
figure
for bucketi = 1:length(bucketlist)
buckets = bucketlist(bucketi);    

plist = logspace(2,7,50);
stdevlist = zeros(size(plist));
biaslist = zeros(size(plist));
firstgrad = zeros(1,50);
secondgrad = zeros(1,50);
thirdgrad = zeros(1,50);
fourthgrad = zeros(1,50);

for i = 1:5
    gradmat = load(sprintf('exp6res/b_%igradpt%i.mat',buckets, i), 'gradsave');
    gradmat = gradmat.gradsave;
    meangrad = mean(gradmat,1);
    devgrad = sqrt(sum((gradmat-meangrad).^2,1));
    stdevlist(i) = norm(devgrad);
    biaslist(i) = norm(meangrad-det_grad);
end
%loglog(plist, stdevlist)
loglog(plist, biaslist)
hold on
loglog(plist, stdevlist)
xval = logspace(2,7,50);
loglog(xval, 0.4./sqrt(xval))
end
