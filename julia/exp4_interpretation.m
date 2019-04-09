close all
clear all
det_grad = load('grad_nauwkeurig.mat', 'det_grad'); % deltat = 0.0001, deltax = 1/500
%det_grad = load('grad_100.mat', 'det_grad');

det_grad = det_grad.det_grad;


%% toenemend aantal deeltjes T = 0:01:1, buckets = 100,200,400,800

bucketlist = [10,20,50,100,200,400,800];
figure
for bucketi = 1:length(bucketlist)
buckets = bucketlist(bucketi);    

plist = logspace(2,7,50);
stdevlist = zeros(size(plist));
biaslist = zeros(size(plist));
for i = 1:50
    gradmat = load(sprintf('exp4res/b_%igradp%i.mat',buckets, i), 'gradsave');
    gradmat = gradmat.gradsave;
    meangrad = mean(gradmat,1);
    devgrad = sqrt(sum((gradmat-meangrad).^2,1));
    stdevlist(i) = norm(devgrad);
    biaslist(i) = norm(meangrad-det_grad);
end
loglog(plist, stdevlist)
hold on
loglog(plist, biaslist)
end
legend('10 std', '10 bias','20 std', '20 bias','50 std', '50 bias','100 std', '100 bias', '200 std', '200 bias','400 std', '400 bias','800 std', '800 bias')



%% toenemend aantal tijdstappen p = 10^6, buckets = 100,200,400,800

bucketlist = [10,20,50,100,200,400,800];

figure
for bucketi = 1:length(bucketlist)
buckets = bucketlist(bucketi);    

tlist = logspace(2,3,20);
stdevlist = zeros(size(tlist));
biaslist = zeros(size(tlist));
for i = 1:20
    gradmat = load(sprintf('exp4res/b_%igradt%i.mat',buckets, i), 'gradsave');
    gradmat = gradmat.gradsave;
    meangrad = mean(gradmat,1);
    devgrad = sqrt(sum((gradmat-meangrad).^2,1));
    stdevlist(i) = norm(devgrad);
    biaslist(i) = norm(meangrad-det_grad);
end
loglog(tlist, stdevlist)
hold on
loglog(tlist, biaslist)
end
loglog(linspace(100,1000), 350./(linspace(100,1000).^2), 'lineWidth', 2)
legend('10 std', '10 bias','20 std', '20 bias','50 std', '50 bias','100 std', '100 bias', '200 std', '200 bias','400 std', '400 bias','800 std', '800 bias')
