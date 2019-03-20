clear all
%J = xt A x + b*x
dimension = 3;

S = randn(dimension);
So = orth(S);
A = eye(dimension);
A(2,2) = 5;
A(3,3) = 7;
A = So*A*So';

%A(2,2) = 4;
%A(3,3) = 5;


b = (1:dimension)';
steps = 100000;
x = zeros(dimension,steps);
x(:,1) = 1*ones(dimension,1);
gradlist = zeros(dimension, size(x,2));
learning_rate = 0.1;
noise = 0.1;
for i = 2:steps
    grad = A*x(:,i-1)+b + noise*randn(dimension,1);
    gradlist(:,i-1) = grad;
    x(:,i) = x(:,i-1) - learning_rate*grad;
    %x(:,i) = x(:,i-1) - learning_rate*grad + 0.01*randn(dimension,1);
end
gradlist(:,end) = A*x(:,end)+b + noise*randn(dimension,1);

figure
plot(x(1,:), x(2,:), 'k.')
hold on
sol = - A\b;
plot(sol(1), sol(2), 'r.')
gemiddelde = mean(x(:,end-1000:end),2);
plot(gemiddelde(1), gemiddelde(2), 'g.')


pointlist = [100:100:10000];
biaslist = zeros(1,length(pointlist));
for i = 1:length(pointlist)
    gemiddelde = mean(x(:,end-pointlist(i):end),2);
    biaslist(i) = norm(sol - gemiddelde);
end
figure
loglog(pointlist, biaslist)



N = 10000;
Alsq = zeros(3*N,9);
gradlsq = zeros(3*N,1);
len = size(x,2);
rij = 0;
for i = len-N+1:len
    xhuidig = x(:,i);
    Alsq(rij+1,1:3) = xhuidig;
    Alsq(rij+2,[2,4,5]) = xhuidig;
    Alsq(rij+3,[3,5,6]) = xhuidig;
    Alsq(rij+1:rij+3,7:9) = eye(3);
    gradlsq(rij+1:rij+3) = gradlist(:,i);
    rij = rij+3;
end

benadering = Alsq \ gradlsq;
Aben = zeros(3);
Aben(1,1:3) = benadering(1:3);
Aben(2,1:3) = benadering([2,4,5]);
Aben(3,1:3) = benadering([3,5,6]);
bben = benadering(7:9);

xben = -Aben \ bben;
gemiddelde = mean(x(:,end-N:end),2);
norm(sol -xben)
norm(sol - gemiddelde)

