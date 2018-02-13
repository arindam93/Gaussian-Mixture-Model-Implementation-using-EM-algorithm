function HW10_Q6()

load('hw10p6_data.mat');
figure(1)
scatter(X(1,:),X(2,:))
hold on
maxIterations = 100;
delta = 10e-6;
[m,n] = size(X);
gaussian = @(x,mu,sigma) (1/(2*pi*sqrt(det(sigma)))).*exp(-0.5*(x-mu)'*inv(sigma)*(x-mu));

% mu = zeros(m,5);
mu = [-0.8,0.125;-0.5,0.05;0,0;0.5,-0.05;0.8,-0.075];
mu = mu';
mu = [-0.9127, -0.5297, 0.03001, 0.5708, 0.9041; 0.1074, 0.08084, 0.005406, -0.02941, -0.06725];
% beta = zeros(1,5);
beta = [1/5,1/5,1/5,1/5,1/5];
sigma = zeros(m,m,5);
constants = [0.1,0.175,0.1,0.125,0.075];
for i=1:1:5
    sigma(:,:,i) = constants(i).*eye(m,m);
end

k=1;
prevlikelihood = 0;
while 1
    Gamma = [];
    for j=1:length(X)
        gamma = [];
        for i=1:1:5
             gamma = [gamma, (beta(i)*gaussian(X(:,j),mu(:,i),sigma(:,:,i)))];
        end
        gamma = gamma/sum(gamma);
        Gamma = [Gamma;gamma];
    end  
    gamma = Gamma;
    size(gamma);
    for i=1:1:5
        beta(i) = mean(gamma(:,i));
        temp = [gamma(:,i),gamma(:,i)];
        mu(:,i) = sum((temp'.*X)')'./sum(gamma(:,i));
        temp_X = ((X-mu(:,i)).*temp')*(X-mu(:,i))';
        sigma(:,:,i) = temp_X./sum(gamma(:,i));
    end
    likelihood = 0;

     if k>maxIterations
         break
     end
     k=k+1;
end
k
mu
sigma

for i=1:1:5
%     G = gmdistribution(mu(:,i)',sigma(:,:,i));
    F = @(x,y) gaussian([x;y],mu(:,i),sigma(:,:,i));
    fcontour(F);
    hold on
end