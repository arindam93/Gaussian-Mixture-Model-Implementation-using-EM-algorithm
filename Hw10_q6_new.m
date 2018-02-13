%% EM algorithm to train a Gaussian Mixture Model of 5 components (Arindam Duttagupta)

clc;
close all;
clear all;

% Load the data
load('hw10p6_data.mat');
figure;
scatter(X(1,:),X(2,:));
xlim([-1.5 1.5])
ylim([-0.5 0.5]); hold on;

% define the component distributions
gaussian = @(x,mu,sigma) (1/(2*pi*sqrt(det(sigma)))).*exp(-0.5*(x-mu)'*pinv(sigma)*(x-mu));

% initialise the variables mean, covariance and beta
mu = [-0.4422, 0.5231, -0.03032, 0.2701, 1.1042; 0.2177, 0.0883, 0.00545, -0.02762, -0.07244];

beta = [1/5, 1/5, 1/5, 1/5, 1/5];

sigma =[];
for i = 1:5
    sigma(:,:,i) = rand*eye(2,2);
end


% Number of iterations
iter = 1;
max_iter = 200;
gamma = [];
Gamma = [];

%% EM algorithm
while 1
    
    for i = 1:size(X,2)
        for j = 1:5
            gamma = [gamma beta(j)*gaussian(X(:,i),mu(:,j),sigma(:,:,j))];
        end
        gamma = gamma/sum(gamma);
        Gamma(i,:) = gamma;
        gamma = [];
    end
    
    for i = 1:5
        beta(i) = mean(Gamma(:,i));
        mu(:,i) = sum(([Gamma(:,i) Gamma(:,i)]'.*X)')./sum(Gamma(:,i));
        sigma(:,:,i) = ((X-mu(:,i)).*[Gamma(:,i),Gamma(:,i)]')*(X-mu(:,i))'./sum(Gamma(:,i));      
    end
    
    if iter > max_iter
        break
    else
        iter = iter + 1;
    end
end
    
    
%% Contour plot of each of five mixture components
for i=1:5
    F = @(x,y) gaussian([x;y],mu(:,i),sigma(:,:,i));
    fcontour(F);
    hold on
end









