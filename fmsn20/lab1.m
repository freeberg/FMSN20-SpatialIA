R1 = [];
R2 = [];
R3 = [];
sigma2 = 1;
kappa = 1;
nu = 1;
dist = [0:80];
% for i = 0:80
%     R1 = [R1 matern_covariance(i, sigma2, kappa, nu)];
%     R2 = [R1 matern_covariance(i, sigma2, kappa + 1, nu)];
%     R3 = [R1 matern_covariance(i, sigma2, kappa, nu + 1)];
% end

R1 = [R1 matern_covariance(dist, sigma2, kappa, nu)];

figure
plot(R1)
% plot(R2)
% plot(R3)

%%


% First use matern_covariance to create a Sigma-covariance matrix.
% and set mu=(a constant mean).
%sz = [50 60];
[u1,u2] = ndgrid(1:50,1:60);
D = distance_matrix([u1(:), u2(:)]);
mu = 0;
sz = [50 60]
kappa = 0.1;
sigma2 = 1;
Sigma = matern_covariance(D, sigma2, kappa, nu);
R = chol(Sigma); % Calculate the Cholesky factorisation
eta = mu+R'*randn(3000,1); % Simulate a sample
eta_image = reshape(eta,sz);  %reshape the column to an imag
figure
imagesc(eta_image)
%%
sigma_epsilon = 0.2;
y=eta+randn(3000,1)*sigma_epsilon;
z=y-mu;
figure
hold on
plot(D, z*z','.k');
%%
Kmax = 100;
Dmax = max(D(:))
[rhat,s2hat,m,n,d,varioest]=covest_nonparametric(D,z,Kmax,Dmax);
figure
plot(d,rhat,'-',0,s2hat,'o')
hold on
plot(d,varioest,'-')

%%

par=covest_ls(rhat, s2hat, m, n, d, 'matern', [0; 0; 0.2; 0]);
par_cov = par(1:3);
Sigma = matern_covariance(D, par_cov(1), par_cov(2), par_cov(3));
R = chol(Sigma); % Calculate the Cholesky factorisation
eta = mu+R'*randn(3000,1); % Simulate a sample
eta_image = reshape(eta,sz);  %reshape the column to an imag
figure
imagesc(eta_image)
%%

y=eta+randn(3000,1)*par(4);
z=y-mu;
figure
hold on
plot(D, z*z','.k');





