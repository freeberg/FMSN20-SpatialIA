%OLS Ordinary Least Squares

sR_red = swissRain(swissRain(:,5)==0,:);
sR_mred = swissRain(swissRain(:,5)==1,:);
Y = sqrt(sR_red(:,1));
X = [ones(size(sR_red(:,2),1),1) sR_red(:,2) sR_red(:,2).^2];   % sR_red(:,2:4);
%X = sR_red(:,2:4);
%plot(X,Y,'*')
B = (X'*X)\(X'*Y);
s2 = norm(Y-X*B)^2/(length(Y)-3);

valueGrid = ~isnan(swissElevation);
elevGrid = swissElevation(valueGrid);

X_hat = [ones(size(elevGrid,1),1) elevGrid(:,1) elevGrid(:,1).^2];% 
%X_hat = [elevGrid(:,1) swissX(valueGrid) swissY(valueGrid)];

Y_hat = nan(size(swissElevation));

e = normpdf(X_hat,0,s2);
Y_hat(valueGrid) = X_hat*B + (e(:,1) + e(:,2) + e(:,3));

imagesc([0 max(swissX(:))], [0 max(swissY(:))], Y_hat, ...
        'alphadata', ~isnan(Y_hat))
hold on
plot(swissBorder(:,1), swissBorder(:,2),'k')
scatter(sR_red(:,3), sR_red(:,4), 20, Y, 'filled','markeredgecolor','r')
axis xy tight; hold off; colorbar
title('sqrt of rainfall and predictions')

%%
%Kriging

D = distance_matrix([sR_red(:, 3), sR_red(:, 4)]);
sz = size(swissElevation);  
coords_all = [sR_red(:, 3:4); sR_mred(:, 3:4); swissGrid(:, 2:3)];
D_all = distance_matrix(coords_all);
I_obs = logical([ones(size(sR_red, 1), 1); zeros(size(sR_mred, 1), 1); zeros(size(swissGrid, 1), 1)]);
I_uobs = logical([zeros(size(sR_red, 1), 1); zeros(size(sR_mred, 1), 1); ones(size(swissGrid, 1), 1)]);

z = Y - X*B;

%%
% Spatial dependecy
Kmax = 90;
[rhat,s2hat,m,n,d,varioest]=covest_nonparametric(D,z,Kmax);
figure
plot(d,rhat,'-',0,s2hat,'o')
hold on

%%
% Bootstrap
% Loop 100 times
% Take a randomperm from Y
rhat_all = [];
for i=1:100
    randz = z(randperm(size(Y,1)));
    Kmax = 90;
    [rhat,s2hat,m,n,d,varioest]=covest_nonparametric(D,randz,Kmax);
    rhat_all = [rhat_all; rhat];
end
rhat_rand = quantile(rhat_all, 4);
plot(d,rhat_rand,'-',0,s2hat,'o')

%%
par1 = covest_ml(D,z,'cauchy');
Sigma_cauchy = cauchy_covariance(D_all,par1(1),par1(2),par1(3));
par2 = covest_ml(D,z,'matern');
Sigma_matern = matern_covariance(D_all,par2(1),par2(2),par2(3));
par3 = covest_ml(D,z,'gaussian');
Sigma_gaussian = gaussian_covariance(D_all,par3(1),par3(2));
%%

%add nugget to the covariance matrix
Sigma_yy = Sigma_gaussian + par3(3)^2*eye(size(Sigma_gaussian));
%and divide into observed/unobserved
Sigma_uu = Sigma_yy(I_uobs, I_uobs);
Sigma_uo = Sigma_yy(I_uobs, I_obs);
Sigma_ou = Sigma_yy(I_uobs, I_obs)';
Sigma_oo = Sigma_yy(I_obs, I_obs);

y_k = Y; 
Beta = ((ones(size(y_k,1),1)')*(Sigma_oo \ ones(size(y_k,1),1)))\(ones(size(y_k,1),1)'*(Sigma_oo \ y_k));
y_u = ones(size(I_obs(I_uobs),1),1) * Beta + Sigma_uo*(Sigma_oo \ (y_k - ones(size(y_k,1),1)*Beta));

y_u = y_u(1:4440);
Y_hat(~isnan(swissElevation)) = y_u;
figure
imagesc([0 max(swissX(:))], [0 max(swissY(:))], Y_hat, ...
        'alphadata', ~isnan(Y_hat))
hold on
plot(swissBorder(:,1), swissBorder(:,2),'k')
scatter(sR_red(:,3), sR_red(:,4), 20, Y, 'filled','markeredgecolor','r')
axis xy tight; hold off; colorbar
title('sqrt of rainfall and predictions')


%%
% Predictions interval
x_u = swissGrid;

% Ordinary Least Squares
var_B = s2 * (X'*X) \ eye(3,3);
var_yhat = x_u * var_B * x_u';

% Kriging
var_yu = Sigma_uu - (Sigma_uo / Sigma_oo) * Sigma_ou + (x_u' - (X' / Sigma_oo) * Sigma_ou)' ...
            / ((X' / Sigma_oo) * X) * (x_u' - (X' / Sigma_oo) * Sigma_ou);
%%
% Validation points comparision

X_val = [ones(size(sR_mred(:,2))) sR_mred(:,2) sR_mred(:,2).^2];
Y_val = sR_mred(:,1);
% Ols
e = normpdf(X_val,0,s2);
Y_val_ols = X_val * B + (e(:,1) + e(:,2) + e(:,3));

% Kriging
I_val = logical([zeros(size(sR_red, 1), 1); ones(size(sR_mred, 1), 1); zeros(size(swissGrid, 1), 1)]);
Y_val_krig = y_u(I_val);

% Errors
err_ols = abs(sR_mred(:,1) - Y_val_ols)
err_krig = abs(sR_mred(:,1) - Y_val_krig)
mean_ols = mean(err_ols)
mean_krig = mean(err_krig)
