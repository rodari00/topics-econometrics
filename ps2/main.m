clear
close all
clc

%% Setup folders
[figures_dir,data_dir,tables_dir] = directorysetup(cd,'Win');
type = 'pdf';

%% PART 0: DATA PRE-PROCESSING 

% =========================================================================
% LOAD AND LABEL DATA

csv_in = strcat(data_dir,'\fredmd.csv');

% Load data from CSV file
dum=importdata(csv_in,',');

% Variable names
series=dum.textdata(1,2:end);

% Transformation numbers
tcode=dum.data(1,:);

% Raw data
rawdata=dum.data(2:end,:);

% Month/year of final observation
final_datevec=datevec(dum.textdata(end,1));
final_month=final_datevec(2);
final_year=final_datevec(1);

% Dates (monthly) are of the form YEAR+MONTH/12
% e.g. March 1970 is represented as 1970+3/12
% Dates go from 1959:01 to final_year:final_month (see above)
dates = (1959+1/12:1/12:final_year+final_month/12)';

% T = number of months in sample
T=size(dates,1);
rawdata=rawdata(1:T,:);


% =========================================================================
% PROCESS DATA

% Transform raw data to be stationary using auxiliary function
% prepare_missing()
yt=prepare_missing(rawdata,tcode);

% Reduce sample to usable dates: remove first two months because some
% series have been first differenced
yt=yt(3:T,:);
dates=dates(3:T,:);

% Remove outliers using auxiliary function remove_outliers(); see function
% or readme.txt for definition of outliers
%   data = matrix of transformed series with outliers removed
%   n = number of outliers removed from each series
[data,n]=remove_outliers(yt);


% =========================================================================
% SELECT DATA 

% Select sample from Jan 1960 to Dec 2019
sample_idx = find(ismember(dates, 1960 + 1/12:1/12:2020));
dates = dates(sample_idx(1):sample_idx(end));
data = data(sample_idx(1):sample_idx(end),:);

% Find index corresponding to December 1989
T0 = find(dates == 1990); 

% Define X and Y matrices
Y = data(:,strcmpi(series,'UNRATE'));
X = data;
X(:,strcmpi(series,'UNRATE')) = [];

%% PART 1.1: AR(p) 

% Information criterion
maxp = 20;
BIC = zeros(maxp,1);
AIC = zeros(maxp,1);
n = size(Y,2);
T = size(Y,1);
for p = 1:20
[beta, eps,~ ] = betaOLS(Y,p,'yes');
SSR = sum(eps.^2);
BIC(p) = log(T^(-1)*(eps'*eps)) + (n^2*p/T)*log(T);
AIC(p) = log(T^(-1)*(eps'*eps)) + 2*(n^2*p/T);
end

%hold on
plot(1:maxp,BIC)
close
%hold off
% chose p = 4
plot(1:maxp,AIC)
close

% =========================================================================
% Compute MSE
[mse_ar pred_ar eps_ar]= MeanSqErr(Y, 4, 1, T0,'yes', 0);

% =========================================================================
% Plot actual vs prediction
hold on
a = plot(dates,Y)
b = plot(dates,[NaN(T0,1);pred_ar])
hold off
xlabel('$t$','interpreter','latex')
ylabel('$\textrm{UNRATE}\; (\upsilon)$','interpreter','latex')
legend([a b], '$\upsilon$', '$\widehat{\upsilon}$', 'interpreter','latex')
close


%% PART 1.2: FACTOR MODEL

% Normalize X
X = X(:,sum(isnan(X),1) == 0);
meanX = mean(X,1);
sdX = std(X,1);
Xnorm = (X-meanX)./sdX;

% Check number of principal components
[V,D] = eig(Xnorm'*Xnorm);
lambda = diag(D);
plot(flipud(lambda(1:20)))
close
r = 20;

% Initialize vector of square errors
eps_pca = zeros(T-T0,1);
pred_pca = zeros(T-T0,1);

for t = 0:T-T0-1
    [V,D] = eig(Xnorm'*Xnorm);
    Xtmp = Xnorm(1:T0+t,:);
    Ytmp = Y(1:T0+t);
    % Factor
    F = Xtmp*V(:,1:r);
    % PCA regression
    betafm = inv(F'*F)*F'*Ytmp;
    % Predict unrate
    pred_pca(t+1) = F(end,:)*betafm;
    % Compute error
    eps_pca(t+1) = (Y(T0+t+1) - pred_pca(t+1))^2;
    
end

% =========================================================================
% Compute MSE
mse_pca = mean(eps_pca)

%% PART 1.3: LASSO

eps_lasso = zeros(T-T0,1);
pred_lasso= zeros(T-T0,1);
k = 10;
% Initialize lasso cross-validated mse
mse_lasso_cv = [];

for t = 0:T-T0 -1
    
    % When to update lambda (5 years)
    if rem(T0 + t, 5*12) == 0
        cv = blockcv([Y(1:T0+t) Xnorm(1:T0+t,:)],k);
        % select grid
        grid_lasso = 0:0.005:0.2;

        % Grid Search over optimal lambda
        [lambda_opt, mse_tmp]= PenaltySearch(cv,grid_lasso,'lasso');

    mse_lasso_cv = [mse_lasso_cv mse_tmp];
    end


    B = lasso(Xnorm(1:T0+t,:),Y(1:T0+t),'Lambda',lambda_opt);

    pred_lasso(t+1) = Xnorm(T0+t,:)*B;
    eps_lasso(t+1) = (Y(T0+t+1) -  pred_lasso(t+1))^2;


end

% =========================================================================
% Compute MSE
mse_lasso = mean(eps_lasso);

%% PART 1.4: RIDGE

eps_ridge = zeros(T-T0,1);
pred_ridge = zeros(T-T0,1);
k = 10;
% Initialize the vector of cross-validated ridges
mse_ridge_cv = [];

for t = 0:T-T0 -1
    
    % When to update lambda (5 years, 5*12)
    if rem(T0 + t, 60) == 0
        cv = blockcv([Y(1:T0+t) Xnorm(1:T0+t,:)],k);
        % select grid
        grid_ridge = 3000:10:10000;

        % Grid Search over optimal lambda
        [lambda_opt, mse_tmp]= PenaltySearch(cv,grid_ridge,'ridge');

        mse_ridge_cv = [mse_ridge_cv mse_tmp];
    end

    B = ridge(Y(1:T0+t), Xnorm(1:T0+t,:),lambda_opt);
    pred_ridge(t+1) = Xnorm(T0+t,:)*B;

    eps_ridge(t+1) = (Y(T0+t+1) - pred_ridge(t+1) )^2;


end

% =========================================================================
% Compute MSE
mse_ridge = mean(eps_ridge);

%% PART 1.5: BAGGING TREES

rng('default')
B = 100;
pred_tree = zeros(T-T0,B);
eps_tree = zeros(T-T0,1);

for t = 0:T-T0 -1
    for b = 1:B
        fprintf('Bootstrap N %d for horizon %d...\n',b,t)
        BootSample = blockbootstrp([Y Xnorm],72,1);
        Yb = BootSample(:,1);
        Xb = BootSample(:,2:end);
        tree = fitrtree(Xb(1:T0+t,:),Yb(1:T0+t));
        pred_tree(t+1,b) = predict(tree,Xb(T0 + t,:));
        
    end
    eps_tree(t+1) = (Y(T0+t+1) - mean(pred_tree(t+1,:)))^2;
end

pred_tree = mean(pred_tree,2);
% =========================================================================
% Compute MSE
mse_tree = mean(eps_tree)

%% PART 1.5: ENSEMBLE

X_ensmbl = [pred_ar pred_pca pred_lasso pred_ridge pred_tree];
num_models = size(X_ensmbl,2);

alpha = inv(X_ensmbl'*X_ensmbl)*X_ensmbl'*Y(T0+1:end);
pred_ensmbl =  X_ensmbl*alpha ;
eps_ensmbl = (Y(T0+1:end) - pred_ensmbl).^2;
mse_ensmbl = mean(eps_ensmbl);


% alpha = zeros(num_models,T-T0);
% eps_ensmbl = zeros(T-T0,1);
% pred_ensmbl = zeros(T-T0,1);
% 
% for t = 1:T-T0
%     
%     alpha(:,t) = inv(X_ensmbl(1:t,:)'* X_ensmbl(1:t,:))* X_ensmbl(t,:)'*Y(T0:T0 +t+-1);
%     pred_ensmbl(t) =  X_ensmbl(t,:)*alpha(:,t) ;
%     eps_ensmbl(t)  = (Y(T0+t) - pred_ensmbl(t))^2;
% end

 
%% PART 1.6: PLOT RESULTS

% =========================================================================
% Plot actual vs prediction
models = {'ar', 'pca', 'lasso', 'ridge', 'tree', 'ensmbl'};
models_titles = {'Autoregression','Factor Model','LASSO','RIDGE','Bagged Tree','Ensemble'};

for i = 1:numel(models)
    
subplot(3,2,i)
hold on
a = plot(datetime(dum.textdata(15:734,1)),Y)
b = plot(datetime(dum.textdata(15:734,1)),[NaN(T0,1);eval(strcat('pred_',string(models(i))))])
hold off
title(string(models_titles(i)), 'interpreter','latex')
xlabel('$t$','interpreter','latex')
ylabel('$\textrm{UNRATE}\; (\upsilon)$','interpreter','latex')
legend([a b], '$\upsilon$', '$\widehat{\upsilon}$', 'interpreter','latex','FontSize',6)

end
saveas(gcf,strcat(figures_dir,'/models_performance'),type)
close

%% PART 1.6: SHOW RESULTS

fileID = fopen(strcat(tables_dir,'\mse_results.txt'),'w');
fprintf(fileID,'Model, MSE\n');
fprintf(fileID,' AR(4),%6.4f\n',mse_ar);
fprintf(fileID,' Factor Model,%6.4f\n',mse_pca);
fprintf(fileID,' LASSO,%6.4f\n',mse_lasso);
fprintf(fileID,' RIDGE,%6.4f\n',mse_ridge);
fprintf(fileID,' Bagged Tree,%6.4f\n',mse_tree);
fprintf(fileID,' Ensemble,%6.4f\n',mse_ensmbl);
fclose(fileID);


%% PART 2: RECESSION/EXPANSION PERFORMANCES
csv_nber = strcat(data_dir,'\USREC.csv');

% Load data from CSV file
nber_data=importdata(csv_nber,',');
recession_dummy = nber_data.data;

% Month/year of final observation
final_datevec_nber=datevec(nber_data.textdata(end,1));
final_month=final_datevec_nber(2);
final_year=final_datevec_nber(1);

% Month/year of first observation
first_datevec_nber=datevec(nber_data.textdata(2,1));
first_month=first_datevec_nber(2);
first_year=first_datevec_nber(1);

% Transform dates
dates_nber = (first_year+first_month:1/12:final_year+final_month/12)';

% Select sample from Jan 1960 to Dec 2019
sample_idx = find(ismember(dates_nber, 1960 + 1/12:1/12:2020));
dates_nber = dates_nber(sample_idx(1):sample_idx(end));
recession_dummy = recession_dummy(sample_idx(1):sample_idx(end),:);
recession_dummy_testset = recession_dummy(T0+1:end);

% =========================================================================
% AR(p)
mse_ar_exp = eps_ar'*(1-recession_dummy_testset)/sum((1-recession_dummy_testset));
mse_ar_rec = eps_ar'*recession_dummy_testset/sum(recession_dummy_testset);

% =========================================================================
% FACTOR MODEL
mse_pca_exp = eps_pca'*(1-recession_dummy_testset)/sum((1-recession_dummy_testset));
mse_pca_rec = eps_pca'*recession_dummy_testset/sum(recession_dummy_testset);

% =========================================================================
% LASSO
mse_lasso_exp = eps_lasso'*(1-recession_dummy_testset)/sum((1-recession_dummy_testset));
mse_lasso_rec = eps_lasso'*recession_dummy_testset/sum(recession_dummy_testset);

% =========================================================================
% RIDGE
mse_ridge_exp = eps_ridge'*(1-recession_dummy_testset)/sum((1-recession_dummy_testset));
mse_ridge_rec = eps_ridge'*recession_dummy_testset/sum(recession_dummy_testset);

% =========================================================================
% BAGGING
mse_tree_exp = eps_tree'*(1-recession_dummy_testset)/sum((1-recession_dummy_testset));
mse_tree_rec = eps_tree'*recession_dummy_testset/sum(recession_dummy_testset);

% =========================================================================
% ENSEMBLE
mse_ensmbl_exp = eps_ensmbl'*(1-recession_dummy_testset)/sum((1-recession_dummy_testset));
mse_ensmbl_rec = eps_ensmbl'*recession_dummy_testset/sum(recession_dummy_testset);

%% PART 2.1: SAVE RESULTS

fileID = fopen(strcat(tables_dir,'\mse_nber_dummies.txt'),'w');
fprintf(fileID,'Model, Expansion, Recession\n');
fprintf(fileID,' AR(4),%6.4f,%6.4f\n',mse_ar_exp, mse_ar_rec);
fprintf(fileID,' Factor Model,%6.4f,%6.4f\n',mse_pca_exp, mse_pca_rec);
fprintf(fileID,' LASSO,%6.4f,%6.4f\n',mse_lasso_exp, mse_lasso_rec);
fprintf(fileID,' RIDGE,%6.4f,%6.4f\n',mse_ridge_exp, mse_ridge_rec);
fprintf(fileID,' Bagged Tree,%6.4f,%6.4f\n',mse_tree_exp, mse_tree_rec);
fprintf(fileID,' Ensemble,%6.4f,%6.4f\n',mse_ensmbl_exp, mse_ensmbl_rec);
fclose(fileID);





