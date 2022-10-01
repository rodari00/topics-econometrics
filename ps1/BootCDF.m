function [mm,ww,mw] =  BootCDF(Data,tau)
quantiles = 2.5:2.5:97.5;
Xnames = {'female',...
           'divorced', 'separated', 'nevermarried', 'married',...
           'hsd911', 'hsg', 'sc', 'cg', 'ad'...
           'so',  'we',  'ne',...
           'exp1', 'exp2', 'exp3', 'exp4'};   
       
X = table2array(Data(:,Xnames));
Y = Data.lnw;

% Female
Xf = X(X(:,1)==1,2:end);
Xf = [ones(size(Xf,1),1) Xf];
% Each column for each gender is the indicator function for the CDF at the
% given quantile
Yf = Y(Data.female==1) < tau;

% Male
Xm =  X(X(:,1)==0,2:end);
Xm = [ones(size(Xm,1),1) Xm];
% Each column for each gender is the indicator function for the CDF at the
% given quantile
Ym = Y(Data.female==0) < tau; %prctile(Y,quantiles)

%% Estimation
Beta = Estimate(Xm,Xf,Ym,Yf,tau);
    
% Report the CDFs
mm = mean(normcdf(Xm*Beta(:,:,1)),1);
ww = mean(normcdf(Xf*Beta(:,:,2)),1);
mw = mean(normcdf(Xf*Beta(:,:,1)),1);
%tau = prctile(Y,quantiles); add it to the function?

end
