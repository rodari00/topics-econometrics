%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% ECON8825 - Problem Set 1 %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup folders
[figures_dir,data_dir,tables_dir] = directorysetup(cd,'Win');
type = 'png';
B = 100; % should be enough
% Estimates and boostrapped estimates are already included in .mat files
%% Setup and Data Import
%T = readtable(strcat(data_dir,'\cps2012.csv'));
T = readtable(strcat(data_dir,'/cps2012.csv'));
T.married  = strcmp(T.married,'TRUE');
T.ne  = strcmp(T.ne,'TRUE');
T.sc  = strcmp(T.sc,'TRUE');


%opts = detectImportOptions(strcat(data_dir,'\cps2012.csv')); % Windows
opts = detectImportOptions(strcat(data_dir,'/cps2012.csv'));
opts = setvartype(opts,{'married','ne','sc'},'double');
VarNames = T.Properties.VariableNames;
Xnames = {'female',...
           'divorced', 'separated', 'nevermarried', 'married',...
           'hsd911', 'hsg', 'sc', 'cg', 'ad'...
           'so',  'we',  'ne',...
           'exp1', 'exp2', 'exp3', 'exp4'};   
VarTypes = opts.VariableTypes;
NumericVars = find(strcmpi(VarTypes,'double'));
NumericVars = NumericVars(NumericVars ~= find(strcmpi(VarNames,'year')));


 
%% Part 1) Summary statistics by gender
expression = '(^|\s).';
replace = '${upper($0)}';

stats_gender = grpstats(T(:,NumericVars),...
                        'female',...
                        {'mean','std','min','max','skewness'});
         
% Adjustments before exporting to LaTeX
str = regexprep(...
       regexprep(...
       stats_gender.Properties.VariableNames(2:end),'_',' '),...
                 expression,replace);

stats = cell2table([str',...
                    num2cell(stats_gender{:,2:end}')],...
                    'VariableNames',{'Summary','Male','Female'});
                
writetable(stats,strcat(tables_dir,'\summary-stats.csv'));                
                
table2latex(stats,strcat(tables_dir,'\summary-stats.tex'))

%% Part 3a) Create Data

% Note: choose the quantiles over the whole distributuion to be able to
% compare men and women (F(y|m) - F(y|f) have the same y for each quantiles
%, to avoid different y for the same quantiles across genders)

% Vector of quantile cutoffs
quantiles = 2.5:2.5:97.5;

% Define X and Y matrices
X = table2array(T(:,Xnames));
Y = T.lnw;
tau = prctile(Y,quantiles);

% Female
Xf = X(X(:,1)==1,2:end);
Xf = [ones(size(Xf,1),1) Xf]; % remove gender dummy and add constant column
% Each column for each gender is the indicator function for the CDF at the
% given quantile
Yf = Y(T.female==1) < tau; 

% Male
Xm =  X(X(:,1)==0,2:end);
Xm = [ones(size(Xm,1),1) Xm]; % remove gender dummy and add constant column
% Each column for each gender is the indicator function for the CDF at the
% given quantile
Ym = Y(T.female==0) < tau;

% Empirical CDF
hold on
a = plot(tau,mean(Ym),'-k');
b = plot(tau,mean(Yf),'-.k');
hold off
xlabel('log(w)','interpreter','latex')
ylabel('$F(\cdot)$','interpreter','latex')
legend([a,b], 'Male','Female', 'interpreter','latex')
saveas(gcf,strcat(figures_dir,'/empirical_cdf'),'pdf')


%% Part 3b) Estimate by gender and quantile
%  Beta is a KxTxG where G is ordered as (M,F)

% Check if Beta is already loaded to avoid running it again
if isfile(strcat(data_dir,'\beta.mat')) 
 
    load(strcat(data_dir,'\beta.mat'))

elseif isfile(strcat(data_dir,'/beta.mat'))
    
    load(strcat(data_dir,'/beta.mat'))
    
else
    
   
    % initialize iterators
    gender = ["m" "f"];
    Beta = zeros(size(X,2),size(quantiles,1),2);

    tic
    for g = 1:length(gender)
        % Define outer/inner loop structure
        tmp = gender(g);
        X = eval(strcat('X',tmp));
        Y = eval(strcat('Y',tmp));

        fprintf('Running Gender %s ...\n\n', tmp)

        for t = 1:length(quantiles)

         tmp_tau = quantiles(t);
         fprintf('Running quantile %.3f ...\n',tmp_tau)

         options = optimset('MaxIter',10000, 'MaxFunEvals',10000);
         fun = @(x)MLE(X,Y(:,t),x);
         [argmin,fval]  = fminsearch(fun,0.01*rand(size(X,2),1),options);  %exitflag,output

         Beta(:,t,g) = argmin;

        end

    end
    toc
    % Save Beta
    save(strcat(data_dir,'/beta'),'Beta')

end


% Plot distribution of betas
% for k = 2:17
% subplot(4,4,k-1)
% hold on
% a =histogram(Beta(k,:,1),20);
% b = histogram(Beta(k,:,2),20);
% hold off
% xlabel('$\beta$','interpreter','latex')
% ylabel('N','interpreter','latex')
% title(string(Xnames(k)),'Interpreter','latex');
% legend([a,b], 'female','male','interpreter','latex','FontSize',5)
% end
% 
% saveas(gcf,strcat(figures_dir,'/beta_distr'),type)


% Plot distribution of betas
for k = 2:17
subplot(4,4,k-1)
hold on
histogram(Beta(k,:,1) -Beta(k,:,2) ,20);
hold off
xlabel('$\beta$','interpreter','latex')
ylabel('N','interpreter','latex')
title(string(Xnames(k)),'Interpreter','latex');
legend('$\beta_m - \beta_f$','interpreter','latex')
end

saveas(gcf,strcat(figures_dir,'/dbeta_distr'),type)


% Plot the counterfactual
hold on
a = plot(tau,mean(normcdf(Xf*Beta(:,:,2)),1),'-k'); %ww
b= plot(tau,mean(normcdf(Xm*Beta(:,:,1)),1),'-.k'); %mm
c = plot(tau,mean(normcdf(Xf*Beta(:,:,1)),1),'-*k'); %mw
hold off
xlabel('log(w)','interpreter','latex')
ylabel('$F(\cdot)$','interpreter','latex')
legend([a,b,c], 'female','male', 'male $\mid$female','interpreter','latex',)

saveas(gcf,strcat(figures_dir,'/counterfactual'),type)

%% Part 4) Decomposition

% Estimated CDFs
mm = mean(normcdf(Xm*Beta(:,:,1)),1);
ww = mean(normcdf(Xf*Beta(:,:,2)),1);
mw = mean(normcdf(Xf*Beta(:,:,1)),1);

% Plot Price and Composition Effect
hold on
a = plot(tau,...
         mm- mw,'-.k');

b = plot(tau,...
         mw - ww,'-ok');
hold off
xlabel('log(w)','interpreter','latex')
ylabel('$Size$','interpreter','latex')
legend([a,b], 'Composition Effect', 'Price Effect','interpreter','latex')

saveas(gcf,strcat(figures_dir,'/price_comp_effect'),type)

%% Part 5) Boostrap and Plot Confidence Bands


% Check if Bootstrap samples are already loaded to avoid running it again
if isfile(strcat(data_dir,'\MM.mat')) 
 
     load(strcat(data_dir,'/MM'));
     load(strcat(data_dir,'/WW'));
     load(strcat(data_dir,'/MW'));

elseif isfile(strcat(data_dir,'/MM.mat')) 
    
     load(strcat(data_dir,'/MM'));
     load(strcat(data_dir,'/WW'));
     load(strcat(data_dir,'/MW'));
    
else
    
    MM= zeros(B,size(quantiles,2));
    WW = zeros(B,size(quantiles,2));
    MW = zeros(B,size(quantiles,2));
    %Tau = zeros(B,size(quantiles,2));


    parfor b = 1:B
        fprintf('B = %d ...\n',b)
        Boot = datasample(T,size(T,1),1);
        [MM(b,:),WW(b,:),MW(b,:)] =  BootCDF(Boot,tau);

    end


     % Save Beta
     save(strcat(data_dir,'/MM'),'MM');
     save(strcat(data_dir,'/WW'),'WW');
     save(strcat(data_dir,'/MW'),'MW');

end


% Compute variance
S = [sum((MM-mm).^2)./B;...
     sum((WW-ww).^2)./B;...
     sum((MW-mw).^2)./B];

% Compute statistic
mb = [max(abs((MM - mm)./S(1,:)),[],2)';...
      max(abs((WW - ww)./S(2,:)),[],2)';...
      max(abs((MW - mw)./S(3,:)),[],2)'];

% Compute critical values
alpha = 0.05;
ci = prctile(mb,1-alpha,2);


specs = ["mm","ww","mw"];
names = ["Male" "Female" "Male | Female"];

% Plot Bootstrapped CDFs
for i = 1:length(specs)

tmp = eval(specs(i));
tmptitle = names(i);
subplot(3,1,i)

% Plot Confidence bands
hold on
p1 = plot(tau, tmp - ci(i)*S(i,:),'-.k') ;
p2 = plot(tau, tmp,'-r');
p3 = plot(tau, tmp + ci(i)*S(i,:),'-.k');
hold off
xlabel('log(w)','interpreter','latex')
ylabel('$F(\cdot)$','interpreter','latex')
legend([p1,p2,p3], '$\hat{F}(y) - c_{1-\alpha}\hat{s}(y)$',...
                 '$\hat{F}(y) $',...
                   '$\hat{F}(y) + c_{1-\alpha}\hat{s}(y)$', 'interpreter','latex')

title(sprintf('F(y) of %s', tmptitle), 'Interpreter','latex')

end

saveas(gcf,strcat(figures_dir,'/bootstrap_cdf'),type)


%% Part 6) Confidennce Bands for Price and Composition Effect

cband_l = zeros(3,size(tau,2));
cband_h = zeros(3,size(tau,2));


% Store confidence bands
for i = 1:length(specs)

tmp = eval(specs(i));

cband_l(i,:) = tmp - ci(i)*S(i,:);
cband_h(i,:) = tmp + ci(i)*S(i,:);

end

% Plot Price and Composition Effect
hold on
% CE
plot(tau,...
     cband_l(1,:) - cband_h(3,:),'-.k');
ce = plot(tau,...
         mm- mw,'-r');
plot(tau,...
     cband_h(1,:) - cband_l(3,:),'-.k');

% PE
plot(tau,...
     cband_l(3,:) - cband_h(2,:),'-.k');
pe = plot(tau,...
         mw - ww,'-or');
plot(tau,...
     cband_h(3,:) - cband_l(2,:),'-.k');
hold off
xlabel('log(w)','interpreter','latex')
ylabel('$Size$','interpreter','latex')
legend([pe,ce], 'Price Effect',...
                 '$Composition Effect$',...
                   'interpreter','latex')

saveas(gcf,strcat(figures_dir,'/price_comp_cband'),type)
