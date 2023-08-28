function [optlambda, MSE_reg] = PenaltySearch(cv,grid, model)

j = 1;
MSE_reg = zeros(numel(grid),1);
% Grid Search over optimal lambda
for lambda = grid
        
disp(sprintf("Evaluating lambda = %d...\n", lambda))

% Initialize temporary MSE for the folds
MSEcv = 0;
k = size(cv,1);

    % For a given lambda, estimate MSE through block-CV
    for n = 1:k

    % Select training and test set for fold (n)
    trainset = cell2mat(cv([1:n-1 n+1:end]));
    testset = cell2mat(cv(n));
    % Estimate the Lasso for folds (-n)
    
    switch(model)
        case 'lasso'
        B = lasso(trainset(:,2:end),trainset(:,1),'Lambda',lambda);
        
        % Compute the out-of-sample performance given lambda for the chosen fold
        Ttest = size(testset,1);
        eps = zeros(Ttest ,1);

        for h = 0:Ttest-1
            if h == 0
            eps(h+1) = (testset(h+1) - trainset(end,2:end)*B)^2;
            else
            eps(h+1) = (testset(h+1) - testset(h,2:end)*B)^2; 
            end
        end
        case 'ridge'
        B = ridge(trainset(:,1),trainset(:,2:end),lambda);
        
        % Compute the out-of-sample performance given lambda for the chosen fold
        Ttest = size(testset,1);
        eps = zeros(Ttest ,1);

        for h = 0:Ttest-1
            if h == 0
            eps(h+1) = (testset(h+1) - trainset(end,2:end)*B)^2;
            else
            eps(h+1) = (testset(h+1) - testset(h,2:end)*B)^2; 
            end
        end
        
        case 'forest'
        tree = fitrtree(trainset(:,2:end),trainset(:,1));    
        pruned = prune(tree,'Level',lambda); 
        
        % Compute the out-of-sample performance given lambda for the chosen fold
        Ttest = size(testset,1);
        eps = zeros(Ttest ,1);
       
        for h = 0:Ttest-1
            
            if h == 0
            eps(h+1) = (testset(h+1) - predict(pruned,trainset(end,2:end)))^2;
            else
            eps(h+1) = (testset(h+1) - predict(pruned,testset(h,2:end)))^2; 
            end
        end
    end

MSEcv = MSEcv + mean(eps);
end

MSE_reg(j) = mean(MSEcv);

j = j+1;

end

optlambda = grid(find(MSE_reg == min(MSE_reg)));
end