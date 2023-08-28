%%% Create function providing OLS estimates for an AR(p) process
function [beta, eps, Y] = betaOLS(dataset,p,int)

% Define dependent and independent variables
Y = dataset(p+1:end);

switch(int)
    case 'yes'
    X =[ones(size(dataset,1)-p,1)];
    case 'no'
    X = [];
    otherwise
    fprintf('Provide a valid input for intercept!\n' );  
end

for i = 1:p
 X = [X dataset(p-i+1:end-i)];
end 

% Compute the vector of parameters
beta = inv(X'*X)*X'*Y;
% Get errors
eps = Y - X*beta;


end 