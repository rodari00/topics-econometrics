function obj = GMM(X,Y,beta)


[N,K] = size(X);
A = diag(normpdf(X*beta)./(normcdf(X*beta).*(1-normcdf(X*beta))));
eps = Y - normcdf(X*beta);
g = inv(N)*X'*A*eps;

obj = g'*eye(size(g,1))*g %*N?

end