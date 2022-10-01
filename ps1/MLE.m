function obj = MLE(X,Y,beta)

obj =-sum(Y.*log(normcdf(X*beta)) + (1-Y).*log(1-normcdf(X*beta)));
end