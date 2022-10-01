function Beta = Estimate(Xm,Xf,Ym,Yf,quantiles)
% initialize iterator
gender = ["m" "f"];
Beta = zeros(size(Xm,2),size(quantiles,1),2);
for g = 1:length(gender)
    % Define outer/inner loop structure
    tmp = gender(g);
    X = eval(strcat('X',tmp));
    Y = eval(strcat('Y',tmp));

    %fprintf('Running Gender %s ...\n\n', tmp)

    for t = 1:length(quantiles)

     tau = quantiles(t);
     %fprintf('Running quantile %.3f ...\n',tau)

     options = optimset('MaxIter',100000, 'MaxFunEvals',100000);
     fun = @(x)MLE(X,Y(:,t),x);
     [argmin,fval]  = fminsearch(fun,0.01*rand(size(X,2),1),options);  %exitflag,output

     Beta(:,t,g) = argmin;

    end

end

end