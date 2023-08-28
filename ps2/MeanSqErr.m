function [mse f se] = MeanSqErr(data, p, h, k,int, m)
% k is the #th observation that splits the sample in 2
% m is the method applies (1 = recursive, 0 = rolling window)
% h is the forecast horizon
% p is the lag of the AR(p) process
% data must be a row vector Tx1

%Define size of the sample
T = size(data,1);
se = zeros(T-k,1);
f = zeros(T-k,1);
if m == 1
    disp 'Applying recursive method...'
   
    for i = k+1:T-h
        datafc = data(1:i-1);
        f(i-k) = Forecast(datafc,p,h-1,i-1,int);
        se(i-k) = (data(i-1+h) - f(end))^2;
    end 
        mse = mean(se);
        disp (['The MSE obtained through the recursive method at horizon ', num2str(h), ' is ', num2str(mse)]) 
elseif m == 0 
       disp 'Applying rolling method...'
    for i = k+1:T-h
        datafc = data(i-k:i-1);
        kk = size(datafc,1);
        f(i-k) = Forecast(datafc,p,h-1,kk,int);
        se(i-k) = (data(i-1+h) - f(end))^2;
    end 
        mse = mean(se);
        disp (['The MSE obtained through the rolling method at horizon ', num2str(h), ' is ', num2str(mse)])
else 
    disp 'An error has occurred. Input m is expected to be either 1 or 0.'

end 

end