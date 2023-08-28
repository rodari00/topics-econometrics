%% Compute the forecast of an AR(p)

function forecast = Forecast(dataset,p,h,k,int)
% p # lags
% h horizon
% k last observation i condition on

beta = betaOLS(dataset,p,int);

F = [beta(2:end)'; eye(p-1) zeros(p-1,1)];
c = beta(1);

FF = zeros(p,p,h+1);
FFF = zeros(p,p,h+1);

for i = 1:h+1
    
    FF(:,:,i) = F^(i-1);
    FFF(:,:,i) = F^(i);
    forecast(i) = c*squeeze(sum(FF(1,1,1:i))) + FFF(1,:,i)*flipud(dataset(k-p+1:k));
end 


end