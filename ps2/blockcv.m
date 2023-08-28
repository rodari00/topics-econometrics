function cv = blocksample(data,k)
% Equal splitting blocked cross-validation, returns a k dimensional cell array
% structure implies the first column is the DV
T = size(data,1);
t = floor(T/k);
cv = cell(k,1);

for b = 0:k-1

cv{b+1}= data(b*t + 1:(b+1)*t,:);

if b == k-1
 cv{b+1} = data(b*t +1:(b+1)*t + rem(T,t),:);
end
    
end

end
