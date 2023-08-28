function B = blockbootstrp(data,blockSize,nSamples)

dataT = data';
numBlocks = size(dataT,2) / blockSize;           % must be integer
blocks = reshape(dataT, [size(dataT,1)*blockSize,numBlocks])';
samples = bootstrp(nSamples, @(x)x', blocks)';

% Initialize Bootstrap
B = zeros(size(data,1),size(data,2),nSamples);

for b = 1:nSamples
B(:,:,b) = reshape(samples(:,b),[size(data,2),size(data,1)])';

end


end
