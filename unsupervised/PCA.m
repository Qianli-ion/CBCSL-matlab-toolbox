function [PCs, varianceX, meanX] = PCA(dataMat,thresh)
% principle component analysis on raw dataMat. dataMat is the data matrix
% with each column as a sample and each row as a variable (feature
% dimension). The returned Vorig and firstVariance are the first several
% principle components and corresponding eigenvalues (variance). The
% selection criterion is the thresh (for example thresh=0.9 as 90% variance)

numSample = size(dataMat,2);
numVariable = size(dataMat,1);
% mean normalize
meanX = mean(dataMat,2);
U = bsxfun(@minus, dataMat, meanX);

% if dataMat has less sample than variables
if numSample < numVariable

    sigmaX = U'*U;
    [V,D] = eig(sigmaX);
    [eigValue,I] = sort(diag(D),'descend'); % descend order for eigenvalues
    V = V(:,I);

%     figure
%     plot(eigValue,'r.','MarkerSize',8);

    totalVar = sum(eigValue);
    if thresh <= 1
        curVar = 0;
        for eigIdx = 1:1:length(eigValue)
            curVar = curVar + eigValue(eigIdx);
            if curVar >= thresh*totalVar
                break
            end
        end
    else
        eigIdx = thresh;
    end

    varianceX = eigValue(1:eigIdx);

    % transform back to original space
    PCs = U*V(:,1:eigIdx);
    PCs = bsxfun(@rdivide,PCs,varianceX.^0.5');
else

    sigmaX = U*U';
    [V,D] = eig(sigmaX);
    [eigValue,I] = sort(diag(D),'descend'); % descend order for eigenvalues
    V = V(:,I);

%     figure
%     plot(eigValue,'r.','MarkerSize',8);

    totalVar = sum(eigValue);
    if thresh < 1
        curVar = 0;
        for eigIdx = 1:1:length(eigValue)
            curVar = curVar + eigValue(eigIdx);
            if curVar >= thresh*totalVar
                break
            end
        end
    else
        eigIdx = thresh;
    end
    
    varianceX = eigValue(1:eigIdx);

    % transform back to original space
    PCs = V(:,1:eigIdx);
    
end

end