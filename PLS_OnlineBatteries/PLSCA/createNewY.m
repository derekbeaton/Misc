% function [meanBootSample,samples] = bootstrapBlockFS(Y,F_sup_blocks,massedY,samples)
%     if ~exist('samples')
%         samples=1000;
%     end
% 
%     [I,fJ,N] = size(F_sup_blocks);
%     [K] = size(Y,2);
%     
%     meanBootSample = zeros(K,fJ,samples);
%     for i=1:samples
%         subjectBoot = randsample(N,N,'true');
%         rowCatBoot = createNewY(Y);
%         meanBootSample(:,:,i) = massedY'*mean(F_sup_blocks(rowCatBoot,:,subjectBoot),3);
%     end
%         
% end

%This should live here, for now. It might get moved later.
function [rowBootIndices] = createNewY(Y)
    [I,K] = size(Y);
    rowBootIndices = zeros(I,1);
    for i=1:K
        theseIndices = find(Y(:,i)==1);
        rowBootIndices(min(theseIndices):max(theseIndices))=randsample((min(theseIndices):max(theseIndices)),size(theseIndices,1),'true');
    end
end