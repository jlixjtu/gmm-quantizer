function [ids] = nn_search(query, base, T)
Nq = size(query,2);
for ii = 1:Nq
    dist = sum(bsxfun(@minus, query(:,ii), base).^2);
    [dummy,idx] = sort(dist,'ascend'); 
    ids(ii,:) = idx(1:T);
    if ~mod(ii,20)
        fprintf('.');
    end
    if ~mod(ii,1000)
        fprintf('%d data is processed \n',ii);
    end
end