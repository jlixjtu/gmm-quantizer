% gnd: n*N matrix: n is the number of nearest neighbors;
%                    N is the number of query vectors
% ids: m*N matrix: m retrieval terms
% K: ap@K K = m for default

function mAP = mAPsearch(gnd, ids, K)

n = size(gnd, 1);  % n nearest neighbors in ground truth
[m N] = size(ids);  % N query; m returned retrieval idx
if nargin == 3
    m = K;
end
AP = zeros(N, 1);
for ii = 1:N
    num = 0;
    for jj = 1:n
        if ~isempty(find(gnd(jj,ii) == ids(1:m,ii),1))
            num = num +1;
            AP(ii) = AP(ii) + num/jj;
        end
    end
    AP(ii) = AP(ii)/m;
end
mAP = mean(AP);