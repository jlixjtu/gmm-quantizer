function [ids] = GMVQ_nn_search(query, codewordg, lij, indexc, indexq, model, T)
% input
% query: the query vectors
% codewordg: universal look-up table
% lij: the allocated intervals 
% indexc: the index for each cluster
% indexq: the code
% T: the number of returned sample
% output
% ids: the returned sample

Nb = size(indexq,1);
Nq = size(query,2);
[p, m] = size(model.mu);

Q = cell(1,m);% Q matrix of cluster
slamda = zeros(p,m);% geometry means for lamdas of cluster
% bit allocation stage1: bit allocation among clusters
for ci=1:m
    covmxi=model.sigma(:,:,ci);
    %[v d u]=svd(covmxi);% AV=VD decomposition
    [v d]=eig(covmxi);% AV=VD decomposition d: already sorted
    Q{ci}=v;
    lamda(:,ci)=diag(d);
    slamda(:,ci)=sqrt(lamda(:,ci));
end

for ii = 1:Nq
%    c = queryc(ii);
%     cand = find(c == indexc); % the data which has same cluster is defined as candidant
%     table = codewordg(quantizer(c,:),:);
    C_table = cell(1,m);
    for ci=1:m % m clusters
        b = zeros(p,1);
        table = zeros(p,size(codewordg,2));
        tmp = query(:,ii) - model.mu(:,ci);
        tmp1 = Q{ci}'*tmp;
        id = find(lij(ci,:)~=1, 1 );
        distortionc(ci) = sum(tmp1(1:id-1).^2);
        for bidx = id:p
            ll = lij(ci,bidx);                        
            C = codewordg(ll,1:ll);
            distmp = ((tmp1(bidx)-C.*slamda(bidx,ci))).^2;
            distortionc(ci) = distortionc(ci) + min(distmp);

            table(bidx,1:ll) = distmp;
        end
        C_table{ci} = table;
    end
    [~, minc] =  sort(distortionc,'ascend');
    
    clear table
    clear cand
    clear dist
    for jj = 1:ceil(m/3.99)
    ci = minc(jj);
    table{jj} = C_table{ci};
    cand{jj} = find(indexc(:,1)==ci);
    id = find(lij(ci,:)~=1, 1 );
    dist{jj} = zeros(length(cand{jj}),1);
    tmp = query(:,ii) - model.mu(:,ci);
    tmp1 = Q{ci}'*tmp;
    dist{jj} = dist{jj} + sum(tmp1(1:id-1).^2);
        for bidx = id:p
            dist{jj} = dist{jj} + table{jj}(bidx,indexq(cand{jj},bidx))';
        end
    end
    dist = single(cat(1, dist{:}));
    cand = single(cat(1, cand{:}));
    [~,idx] = sort(dist,'ascend');
    
    ids(ii,:) = cand(idx(1:T));
    %dis(ii,:) = dummy(1:T);
    if ~mod(ii,20)
        fprintf('.');
    end
    if ~mod(ii,1000)
        fprintf('%d data is processed \n',ii);
    end
end
