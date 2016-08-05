function [distortion indexc indexq codewordg lij] = GMVQ(datay, bits,  model)
% input
% datay is p*N matrix, p is dimension of data, N is number of data
% bits is the total bits for each vector
% model:
%        model.sigma: covariance
%        model.mu: mean vector
%        model.weight: alpha in GMM

% output
% distortion: the mean quantization distortion
% indexc: the index of cluster 
% indexq: the code
% codewordg: the universal codewords look-up table
% lij: intervals for each conponents in each cluster
% authors: Jin Li & Xiangwei Li

%% bit allocation with PDF
nt=size(model.mu,1);% signal dimension
m=size(model.weight,2);% components number of GMM
btot=bits;
p=nt;% p-dimensional Gaussian case
lamda=zeros(nt,m);% lamdas of cluster
Q=cell(1,m);% Q matrix of cluster
clamda=[];% geometry means for lamdas of cluster
saxc=0;% for bit allocation stage2
slamda=zeros(nt,m);%sqrted lamda;
% bit allocation stage1: bit allocation among clusters
for ci=1:m
    covmxi=model.sigma(:,:,ci);
    %[v d u]=svd(covmxi);% AV=VD decomposition
    [v d]=eig(covmxi);% AV=VD decomposition d: already sorted
    Q{ci}=v;
    lamda(:,ci)=diag(d);
    slamda(:,ci)=sqrt(lamda(:,ci));
    clamda(ci)=prod(abs(lamda(p:-1:1,ci)).^(1/p));% (lamda1*lamda2...lamdan)^(1/p)
    saxc = saxc+(model.weight(ci)*clamda(ci))^(p/(p+2));
end

bc=zeros(m,1);% bits used to quantize clusters
bij=zeros(m,p);
%bit allocation stage2: bit allocation within a cluster
for ci = 1:m
    bc(ci) = btot+(p/(p+2))*log2(model.weight(ci)*clamda(ci))-log2(saxc);
    pp = p;
    for bidx = 1:p
        bij(ci,bidx) = bc(ci)/pp+0.5*log2(lamda(bidx,ci)/clamda(ci));
        if bij(ci,bidx) <0;
            bij(ci,bidx) = 0;
            pp = pp-1;
            clamda(ci) = prod(abs(lamda(bidx+1:p,ci)).^(1/pp));
        end
    end%bidx
end%ci

% make nonnegative integer number of levels
lij=ceil(2.^bij);
ppi=zeros(m,1);
nik=zeros(m,p);

for ci=1:m
    for k=1:p;
        ppi(ci)=prod(lij(ci,:));
        tmp=1-2^bc(ci)/ppi(ci);
        nik(ci,k)=floor(lij(ci,k)*tmp);
        lij(ci,k)=lij(ci,k)-nik(ci,k);
        if(lij(ci,k)<1)
            lij(ci,k)=1;
        end
    end%k
end%ci


%% codewords generating
lmax=max(max(lij));
ls=unique(lij);
codewordg=zeros(lmax,lmax);
%partitiong=codewordg;
% 1.codebook design
% optimal quantizer design for gaussian distribution
pg=img_csv_gaupdfexp3(0,-3.4,3.4,4.5,1);% ymax=-ymin=5;optimal
% quantizers for the gaussian ensemble with the same Y/epsilon;
clamdaig=pg.cdf3;
xpg=pg.x;
pdf = pg.pdf3;
ls = reshape(ls,1,length(ls));
for ll=ls
  codewordt{1}=[];
  for i=1:ll
    cdfv=(i-1/2)/ll;
    codewordt{1}=[codewordt{1};xpg(find(clamdaig>=cdfv,1,'first'))];
  end
  codewordg(ll,1:ll)=codewordt{1};
end

third_term = cell(1,m);
for ci = 1:m
    third_term{ci} = zeros(p,size(codewordg,1));
    for bidx = 1:p
        ll = lij(ci,bidx);
        C = codewordg(ll,1:ll);
        third_term{ci}(bidx,1:ll) = (slamda(bidx,ci)*C).^2;
    end
end

%% encoding phase
N = size(datay,2);
batchsize = 10000;
batchnum = N/batchsize;
indexq = uint8(zeros(N,p));
%indexc = uint8(zeros(N,1));

for batch = 1:batchnum
indextmp = uint8(zeros(batchsize,p));
distortionc = zeros(batchsize,m);
databatch = datay(:,(batch-1)*batchsize+1:batch*batchsize);
% select a cluster with minimal distortion
for ci = 1:m % m clusters
    tmp = bsxfun(@minus, databatch , model.mu(:,ci)); % remove mean vector
    id = find(lij(ci,:)~=1, 1 );
    distortionc(:,ci) = sum(tmp.^2);
    for bidx = id:p
        ll = lij(ci,bidx);
        C = codewordg(ll,1:ll);
        distmat = bsxfun(@plus, -2*slamda(bidx,ci)*(tmp'*Q{ci}(:,bidx))*C, third_term{ci}(bidx,1:ll));
        distortionc(:,ci) = distortionc(:,ci) + min(distmat,[],2);
    end
end
[dist((batch-1)*batchsize+1:batch*batchsize) minc] = min(distortionc,[],2);
% encoding with the selected cluster
for ci = 1:m
    cand = find(minc(:,1)==ci);
    tmp = bsxfun(@minus, databatch(:,cand) , model.mu(:,ci));
    id = find(lij(ci,:)~=1, 1 );
    for bidx = id:p
        ll = lij(ci,bidx);
        C = codewordg(ll,1:ll);
        distmat = bsxfun(@plus, -2*slamda(bidx,ci)*(tmp'*Q{ci}(:,bidx))*C, third_term{ci}(bidx,1:ll));
        [~, mni] = min(distmat,[],2);
        indextmp(cand,bidx) = uint8(mni);
    end
end
indexc((batch-1)*batchsize+1:batch*batchsize,:) = minc;
indexq((batch-1)*batchsize+1:batch*batchsize,:) = indextmp;
end
 
distortion = mean(dist);
