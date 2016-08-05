% this codes only include .m file that run directly
% test dataset is siftsmall which is a part of sift 1M
% author: 
%   Jin Li    
% email:
% lijin199102@gmail.com

clear
dataset = 'siftsmall';
pq_test_load_vectors;

%% off-line training Gaussian kernel
for nGMM =[1 4 16];    % the number of Gaussian kernels
 fprintf('================== %d cluster GMM =========================\n',nGMM)   
modelfile = [dataset '_gmm_model_' int2str(nGMM) '.mat'];
if exist(modelfile,'file') 
    load(modelfile);
else
    fprintf('parameter is not exist\n');
    fprintf('======================training Gaussian kernel=======================\n')
    tic;
    options = statset('Display','final','MaxIter',100);
    obj = gmdistribution.fit(vtrain',nGMM,'Regularize',10^-3,'Options',options);%GMM TOOLS
    pai = obj.PComponents;
    Mu = obj.mu';
    Sig = obj.Sigma;
    model.weight = pai;
    model.sigma =Sig ;
    model.mu = Mu;
    save(modelfile,'pai','Mu','Sig','model');
    t1 = toc;
    fprintf('%d Gaussian kernel has been trained\n', nGMM)
    fprintf('training time is %.2f\n',t1);
end


for nbit = [16 32 64 128];     % bit length of each code
%% GMPQ encoding for arbitrary given code length
fprintf('=================encoding vectors into %d bits ===================\n', nbit)
tic;
[distortion indexc indexq codewordg lij] = GMVQ(vbase,nbit,model);
t2 = toc;
fprintf('encoding time is: %.2f\n',t2);
fprintf('distortion is: %.2f\n ', distortion)

%% GMPQ searching 
fprintf('=============searching nearest neighbors for query ==============\n')
tic;
[ids] = GMVQ_nn_search(vquery(:,1:100), codewordg, lij, indexc, indexq, model, 100);
t3 = toc;
fprintf('searching time is: %.2f\n',t3);

T = [ 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000];
for t = 1:7
    acc(t) = 0;
    for ii = 1:size(ids,1)
        if ~isempty(find(ids(ii,1:T(t))==ids_gnd(ii), 1))
            acc(t) = acc(t)+1;
        end
    end
    
end
acc
end
end

 
