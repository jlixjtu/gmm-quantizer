% VQ_train.m
%
% get GMM model for VQ from training data
% 
% Written by: xwli
% Created: Nov. 2014; May. 2015
%
clear all
%
% % include path
% % addpath ('D:\Xiangwei\cs recovery model\Utils')
% % addpath ('D:\Xiangwei\cs recovery model\Data\timg')
dataset = 'sift';
pq_test_load_vectors;
for nGMM = [16 32 64 128];
tic;
options = statset('Display','final','MaxIter',100);
obj = gmdistribution.fit(vtrain',nGMM,'Regularize',10^-3,'Options',options);%GMM TOOLS
pai = obj.PComponents;
Mu = obj.mu';
Sig = obj.Sigma;
model.weight = pai;
model.sigma =Sig ;
model.mu = Mu;
save([dataset '_gmm_model_' int2str(nGMM)],'pai','Mu','Sig','model');
toc
end
