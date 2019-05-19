function [pai, Mu, Sig] = GMM_online_update_kappa(Xt,para,Nt,m)

% Update GMM prior  下载训练模型的原始数据
load(strcat(para.filename_training,para.training_data0));
nd = size(data0,2);    

tt=size(Xt,2);

p = randperm(nd); 
data1 = data0(:,p(1:nd-tt));   %作用：随机选取data0的一部分列，剔除列数和前m次恢复次数有关，和下降系数有关。

data1 = [data1 Xt];

data1=data1(:,p);

clear data0 Xpre

options = statset('Display','final','MaxIter',200);
obj = gmdistribution.fit(data1',para.C,'Regularize',10^-3,'Options',options);
clear data1
pai = obj.PComponents;
Mu = obj.mu';
Sig = obj.Sigma;