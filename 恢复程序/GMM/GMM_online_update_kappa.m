function [pai, Mu, Sig] = GMM_online_update_kappa(Xt,para,Nt,m)

% Update GMM prior  下载训练模型的原始数据
load(strcat(para.filename_training,para.training_data0));
nd = size(data0,2);    

for i = 0:m
    xi(i+1) = (i+1)^para.kappa;
end

xi = xi./sum(xi);
ndt = Nt/xi(end);


nm = round(ndt*xi); 
p = randperm(nd); 
data1 = data0(:,p(1:nd-sum(nm(2:end))));   %作用：随机选取data0的一部分列，剔除列数和前m次恢复次数有关，和下降系数有关。

for i = 1:m
    p = randperm(Nt);
    temp = Xt(:,(i-1)*Nt+1:i*Nt);
    if nm(i+1)<Nt
        data1 = [data1 temp(:,p(1:nm(i+1)))];  %将上述data1和前（m-1）次恢复得到的X的部分列结合
    else
        data1 = [data1 temp];                  %最后将最后一次（第m次）全部列随意排列，和data1结合。
    end
end
%%data1结合了恢复的数据X和data0。

clear data0 Xpre

options = statset('Display','final','MaxIter',200);
obj = gmdistribution.fit(data1',para.C,'Regularize',10^-3,'Options',options);
clear data1
pai = obj.PComponents;
Mu = obj.mu';
Sig = obj.Sigma;