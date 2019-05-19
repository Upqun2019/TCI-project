function Xpre = GMM_3D_JBY(y, CS, para,model)

Mu = model.Mu;
Sig = model.Sig;
pai = model.pai;

[d,Nt] = size(y);
Ri = para.R*eye(d);   %Ri是论文中的R矩阵，是假定的0均值噪声的协方差。


% GMM Inversion
Xpre = zeros(size(CS{1},2),Nt); %用来存用一帧Y恢复得到T帧X的数据。每一行为相应块的数据。

% for i = 1:Nt
parfor i = 1:Nt  % (parallel computing for multicore)
    if mod(i,1000) == 1
        fprintf('i=%d\n',i);
    end
    [pai1,temp,~] = GMM_Inference(y(:,i), CS{i}, Ri, Mu, Sig, pai);     %一块一块恢复。
    Xpre(:,i) = temp;
end


