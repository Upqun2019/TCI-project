function Xpre = GMM_3D_JBY(y, CS, para,model)

Mu = model.Mu;
Sig = model.Sig;
pai = model.pai;

[d,Nt] = size(y);
Ri = para.R*eye(d);   %Ri�������е�R�����Ǽٶ���0��ֵ������Э���


% GMM Inversion
Xpre = zeros(size(CS{1},2),Nt); %��������һ֡Y�ָ��õ�T֡X�����ݡ�ÿһ��Ϊ��Ӧ������ݡ�

% for i = 1:Nt
parfor i = 1:Nt  % (parallel computing for multicore)
    if mod(i,1000) == 1
        fprintf('i=%d\n',i);
    end
    [pai1,temp,~] = GMM_Inference(y(:,i), CS{i}, Ri, Mu, Sig, pai);     %һ��һ��ָ���
    Xpre(:,i) = temp;
end


