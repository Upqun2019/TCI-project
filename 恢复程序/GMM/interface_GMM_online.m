function Xpre = interface_GMM_online(Y,C,para)

n = prod(para.patchSize);
delta1 = para.patchSize;
delta2 = para.patchSize;
para.R = 1e-6;                                                              %��˹����
para.C =13;                                                                 %�����趨�ٴ�ѵ��ģ���и�˹�˸�����Ҳ�ɺ�ԭʼGMMģ���и�˹�˸�����ͬ

A = Phi2patches_fast(C, para.patchSize(1), para.patchSize(2), para.T, delta1, delta2);

%���س�ʼѵ���ĸ�˹���ģ�Ͳ���
load(strcat(para.filename_training,para.training_model));
model.Mu = Mu;
model.Sig = Sig;
model.pai = pai;
% para.C =size(pai,2);                                                      %��ԭѵ��ģ���и�˹�˵ĸ���

MODEL=cell(1,para.M);
Xt = []; para.kappa = 2;                                                     %kappa��decay rate
for m = 1:para.M
    MODEL{m}=model;
    y = video2patches_fast(Y(:,:,m), para.patchSize(1), para.patchSize(2), 1, delta1, delta2);
    X = GMM_3D_JBY(y, A, para, model);

    Xt = [Xt X];
    Xpre(:,:,(m-1)*para.T+1:m*para.T) = patches2video_fast(X, para.row, para.col, para.patchSize(1), para.patchSize(1), para.T, delta1, delta2);
    
    if m < para.M
        [pai, Mu, Sig] = GMM_online_update_kappa(Xt,para,size(X,2),m);      %�����ع�֡������ѵ��ģ�����ݣ��ٴ�ѵ��ģ��
        model.Mu = Mu;
        model.Sig = Sig;
        model.pai = pai;
    end
end
save([para.filename_training 'model_new.mat'],'MODEL');                     %����ÿһ�θ��º��ģ�Ͳ���







