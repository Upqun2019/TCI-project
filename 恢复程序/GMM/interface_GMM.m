function Xpre = interface_GMM(Y,C,para)

n = prod(para.patchSize);
delta1 = para.patchSize;
delta2 = para.patchSize;
para.R = 1e-6;%��˹����

A = Phi2patches_fast(C, para.patchSize(1), para.patchSize(2), para.T, delta1, delta2);  %�õ��Կ鴦��Ĳ�������

%���س�ʼѵ���ĸ�˹���ģ�Ͳ���
load(strcat(para.filename_training,para.training_model));
model.Mu = Mu;
model.Sig = Sig;
model.pai = pai;
para.C =size(pai,2);   %ѵ��ģ���к˵ĸ���

for m = 1:para.M
    y = video2patches_fast(Y(:,:,m), para.patchSize(1), para.patchSize(2), 1, delta1, delta2);    %���۲⵽����Ƶ��һ֡һ֡�ָ���
    %%��Y(:,:,m)������һ֡���ֿ�������������һ����ά����
    X = GMM_3D_JBY(y, A, para, model);
    Xpre(:,:,(m-1)*para.T+1:m*para.T) = patches2video_fast(X, para.row, para.col, para.patchSize(1), para.patchSize(1), para.T, delta1, delta2); 
    %%����ԭ�����Ƶ�洢��Xpre�����С�1֡�ָ���T֡��
end







