function Xrecon = interface_wiener(Y,C,para)

n = prod(para.patchSize);
delta1 = para.patchSize;
delta2 = para.patchSize;

load('..\..\ѵ��\ѵ����Ƶ\X_train.mat');                                    %����ʹ�ò���ص�ѵ����Ƶ����Ҳ��ֱ��ʹ��ԭʼĿ����Ƶ
Xtst=X_train;

R = R_wiener(Xtst,para);                                                   %���ú���R_wiener����������ؾ���

addpath('..\')
A = Phi2patches_fast(C, para.patchSize(1), para.patchSize(2), para.T, delta1, delta2);  %�õ��Կ鴦��Ĳ�������

for m=1:para.M
    y = measurement2patches_fast(Y(:,:,m), para.patchSize);                %��ÿһ�Ź۲�ͼ���Կ����ʽ��Ϊһ��һ�еľ���
    y_all(:,:,m)=y;
    
    [d,Nt] = size(y);                                                       %NtΪ��ĸ���
    I=eye(d);
    for i = 1:Nt
        F=A{i};
        W=R*F'*pinv(F*R*F'+I*10^(-8));        
        recon_X(:,i)=W*y(:,i); 
    end
    Xrecon(:,:,(m-1)*para.T+1:m*para.T) = patches2video_fast(recon_X, para.row, para.col, para.patchSize(1), para.patchSize(1), para.T, delta1, delta2);   
end
% save(['..\�ָ����\wiener\block',num2str(patchSize_w(1)),'_T',num2str(para.T),'_F',num2str(para.M),'_wiener'],'Xrecon');