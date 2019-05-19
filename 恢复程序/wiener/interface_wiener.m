function Xrecon = interface_wiener(Y,C,para)

n = prod(para.patchSize);
delta1 = para.patchSize;
delta2 = para.patchSize;

load('..\..\训练\训练视频\X_train.mat');                                    %可以使用不相关的训练视频集，也可直接使用原始目标视频
Xtst=X_train;

R = R_wiener(Xtst,para);                                                   %利用函数R_wiener，计算自相关矩阵

addpath('..\')
A = Phi2patches_fast(C, para.patchSize(1), para.patchSize(2), para.T, delta1, delta2);  %得到以块处理的测量矩阵

for m=1:para.M
    y = measurement2patches_fast(Y(:,:,m), para.patchSize);                %将每一张观测图像，以块的形式变为一列一列的矩阵。
    y_all(:,:,m)=y;
    
    [d,Nt] = size(y);                                                       %Nt为块的个数
    I=eye(d);
    for i = 1:Nt
        F=A{i};
        W=R*F'*pinv(F*R*F'+I*10^(-8));        
        recon_X(:,i)=W*y(:,i); 
    end
    Xrecon(:,:,(m-1)*para.T+1:m*para.T) = patches2video_fast(recon_X, para.row, para.col, para.patchSize(1), para.patchSize(1), para.T, delta1, delta2);   
end
% save(['..\恢复结果\wiener\block',num2str(patchSize_w(1)),'_T',num2str(para.T),'_F',num2str(para.M),'_wiener'],'Xrecon');