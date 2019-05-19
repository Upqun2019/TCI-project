function Xpre = interface_GMM(Y,C,para)

n = prod(para.patchSize);
delta1 = para.patchSize;
delta2 = para.patchSize;
para.R = 1e-6;%高斯噪声

A = Phi2patches_fast(C, para.patchSize(1), para.patchSize(2), para.T, delta1, delta2);  %得到以块处理的测量矩阵

%下载初始训练的高斯混合模型参数
load(strcat(para.filename_training,para.training_model));
model.Mu = Mu;
model.Sig = Sig;
model.pai = pai;
para.C =size(pai,2);   %训练模型中核的个数

for m = 1:para.M
    y = video2patches_fast(Y(:,:,m), para.patchSize(1), para.patchSize(2), 1, delta1, delta2);    %将观测到的视频，一帧一帧恢复。
    %%将Y(:,:,m)即其中一帧，分块成列向量，组成一个二维矩阵。
    X = GMM_3D_JBY(y, A, para, model);
    Xpre(:,:,(m-1)*para.T+1:m*para.T) = patches2video_fast(X, para.row, para.col, para.patchSize(1), para.patchSize(1), para.T, delta1, delta2); 
    %%将复原后的视频存储到Xpre矩阵中。1帧恢复成T帧。
end







