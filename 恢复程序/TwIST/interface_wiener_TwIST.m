%% 先用维纳算法小块恢复，然后将其作为初值，TwIST大块恢复
function Xpre = interface_wiener_TwIST(Y,C,para,orig_str)

n = prod(para.patchSize);
delta1 = para.patchSize;
delta2 = para.patchSize;
A = Phi2patches_fast(C, para.patchSize(1), para.patchSize(2), para.T, delta1, delta2);                %最外层TwIST算法进行重构的测量矩阵

patchSize_w = para.patchSize_w;

load(['..\恢复结果\' orig_str '\wiener\block',num2str(patchSize_w(1)),'_T',num2str(para.T),'_F',num2str(para.T),'_二值_wiener.mat']); %para.M

X_wiener=Xrecon;

for m = 1:para.M    
    X_col=video2patches_fast(X_wiener(:,:,(m-1)*para.T+1:m*para.T), para.patchSize(1), para.patchSize(2), para.T, delta1, delta2); 
    y = video2patches_fast(Y(:,:,m), para.patchSize(1), para.patchSize(2), 1, delta1, delta2);    
    
    [d,Nt] = size(y);
    X=[];
    addpath('TwIST\TwIST_TV_CS');      
     for i = 1:Nt
        RR=A{i};
        y_column=y(:,i);       
        hR = @(x) RR*x;
        hRt = @(x) RR'*x;
        
        % denoising function;
        tv_iters = 5;
        Psi = @(x,th,T)  tvdenoise_3d(x,2/th,T,tv_iters);
        % TV regularizer;
        Phi = @(x,T) TVnorm_3d(x,T);
        
        % regularization parameter
        % tau = 0.1*max(abs(hRt(y)));
        tau =0.001;       
        % TwIST parameters 
        lambda1 =0.00001; %0.001
        
        % stopping theshold
        tolA = 1e-5;
        
%         迭代的初值
        x_pre=X_col(:,i);
        
        % -- TwIST ---------------------------
        % stop criterium:  the relative change in the objective function
        % falls below 'ToleranceA'
        [x_twist,x_debias_twist,obj_twist,...
            times_twist,debias_start_twist,mse]= ... 
            TwIST_convert(y_column,hR,tau,para.T,x_pre, ...
            'Lambda', lambda1, ...
            'Debias',0,...
            'AT', hRt, ...
            'Psi', Psi, ...
            'Phi',Phi, ...
            'Monotone',1,...
            'Sparse', 1,...
            'Initialization',0,...
            'StopCriterion',1,...
            'ToleranceA',tolA,...
            'Verbose', 1);
%         x_pre=x_twist;
        x_r=x_twist;
        X=[X x_r];   
     end
     before_X=X;
    Xpre(:,:,(m-1)*para.T+1:m*para.T) = patches2video_fast(X, para.row, para.col, para.patchSize(1), para.patchSize(1), para.T, delta1, delta2); 
    %%将复原后的视频存储到Xpre矩阵中。1帧恢复成T帧。
end