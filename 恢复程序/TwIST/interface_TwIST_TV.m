function Xpre = interface_TwIST_TV(Y,C,para)

ini_md=3;                                                                                         %choose initialization method

n = prod(para.patchSize);
delta1 = para.patchSize;
delta2 = para.patchSize;
A = Phi2patches_fast(C, para.patchSize(1), para.patchSize(2), para.T, delta1, delta2);             %得到以块处理的测量矩阵

for m = 1:para.M
    y = video2patches_fast(Y(:,:,m), para.patchSize(1), para.patchSize(2), 1, delta1, delta2);     %将观测到的视频，一帧一帧恢复。  
    addpath('TwIST\TwIST_TV_CS');
    [d,Nt] = size(y);
    X=[];
    
    for i = 1:Nt
        R=A{i};
        y_column=y(:,i);
        hR = @(x) R*x;
        hRt = @(x) R'*x;
        
        % denoising function;
        tv_iters = 5;
        Psi = @(x,th,T)  tvdenoise_3d(x,2/th,T,tv_iters);
        % TV regularizer;
        Phi = @(x,T) TVnorm_3d(x,T);
        
        % regularization parameter
        % tau = 0.1*max(abs(hRt(y)));
        tau =0.009;   %0.0015
        % TwIST parameters
        lambda1 =0.00001;
        % stopping theshold
        tolA = 1e-5;
        
        %迭代的初值
        if ini_md==1                                               %方式1
            x_pre=zeros(n,1);
        elseif ini_md==2                                           %方式2
            x_pre=hRt(y_column);
        elseif ini_md==3                                           %方式3
            if m==1
                x_pre=zeros(n,1);
            else
                x_pre=before_X(:,i);
            end
        end
        
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
        x_r=x_twist;
        X=[X x_r];
    end
    if  ini_md==3                                                   
        before_X=X;
    end
    Xpre(:,:,(m-1)*para.T+1:m*para.T) = patches2video_fast(X, para.row, para.col, para.patchSize(1), para.patchSize(1), para.T, delta1, delta2);
    %%将复原后的视频存储到Xpre矩阵中。1帧恢复成T帧。
end

