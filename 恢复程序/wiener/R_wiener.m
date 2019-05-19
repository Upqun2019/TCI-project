%% 得到相应视频的自相关矩阵

function R = R_wiener(X,para)

n = prod( para.patchSize);
delta1 =  para.patchSize;
delta2 =  para.patchSize;
% delta1 =  para.patchSize/2;                                               % overlapping block
% delta2 =  para.patchSize/2;

M=20;                                                                       %训练视频的帧数为(M*para.T)。在一定范围内，该值越大越好。
X_column=[];
for m=1:M
    %由第m帧对应的原视频的T帧X，将其按照视频块的分列成矩阵
    X_temp=X(:,:,1+(m-1)*para.T:m*para.T);
    X_col=[];
    for i=1:para.T
        temp=measurement2patches_fast(X_temp(:,:,i),para.patchSize);
        X_col=[X_col;temp];
    end
    X_column=[X_column X_col];
end
[N num]=size(X_column);

temp=mean(X_column,2);
mean_x=repmat(temp,1,num);
X_center=X_column-mean_x;

temp=zeros(N,N);
for i=1:num
    temp=temp+X_center(:,i)*X_center(:,i)';
end
R=temp./num;
% mkdir(strcat('R_mask\R',num2str(para.patchSize(1)),'_T',num2str(para.T)));
% save(strcat('R_mask\R',num2str(para.patchSize(1)),'_T',num2str(para.T)),'R');
