%% �õ���Ӧ��Ƶ������ؾ���

function R = R_wiener(X,para)

n = prod( para.patchSize);
delta1 =  para.patchSize;
delta2 =  para.patchSize;
% delta1 =  para.patchSize/2;                                               % overlapping block
% delta2 =  para.patchSize/2;

M=20;                                                                       %ѵ����Ƶ��֡��Ϊ(M*para.T)����һ����Χ�ڣ���ֵԽ��Խ�á�
X_column=[];
for m=1:M
    %�ɵ�m֡��Ӧ��ԭ��Ƶ��T֡X�����䰴����Ƶ��ķ��гɾ���
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
