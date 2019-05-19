function y = TVnorm_3d(x,T)

x=convert_to3d(x,T) ;
y=0;
for i=1:T
    temp = sum(sum(sqrt(diffh(x(:,:,i)).^2+diffv(x(:,:,i)).^2)));
    y=y+temp;
end   

function x_tv=convert_to3d(x,T)     %% 将一维稀疏信号转变成空间域，并变为三维信号
size_block=size(x,1)*size(x,2)/T;
n=sqrt(size_block);
x_tv=zeros(n,n,T);

for i=1:T
   temp=x(size_block*(i-1)+1:size_block*i);
     x_tv(:,:,i)=reshape(temp,n,n);
end