function y = TVnorm(x,w)

if size(x,2)==1
    x=convert(x,w);
end

y = sum(sum(sqrt(diffh(x).^2+diffv(x).^2)));



function x_tv=convert(x,w)     %% ��һάϡ���ź�ת��ɿռ��򣬲���Ϊ��ά�ź�
n=sqrt(size(x,1)*size(x,2));
temp=w'*x; 
temp=reshape(temp,n,n);
x_tv=temp';