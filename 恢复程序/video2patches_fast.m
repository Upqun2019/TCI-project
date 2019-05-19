%%���۲���ƵY�Էֿ��Ĺ۲�ֵ˳������
%%�õ���X��һ����ά����ÿһ����Y����֡�ж�Ӧ��Ԫ����ɡ�

function X = video2patches_fast(V, n1, n2, T, delta1, delta2)
% This function convert the 3D cube into 2D patches

temp = image2patches_fast(V(:,:,1), n1, n2, delta1, delta2);
[d1, d2] = size(temp);
M_patches = zeros(d1,d2,T); 
for t = 1:T
    M_patches(:,:,t) = image2patches_fast(V(:,:,t), n1, n2, delta1, delta2);
end
n = size(M_patches,2);
X = zeros(n1*n2*T,n);
for i = 1:n
    temp = squeeze(M_patches(:,i,:));    % squeeze������һ����Ϻ���������������еĶ�Ӧ�ĵ�i��ȡ��˳����ϡ�
    X(:,i) = temp(:);
end
