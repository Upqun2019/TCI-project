%% ������ģ��C��8*8��С�ֿ鴦����ȡ��ÿһ��ģ���ж�Ӧ���Ԫ�أ��Խǻ���ɲ��������п��ǿ���غϡ�
%% Phi��ԭʼ����ģ�����ݣ�n1��n2�Ƿֿ�Ĵ�С��delta1, delta2�ǿ������ص��Ĵ�С�����þ���Phi2��һ����ά����
%% ÿһ����������Ӧ��Խǻ�����ƴ�Ӻ���ɵĲ�������

function Phi2 = Phi2patches_fast(Phi, n1, n2, T, delta1, delta2)
% This function convert the 3D cube into 2D patches
% Jianbo Yang
% 06/19/2013

temp = image2patches_fast(Phi(:,:,1), n1, n2, delta1, delta2);      %temp��СΪ64*3969��һ�ֽ�ģ��ȡ�����еķ�ʽ�����ǵ��˷ֿ��Լ�����غϡ�
[d1, d2] = size(temp);
M_patches = zeros(d1,d2,T);
for t = 1:T
    M_patches(:,:,t) = image2patches_fast(Phi(:,:,t), n1, n2, delta1, delta2);
end
n = size(M_patches,2);
Phi2 = cell(n,1); 
for i = 1:n
    temp = [];
    for t = 1:T
        temp = [temp sparse(diag(M_patches(:,i,t)))];     %ά��̫�����ڴ治��������ʹ��ϡ�����Ҳ���Խ��;���
    end
    Phi2{i} = temp;
end
return;
