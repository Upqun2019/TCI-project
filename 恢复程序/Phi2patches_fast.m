%% 将调制模板C以8*8大小分块处理，并取出每一张模板中对应块的元素，对角化组成测量矩阵。有考虑块的重合。
%% Phi是原始测量模板数据；n1，n2是分块的大小；delta1, delta2是块左右重叠的大小。所得矩阵Phi2是一个多维矩阵。
%% 每一个矩阵是相应块对角化横向拼接后组成的测量矩阵。

function Phi2 = Phi2patches_fast(Phi, n1, n2, T, delta1, delta2)
% This function convert the 3D cube into 2D patches
% Jianbo Yang
% 06/19/2013

temp = image2patches_fast(Phi(:,:,1), n1, n2, delta1, delta2);      %temp大小为64*3969，一种将模板取出排列的方式。考虑到了分块以及块的重合。
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
        temp = [temp sparse(diag(M_patches(:,i,t)))];     %维度太大导致内存不够，所以使用稀疏矩阵，也可以降低精度
    end
    Phi2{i} = temp;
end
return;
