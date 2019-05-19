%%将每一帧观测图像Y以分块后的观测值顺序排列
%%得到的X是一个二维矩阵，每一列由Y该帧中对应块元素组成。

function y = measurement2patches_fast(Y,patchSize)
n1=patchSize(1);n2=patchSize(2);T=patchSize(3);
delta1 = patchSize;
delta2 = patchSize;
y = image2patches_fast(Y, n1, n2, delta1, delta2);


