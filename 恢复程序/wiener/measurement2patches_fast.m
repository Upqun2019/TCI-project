%%��ÿһ֡�۲�ͼ��Y�Էֿ��Ĺ۲�ֵ˳������
%%�õ���X��һ����ά����ÿһ����Y��֡�ж�Ӧ��Ԫ����ɡ�

function y = measurement2patches_fast(Y,patchSize)
n1=patchSize(1);n2=patchSize(2);T=patchSize(3);
delta1 = patchSize;
delta2 = patchSize;
y = image2patches_fast(Y, n1, n2, delta1, delta2);


