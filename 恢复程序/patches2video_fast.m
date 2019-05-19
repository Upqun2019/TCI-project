function V = patches2video_fast(X, Row,Col, n1, n2, T, delta1, delta2)
% This function convert the 2D patches into 3D cube
% Jianbo Yang
% 06/19/2013

d = n1*n2;
V = zeros(Row,Col,T);
for t = 1:T
    temp = X((t-1)*d+1:t*d,:);   %大小为（8^2）*T，每一帧恢复后X的T帧中，对应块变为列向量后组成的矩阵
    V(:,:,t) = patches2image_fast(temp, Row, Col, delta1, delta2);
end
