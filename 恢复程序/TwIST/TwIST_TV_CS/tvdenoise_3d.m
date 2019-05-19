function u = tvdenoise_3d(f,lambda,T,iters)
%TVDENOISE  Total variation grayscale and color image denoising
%   u = TVDENOISE(f,lambda) denoises the input image f.  The smaller
%   the parameter lambda, the stronger the denoising.
%
%   The output u approximately minimizes the Rudin-Osher-Fatemi (ROF)
%   denoising model
%
%       Min  TV(u) + lambda/2 || f - u ||^2_2,
%        u
%
%   where TV(u) is the total variation of u.  If f is a color image (or any
%   array where size(f,3) > 1), the vectorial TV model is used,
%
%       Min  VTV(u) + lambda/2 || f - u ||^2_2.
%        u
%
%   TVDENOISE(...,Tol) specifies the stopping tolerance (default 1e-2).
%
%   The minimization is solved using Chambolle's method,
%      A. Chambolle, "An Algorithm for Total Variation Minimization and
%      Applications," J. Math. Imaging and Vision 20 (1-2): 89-97, 2004.
%   When f is a color image, the minimization is solved by a generalization
%   of Chambolle's method,
%      X. Bresson and T.F. Chan,  "Fast Minimization of the Vectorial Total
%      Variation Norm and Applications to Color Image Processing", UCLA CAM
%      Report 07-25.
%
%   Example:
%   f = double(imread('barbara-color.png'))/255;
%   f = f + randn(size(f))*16/255;
%   u = tvdenoise(f,12);
%   subplot(1,2,1); imshow(f); title Input
%   subplot(1,2,2); imshow(u); title Denoised

% Pascal Getreuer 2007-2008
%  Modified by Jose Bioucas-Dias  & Mario Figueiredo 2010 
%  (stopping rule: iters)
%   

% iscon=0;
% if size(f,2)==1    
%     f=convert(f,w);
%     iscon=1;
% end

iscon=0;
if size(f,2)==1    
    f=convert_to3d(f,T) ;
    iscon=1;
end

if nargin < 5
    Tol = 1e-2;
end

if lambda < 0
    error('Parameter lambda must be nonnegative.');
end

dt = 0.25;

N = size(f);
id = [2:N(1),N(1)];
iu = [1,1:N(1)-1];
ir = [2:N(2),N(2)];
il = [1,1:N(2)-1];
p1 = zeros(size(f));
p2 = zeros(size(f));
divp = zeros(size(f));
lastdivp = ones(size(f));


if length(N) == 2           % TV denoising
    %while norm(divp(:) - lastdivp(:),inf) > Tol
    for i=1:iters
        lastdivp = divp;
        z = divp - f*lambda;
        z1 = z(:,ir) - z;
        z2 = z(id,:) - z;
        denom = 1 + dt*sqrt(z1.^2 + z2.^2);
        p1 = (p1 + dt*z1)./denom;
        p2 = (p2 + dt*z2)./denom;
        divp = p1 - p1(:,il) + p2 - p2(iu,:);
    end
elseif length(N) == 3       % Vectorial TV denoising
    
    repchannel = ones(N(3),1);

    %while norm(divp(:) - lastdivp(:),inf) > Tol
    for i=1:iters
        lastdivp = divp;
        z = divp - f*lambda;
        z1 = z(:,ir,:) - z;
        z2 = z(id,:,:) - z;
        denom = 1 + dt*sqrt(sum(z1.^2 + z2.^2,3));
        denom = denom(:,:,repchannel);
        p1 = (p1 + dt*z1)./denom;
        p2 = (p2 + dt*z2)./denom;
        divp = p1 - p1(:,il,:) + p2 - p2(iu,:,:);
    end
end

u = f - divp/lambda;

if iscon==1
    u=convert_ni_to3d(u,T);
end

% function x_tv=convert(x,w)     %% 将一维稀疏信号转变成空间域，并变为二维信号
% n=sqrt(size(x,1)*size(x,2));
% temp=w'*x; 
% temp=reshape(temp,n,n);
% x_tv=temp';
% 
% function x_sparse=convert_ni(x_tv,w)    %%与上面的convert函数相反
% temp=x_tv';x=temp(:);
% x_sparse=w*x; 

function x_tv=convert_to3d(x,T)     %% 将一维稀疏信号转变成空间域，并变为三维信号
size_block=size(x,1)*size(x,2)/T;
n=sqrt(size_block);
x_tv=zeros(n,n,T);

for i=1:T
   temp=x(size_block*(i-1)+1:size_block*i);
     x_tv(:,:,i)=reshape(temp,n,n);
end

function  x_column=convert_ni_to3d(x_tv,T)    %%与上面的convert_to3d函数相反
x_column=[];
for i=1:T
    temp=x_tv(:,:,i);
    temp=temp(:);
    x_column=[x_column;temp];
end




