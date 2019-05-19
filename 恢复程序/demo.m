%% 重构算法的主函数
%Author:Qun Zhou  E-mail:zhouuqn817@126.com
%%2019/05/16
%Reference:Zhou Q, Ke J, Lam E Y. Near-infrared temporal compressive imaging for video[J]. Optics Letters, 2019,44(7):1702-1705.

clear all;clc;

para.mask='二值';  %二值    高斯
orig_str='Billiards';                                                                       %gun     Lighter       Billiards 
training_str='combine';                                                                     %gun     Lighter      Billiards    combine

%下载测量视频和测量模板
load(strcat('..\..\测量模板及观测视频\',orig_str,'\data_global_',orig_str,'_',para.mask));

[para.row,para.col,para.M]=size(Y);
[para.row, para.col, para.T] = size(C);

b_size=128;                                                                                 %块大小，GMM算法中选取块大小为：8；TwIST算法中选取块大小为：128  
para.M=8;            
para.patchSize = [b_size b_size para.T];
 
%训练模型参数
para.filename_training=strcat('..\..\训练\训练结果\combine\block',num2str(b_size),'_T',num2str(para.T),'\');      %选择训练模型路径
para.training_model=strcat('model');  
                         

% 恢复算法
Algorithm =3;
tic
switch Algorithm
    case 1 % learn the GMM offline
        method='GMM';
        para.method = 'GMM_offline';                                                %GMM_offline   wiener_GMM_offline
        if strcmp(para.method,'wiener_GMM_offline')                                 %选择该算法时，需要对重选模型参数（应该是使用wiener算法预重构后的视频作为训练视频训练后的结果)
            para.filename_training=strcat('..\..\训练\训练结果\wiener\block',num2str(b_size),'_T',num2str(para.T),'\');
        end
        addpath('GMM/')
        Xrecon = interface_GMM(Y,C,para);
        rmpath('GMM/')
    case 2 % learn the GMM offline and update the GMM online. Caution: it takes longer time.
        method='GMM';
        para.method = 'GMM_online';
        para.training_data0='data0';                                                   % 训练改变模型参数，选择模型的原始数据集
        addpath('GMM/')
        Xrecon = interface_GMM_online(Y,C,para);
        rmpath('GMM/')
    case 3
        method='TwIST';
        addpath('TwIST/')
        para.method = 'TwIST_online';                                                   %可以改变初始化方式：1:TwIST_zero   2: x_0=A'*y    3:TwIST_online
        Xrecon =interface_TwIST_TV(Y,C,para);
        rmpath('TwIST/')
    case 4                                                                              %单独使用wiener线性算法进行重构
        method='wiener';
        addpath('wiener/')
        para.method = 'wiener';                
        Xrecon =interface_wiener(Y,C,para);
        rmpath('wiener/')
    case 5                                                                               %使用wiener算法预重构的结果作为迭代的初始值
        method='wiener_TwIST';
        addpath('TwIST/')
        para.method = 'wiener_TwIST';                        
        b_w=8;                                                                           %b_w的值大小是使用wiener算法进行预重构时的块大小
        para.patchSize_w=[b_w b_w para.T];
        Xrecon = interface_wiener_TwIST(Y,C,para,orig_str);
        rmpath('TwIST/')
end
%% 直接将恢复结果显示、保存
para.time=toc;
filename=strcat('..\恢复结果\',orig_str,'\',method,'\');     %保存结果路径
[PSNR, SSIM,ave_psnr] = saveResults(Xrecon,Xtst,Y,para,filename); 

temp=mean(PSNR);
figure;plot(PSNR);
figname=strcat('PSNR of different reconstructed frames(mean=',num2str(temp),')');
title(figname,'fontsize',20); 
xlabel('Frames');ylabel('PSNR/dB');
set(gca,'fontsize',20);
scrsz = get(0,'ScreenSize');

