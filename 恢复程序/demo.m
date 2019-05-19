%% �ع��㷨��������
%Author:Qun Zhou  E-mail:zhouuqn817@126.com
%%2019/05/16
%Reference:Zhou Q, Ke J, Lam E Y. Near-infrared temporal compressive imaging for video[J]. Optics Letters, 2019,44(7):1702-1705.

clear all;clc;

para.mask='��ֵ';  %��ֵ    ��˹
orig_str='Billiards';                                                                       %gun     Lighter       Billiards 
training_str='combine';                                                                     %gun     Lighter      Billiards    combine

%���ز�����Ƶ�Ͳ���ģ��
load(strcat('..\..\����ģ�弰�۲���Ƶ\',orig_str,'\data_global_',orig_str,'_',para.mask));

[para.row,para.col,para.M]=size(Y);
[para.row, para.col, para.T] = size(C);

b_size=128;                                                                                 %���С��GMM�㷨��ѡȡ���СΪ��8��TwIST�㷨��ѡȡ���СΪ��128  
para.M=8;            
para.patchSize = [b_size b_size para.T];
 
%ѵ��ģ�Ͳ���
para.filename_training=strcat('..\..\ѵ��\ѵ�����\combine\block',num2str(b_size),'_T',num2str(para.T),'\');      %ѡ��ѵ��ģ��·��
para.training_model=strcat('model');  
                         

% �ָ��㷨
Algorithm =3;
tic
switch Algorithm
    case 1 % learn the GMM offline
        method='GMM';
        para.method = 'GMM_offline';                                                %GMM_offline   wiener_GMM_offline
        if strcmp(para.method,'wiener_GMM_offline')                                 %ѡ����㷨ʱ����Ҫ����ѡģ�Ͳ�����Ӧ����ʹ��wiener�㷨Ԥ�ع������Ƶ��Ϊѵ����Ƶѵ����Ľ��)
            para.filename_training=strcat('..\..\ѵ��\ѵ�����\wiener\block',num2str(b_size),'_T',num2str(para.T),'\');
        end
        addpath('GMM/')
        Xrecon = interface_GMM(Y,C,para);
        rmpath('GMM/')
    case 2 % learn the GMM offline and update the GMM online. Caution: it takes longer time.
        method='GMM';
        para.method = 'GMM_online';
        para.training_data0='data0';                                                   % ѵ���ı�ģ�Ͳ�����ѡ��ģ�͵�ԭʼ���ݼ�
        addpath('GMM/')
        Xrecon = interface_GMM_online(Y,C,para);
        rmpath('GMM/')
    case 3
        method='TwIST';
        addpath('TwIST/')
        para.method = 'TwIST_online';                                                   %���Ըı��ʼ����ʽ��1:TwIST_zero   2: x_0=A'*y    3:TwIST_online
        Xrecon =interface_TwIST_TV(Y,C,para);
        rmpath('TwIST/')
    case 4                                                                              %����ʹ��wiener�����㷨�����ع�
        method='wiener';
        addpath('wiener/')
        para.method = 'wiener';                
        Xrecon =interface_wiener(Y,C,para);
        rmpath('wiener/')
    case 5                                                                               %ʹ��wiener�㷨Ԥ�ع��Ľ����Ϊ�����ĳ�ʼֵ
        method='wiener_TwIST';
        addpath('TwIST/')
        para.method = 'wiener_TwIST';                        
        b_w=8;                                                                           %b_w��ֵ��С��ʹ��wiener�㷨����Ԥ�ع�ʱ�Ŀ��С
        para.patchSize_w=[b_w b_w para.T];
        Xrecon = interface_wiener_TwIST(Y,C,para,orig_str);
        rmpath('TwIST/')
end
%% ֱ�ӽ��ָ������ʾ������
para.time=toc;
filename=strcat('..\�ָ����\',orig_str,'\',method,'\');     %������·��
[PSNR, SSIM,ave_psnr] = saveResults(Xrecon,Xtst,Y,para,filename); 

temp=mean(PSNR);
figure;plot(PSNR);
figname=strcat('PSNR of different reconstructed frames(mean=',num2str(temp),')');
title(figname,'fontsize',20); 
xlabel('Frames');ylabel('PSNR/dB');
set(gca,'fontsize',20);
scrsz = get(0,'ScreenSize');

