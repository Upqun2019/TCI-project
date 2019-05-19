function [PSNR,SSIM,ave_psnr] = saveResults(Xrecon,Xtst,Y,para,filename)
Row=para.row; Col=para.col;
M=para.M; T=para.T;
b_size=para.patchSize(1);

% Xrecon = zeros(Row,Col,M*T);
temp=0;
for k=1:(M*T)
    PSNR(k) = psnr(Xrecon(:,:,k), Xtst(:,:,k));
    SSIM(k) = ssim(Xrecon(:,:,k), Xtst(:,:,k));
    temp=temp+(Xrecon(:,:,k)-Xtst(:,:,k)).^2;
%     MSE2=sum(sum((Xrecon(:,:,k)-Xtst(:,:,k)).^2))/(Row*Col);
%     ave_psnr2(k)=10*log10(1/MSE2);
end
temp=temp/(Row*Col)/(M*T);
MSE=sum(sum(temp));
ave_psnr=10*log10(1/MSE);
time=para.time;
Algorim=para.method;
savename = [filename 'block' num2str(b_size) '_T' num2str(T) '_F' num2str(M) '_' para.mask '_' Algorim];
format short
% save(savename, '-v7.3', 'Xrecon','Xpre','PSNR','time','SSIM');
save(savename, '-v7.3', 'Xtst','Xrecon','PSNR','SSIM','ave_psnr','time');

writerObj = VideoWriter([savename '.mp4'],'MPEG-4');
writerObj.FrameRate = 12;
open(writerObj);
scrsz = get(0,'ScreenSize');
fig=figure('Position',[50 100 floor(scrsz(3)*0.8) floor(scrsz(4)*0.6)]);
for nF=1:(M*T)
    subplot(1,3,1);
%     imshow(uint8(Xtst(:,:,nF)));
    imshow(Xtst(:,:,nF));
    title(['Original Video, Frame:' num2str(nF)],'fontsize',20);
    
    
    subplot(1,3,2);
    imshow(Y(:,:,ceil(nF/T))/max(max(Y(:,:,ceil(nF/T)))));
    title(['Measurement video, Frame: ' num2str(ceil(nF/T))],'fontsize',20);
    
    subplot(1,3,3);
%     imshow(uint8(Xrecon(:,:,nF)));
    imshow(Xrecon(:,:,nF));
    title(['Recon video, PSNR: ' num2str(PSNR(nF))],'fontsize',20);
        
    %pause(0.02);
    frame = getframe(fig);
    writeVideo(writerObj,frame);
end
close(writerObj);

% figure;
% i=M*(T/2)
% d=i+4;
% % subplot(1,3,1);
% %     imshow(uint8(Xtst(:,:,nF)));
% imshow(Xtst(:,:,d));
% title(['Original Video, Frame:' num2str(d)],'fontsize',20);
% rectangle('position',[95,125,40,40],'edgecolor','r');
% 
% figure;
% % subplot(1,3,2);
% imshow(Y(:,:,ceil(i/T))/max(max(Y(:,:,ceil(i/T)))));
% title(['Measurement video, Frame: ' num2str(ceil(i/T))],'fontsize',20);
% figure;
% % subplot(1,3,3);
% %     imshow(uint8(Xrecon(:,:,nF)));
% imshow(Xrecon(:,:,d));
% % title(['Recon video, PSNR: ' num2str(PSNR(d))],'fontsize',20);
% %  title(['Reconstructed frame of Offline, PSNR: ' num2str(PSNR(d))],'fontsize',20);
%  title(['Reconstructed frame of Online, PSNR: ' num2str(PSNR(d))],'fontsize',20);
%  rectangle('position',[95,125,40,40],'edgecolor','r');
%  
%  figure;
%  imshow(Xrecon(125:125+50,95:95+50,d)*5);
% %  title('Local pattern of Offline');
%  title('Local pattern of Online');

       


