close all; clear; clc;

load('err3'); load('err6'); load('err10'); load('err15'); load('err21');

dist_1 = [norm(err3(:,500)) norm(err6(:,500)) norm(err10(:,500)) ...
        norm(err15(:,500)) norm(err21(:,500))];
dist_2 = [norm(err3(:,750)) norm(err6(:,750)) norm(err10(:,750)) ...
        norm(err15(:,750)) norm(err21(:,750))];
dist_3 = [norm(err3(:,end-1)) norm(err6(:,end-1)) norm(err10(:,end-1)) ...
        norm(err15(:,end-1)) norm(err21(:,end-1))];

P = [2 5 9 14 20];    
figure, hold on;
plot(P,log(dist_1),'b-*','Linewidth',2);
plot(P,log(dist_2),'g-*','Linewidth',2);
plot(P,log(dist_3),'r-*','Linewidth',2);
hold off; xlabel('P'); setleg(legend('i = 0.5*N','i = 0.75*N','i = N'));
    
figure, 
subplot(2,2,1)
imagesc(abs(err6)); colorbar; xlabel('i'); 
title('P = 5'); ylabel('k','Rotation',0);
subplot(2,2,2)
imagesc(abs(err10)); colorbar; xlabel('i'); 
title('P = 9'); ylabel('k','Rotation',0);
subplot(2,2,3)
imagesc(abs(err15)); colorbar; xlabel('i'); 
title('P = 14'); ylabel('k','Rotation',0);
subplot(2,2,4)
imagesc(abs(err21)); colorbar; xlabel('i'); 
title('P = 20'); ylabel('k','Rotation',0);