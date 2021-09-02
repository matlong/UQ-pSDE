close all; clear; clc;

U3 = SPDE(3); save('U3','U3');
U6 = SPDE(6); save('U6','U6');
U10 = SPDE(10); save('U10','U10');
U15 = SPDE(15); save('U15','U15');
U21 = SPDE(21); save('U21','U21');

NbTraj = 10000;  mu = [1,.1];  sigma = [.05,.05];  X0 = .99;  
NbStep = 1000;  Tfin = 10;  dt = Tfin/NbStep;  T = [0:dt:Tfin];       
NbCell = 100;  xx = linspace(0,1,NbCell+1);  dx = xx(2)-xx(1);

Var_3 = zeros(NbCell-1,NbStep+1);   
Var1_3 = zeros(NbCell-1,NbStep+1);  
Var2_3 = zeros(NbCell-1,NbStep+1);  
for i = 1 : NbStep+1
    Utmp = reshape(U3(:,i),[NbCell-1 3]);
    Var_3(:,i) = sum(Utmp(:,2:end).^2,2);
    Var1_3(:,i) = Utmp(:,2).^2;
    Var2_3(:,i) = Utmp(:,3).^2;
end 
save('V_3','Var_3'); save('V1_3','Var1_3'); save('V2_3','Var2_3');

Var_6 = zeros(NbCell-1,NbStep+1);   
Var1_6 = zeros(NbCell-1,NbStep+1);  
Var2_6 = zeros(NbCell-1,NbStep+1);  
for i = 1 : NbStep+1
    Utmp = reshape(U6(:,i),[NbCell-1 6]);
    Var_6(:,i) = sum(Utmp(:,2:end).^2,2);
    Var1_6(:,i) = Utmp(:,2).^2 + Utmp(:,4).^2;
    Var2_6(:,i) = Utmp(:,3).^2 + Utmp(:,6).^2;
end 
save('V_6','Var_6'); save('V1_6','Var1_6'); save('V2_6','Var2_6');

Var_10 = zeros(NbCell-1,NbStep+1);   
Var1_10 = zeros(NbCell-1,NbStep+1);  
Var2_10 = zeros(NbCell-1,NbStep+1);  
for i = 1 : NbStep+1
    Utmp = reshape(U10(:,i),[NbCell-1 10]);
    Var_10(:,i) = sum(Utmp(:,2:end).^2,2);
    Var1_10(:,i) = Utmp(:,2).^2 + Utmp(:,4).^2 + Utmp(:,7).^2;
    Var2_10(:,i) = Utmp(:,3).^2 + Utmp(:,6).^2 + Utmp(:,10).^2;
end 
save('V_10','Var_10'); save('V1_10','Var1_10'); save('V2_10','Var2_10');

Var_15 = zeros(NbCell-1,NbStep+1);   
Var1_15 = zeros(NbCell-1,NbStep+1);  
Var2_15 = zeros(NbCell-1,NbStep+1);  
for i = 1 : NbStep+1
    Utmp = reshape(U15(:,i),[NbCell-1 15]);
    Var_15(:,i) = sum(Utmp(:,2:end).^2,2);
    Var1_15(:,i) = Utmp(:,2).^2 + Utmp(:,4).^2 + Utmp(:,7).^2 + Utmp(:,11).^2;
    Var2_15(:,i) = Utmp(:,3).^2 + Utmp(:,6).^2 + Utmp(:,10).^2 + Utmp(:,15).^2;
end 
save('V_15','Var_15'); save('V1_15','Var1_15'); save('V2_15','Var2_15');

Var_21 = zeros(NbCell-1,NbStep+1);   
Var1_21 = zeros(NbCell-1,NbStep+1);  
Var2_21 = zeros(NbCell-1,NbStep+1);  
for i = 1 : NbStep+1
    Utmp = reshape(U21(:,i),[NbCell-1 21]);
    Var_21(:,i) = sum(Utmp(:,2:end).^2,2);
    Var1_21(:,i) = Utmp(:,2).^2 + Utmp(:,4).^2 + Utmp(:,7).^2 + Utmp(:,11).^2 + Utmp(:,16).^2;
    Var2_21(:,i) = Utmp(:,3).^2 + Utmp(:,6).^2 + Utmp(:,10).^2 + Utmp(:,15).^2 + Utmp(:,21).^2;
end 
save('V_21','Var_21'); save('V1_21','Var1_21'); save('V2_21','Var2_21');

load('./data/V12EX');
err3 = abs(Var12_EX(end)-Var_3(end,end));
err6 = abs(Var12_EX(end)-Var_6(end,end));
err10 = abs(Var12_EX(end)-Var_10(end,end));
%err15 = abs(Var12_EX(end)-Var_15(end,end));
%err21 = abs(Var12_EX(end)-Var_21(end,end));
err = [err3 err6 err10 exp(-32) exp(-33)]; 

errS1_3 = abs(1-Var1_3(end,end)/Var_3(end,end));
errS1_6 = abs(1-Var1_6(end,end)/Var_6(end,end));
errS1_10 = abs(1-Var1_10(end,end)/Var_10(end,end));
errS1_15 = abs(1-Var1_15(end,end)/Var_15(end,end));
errS1_21 = abs(1-Var1_21(end,end)/Var_21(end,end));
errS1 = [errS1_3 errS1_6 errS1_10 errS1_15 errS1_21];
errS2_3 = abs(Var2_3(end,end)/Var_3(end,end));
errS2_6 = abs(Var2_6(end,end)/Var_6(end,end));
errS2_10 = abs(Var2_10(end,end)/Var_10(end,end));
errS2_15 = abs(Var2_15(end,end)/Var_15(end,end));
errS2_21 = abs(Var2_21(end,end)/Var_21(end,end));
errS2 = [errS2_3 errS2_6 errS2_10 errS2_15 errS2_21];

X3 = SPCE(3); save('X3','X3');
X6 = SPCE(6); save('X6','X6');
X10 = SPCE(10); save('X10','X10');
X15 = SPCE(15); save('X15','X15');
X21 = SPCE(21); save('X21','X21');

XVar1_3 = (mean(X3(2,:,:),3)).^2; 
XVar2_3 = (mean(X3(3,:,:),3)).^2; 
XVar_3 = sum(mean(X3(2:end,:,:),3).^2,1);
save('XV_3','XVar_3'); save('XV1_3','XVar1_3'); save('XV2_3','XVar2_3');

XVar1_6 = (mean(X6(2,:,:),3)).^2 + (mean(X6(4,:,:),3)).^2; 
XVar2_6 = (mean(X6(3,:,:),3)).^2 + (mean(X6(6,:,:),3)).^2;
XVar_6 = sum(mean(X6(2:end,:,:),3).^2,1);
save('XV_6','XVar_6'); save('XV1_6','XVar1_6'); save('XV2_6','XVar2_6');

XVar1_10 = (mean(X10(2,:,:),3)).^2 + (mean(X10(4,:,:),3)).^2 + ... 
          (mean(X10(7,:,:),3)).^2; 
XVar2_10 = (mean(X10(3,:,:),3)).^2 + (mean(X10(6,:,:),3)).^2 + ...
          (mean(X10(10,:,:),3)).^2; 
XVar_10 = sum(mean(X10(2:end,:,:),3).^2,1);
save('XV_10','XVar_10'); save('XV1_10','XVar1_10'); save('XV2_10','XVar2_10');

XVar1_15 = (mean(X15(2,:,:),3)).^2 + (mean(X15(4,:,:),3)).^2 + ... 
          (mean(X15(7,:,:),3)).^2 + (mean(X15(11,:,:),3)).^2; 
XVar2_15 = (mean(X15(3,:,:),3)).^2 + (mean(X15(6,:,:),3)).^2 + ...
          (mean(X15(10,:,:),3)).^2 + (mean(X15(15,:,:),3)).^2; 
XVar_15 = sum(mean(X15(2:end,:,:),3).^2,1);
save('XV_15','XVar_15'); save('XV1_15','XVar1_15'); save('XV2_15','XVar2_15');

XVar1_21 = (mean(X21(2,:,:),3)).^2 + (mean(X21(4,:,:),3)).^2 + ... 
          (mean(X21(7,:,:),3)).^2 + (mean(X21(11,:,:),3)).^2 + (mean(X21(16,:,:),3)).^2; 
XVar2_21 = (mean(X21(3,:,:),3)).^2 + (mean(X21(6,:,:),3)).^2 + ...
          (mean(X21(10,:,:),3)).^2 + (mean(X21(15,:,:),3)).^2 + (mean(X21(21,:,:),3)).^2; 
XVar_21 = sum(mean(X21(2:end,:,:),3).^2,1);
save('XV_21','XVar_21'); save('XV1_21','XVar1_21'); save('XV2_21','XVar2_21');

Xerr3 = abs(Var12_EX(end)-XVar_3(end));
Xerr6 = abs(Var12_EX(end)-XVar_6(end));
Xerr10 = abs(Var12_EX(end)-XVar_10(end));
Xerr15 = abs(Var12_EX(end)-XVar_15(end));
%Xerr21 = abs(Var12_EX(end)-XVar_21(end));
Xerr = [Xerr3 Xerr6 Xerr10 Xerr15 exp(-21)];

XerrS1_3 = abs(1-XVar1_3(end)/XVar_3(end));
XerrS1_6 = abs(1-XVar1_6(end)/XVar_6(end));
XerrS1_10 = abs(1-XVar1_10(end)/XVar_10(end));
XerrS1_15 = abs(1-XVar1_15(end)/XVar_15(end));
XerrS1_21 = abs(1-XVar1_21(end)/XVar_21(end));
XerrS1 = [XerrS1_3 XerrS1_6 XerrS1_10 XerrS1_15 XerrS1_21];
XerrS2_3 = abs(XVar2_3(end)/XVar_3(end));
XerrS2_6 = abs(XVar2_6(end)/XVar_6(end));
XerrS2_10 = abs(XVar2_10(end)/XVar_10(end));
XerrS2_15 = abs(XVar2_15(end)/XVar_15(end));
XerrS2_21 = abs(XVar2_21(end)/XVar_21(end));
XerrS2 = [XerrS2_3 XerrS2_6 XerrS2_10 XerrS2_15 XerrS2_21];

save('Xerr','Xerr'); save('err','err'); save('errS1','errS1'); save('errS2','errS2');

load('Xerr'); load('err'); load('errS1'); load('errS2');
P = [2 5 9 14 20]; 
figure, hold on;
plot(P,log(Xerr),'b-*','Linewidth',2);
plot(P,log(err),'r-*','Linewidth',2);
hold off; xlabel('P'); setleg(legend('PC SDE','PC PDE')); 

figure, hold on;
plot(P,log(errS1),'b-*','Linewidth',2);
plot(P,log(errS2),'r-*','Linewidth',2);
hold off; xlabel('P'); setleg(legend('S1','S2'));