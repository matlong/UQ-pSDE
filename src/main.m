% Estimate Sobol' indices for parameterized Ornstein-Uhlenbeck (OU)
% process using different methods: 
% 1) Double Monte Carlo method (MC);
% 2) Polynomial Chaos expansion (PC);
% 3) First transform the parameterized stochastic differential equation 
%    (SDE) into a parametrized (deterministic) partial differential 
%    equation (PDE) using the Feynman-Kac formula, then apply either the MC
%    or the PC method.
%
% Reference: 
% Sobol' sensitivity analysis for a parametrized diffusion 
% process, MSc Thesis, L.Li, 2017.
%
% Written by L. Li - May 2017.
%

close all; clear all; clc;

%% Set parameters
%
% Parametrized SDE (OU process):
%      dX_t = - alpha X_t dt + beta dW_t,
% where the drift coef. 'alpha' and the diffusion coef. 'beta' are two
% independent uniform random variables, and 'W' is a Wiener process.
% 
param = [];

param.nTraj = 10000; % Number of trajectories for 'W'      
param.nStep = 1000; % Number of time steps for SDE & PDE       
param.nCell = 100; % Number of spatial cells 
param.order_PC = 10; % Order of PC
param.Tend = 10; % SDE simulation duration
param.mu = [1,.1]; % mean of 'alpha' and 'beta'
param.sigma = [.05,.05]; % std of alpha/sqrt(3) and beta//sqrt(3)
param.X0 = .99; % Initial condition of X

% Derived param.
param.dt = param.Tend/param.nStep; % Time step for SDE & PDE 
param.t = 0:param.dt:param.Tend; % Time axis
param.xm = linspace(0,1,param.nCell+1); % Spatial mesh  
param.dx = param.xm(2)-param.xm(1); % Grid spacing

% Define normalized bivariate Legendre polynomials 
param.PC = Legendre2D(param.order_PC);

%% Double-MC method 

n = 10;
xi_A = [rand(NbTraj,1),rand(NbTraj,1)];
xi_B = [rand(NbTraj,1),rand(NbTraj,1)];
xi_C1 = [xi_B(:,1),xi_A(:,2)];
xi_C2 = [xi_A(:,1),xi_B(:,2)];
YB = titi_OU(NbTraj,NbStep,Tfin,X0,mu,sigma,xi_B,n);
YC1 = titi_OU(NbTraj,NbStep,Tfin,X0,mu,sigma,xi_C1,n);
YC2 = titi_OU(NbTraj,NbStep,Tfin,X0,mu,sigma,xi_C2,n);
var1 = mean(YB.*YC1,1) - mean(.5*(YB+YC1),1).^2;
var2 = mean(YB.*YC2,1) - mean(.5*(YB+YC2),1).^2;
Var1 = mean(.5*(YB.^2+YC1.^2),1) - mean(.5*(YB+YC1),1).^2;
Var2 = mean(.5*(YB.^2+YC2.^2),1) - mean(.5*(YB+YC2),1).^2;
figure, hold on;
plot(T,var1./Var1,'Linewidth',2);
plot(T,var2./Var2,'Linewidth',2);
hold off; legend('xi1','xi2');

Xi_A = [rand(NbTraj,1),rand(NbTraj,1),randn(NbTraj,NbStep)]; % sampling1
Xi_B = [rand(NbTraj,1),rand(NbTraj,1),randn(NbTraj,NbStep)]; % sampling2
%A = [LHS(NbTraj,1),LHS(NbTraj,1),randn(NbTraj,NbStep)];
%B = [LHS(NbTraj,1),LHS(NbTraj,1),randn(NbTraj,NbStep)];
Xi_C1 = [Xi_B(:,1),Xi_A(:,2),Xi_A(:,3:end)]; % for xi_1
Xi_C2 = [Xi_A(:,1),Xi_B(:,2),Xi_A(:,3:end)]; % for xi_2plot(T,Var12_MC,'Linewidth',2);
Xi_C3 = [Xi_A(:,1),Xi_A(:,2),Xi_B(:,3:end)]; % for z 
Xi_C12 = [Xi_B(:,1),Xi_B(:,2),Xi_A(:,3:end)];

tic;  
XA = TrajsOU(NbTraj,NbStep,Tfin,X0,Xi_A,mu,sigma); % output using parameter A
XB = TrajsOU(NbTraj,NbStep,Tfin,X0,Xi_B,mu,sigma); % using B
XC1 = TrajsOU(NbTraj,NbStep,Tfin,X0,Xi_C1,mu,sigma); % using C1
XC2 = TrajsOU(NbTraj,NbStep,Tfin,X0,Xi_C2,mu,sigma); % using C2
XC3 = TrajsOU(NbTraj,NbStep,Tfin,X0,Xi_C3,mu,sigma); % using C3
XC12 = TrajsOU(NbTraj,NbStep,Tfin,X0,Xi_C12,mu,sigma);

% numerator of SobolEff() -> estimation of partial variances 
Var1_MC = mean(XB.*XC1,1) - mean(.5*(XB+XC1),1).^2; %save('V1MC','Var1_MC');
Var2_MC = mean(XB.*XC2,1) - mean(.5*(XB+XC2),1).^2; %save('V2MC','Var2_MC');
Var3_MC = mean(XB.*XC3,1) - mean(.5*(XB+XC3),1).^2; %save('V3MC','Var3_MC');
Var12_MC = mean(XB.*XC12,1) - mean(.5*(XB+XC12),1).^2; %save('V12MC','Var12_MC');
Var_MC = mean(.5*(XB.^2+XC12.^2),1) - mean(.5*(XB+XC12),1).^2; %save('VMC','Var_MC');
t_MC = toc;
disp(['Execution time by Sobol based method = ' num2str(t_MC) 's']);


%% PC_SDE method

A1 = mu(1)*IP_Legendre2D(psi(1:end-1)',psi(1:end-1),x,y) + ...
     sigma(1)*IP_Legendre2D(psi(2)*psi(1:end-1)',psi(1:end-1),x,y); 
B1 = mu(2)*IP_Legendre2D(psi(1),psi(1:end-1)',x,y) + ...
     sigma(2)*IP_Legendre2D(psi(3),psi(1:end-1)',x,y);

tic;
% sol of NbTraj-systems of order-coupled EDS (cf. TrajsPCE.m)
X = TrajsPCE(NbTraj,order,NbStep,Tfin,X0,A1,B1); %save('XPCE','X');

% estmation of partial variances
Var1_PCE = (mean(X(2,:,:),3)).^2 + (mean(X(4,:,:),3)).^2 + ... 
           (mean(X(7,:,:),3)).^2; %save('V1PCE','Var1_PCE');
Var2_PCE = (mean(X(3,:,:),3)).^2 + (mean(X(6,:,:),3)).^2 + ...
           (mean(X(10,:,:),3)).^2; %save('V2PCE','Var2_PCE');
Var3_PCE = mean((X(1,:,:)).^2,3) - (mean(X(1,:,:),3)).^2; %save('V3PCE','Var3_PCE');
Var12_PCE = sum(mean(X(2:end,:,:),3).^2,1); %save('V12PCE','Var12_PCE');
Var123_PCE = sum(mean(X(2:end,:,:).^2,3)-mean(X(2:end,:,:),3).^2,1); %save('V123PCE','Var123_PCE');
Var_PCE = Var12_PCE + Var3_PCE + Var123_PCE; %save('VPCE','Var_PCE');
t_PCE = toc;
disp(['Execution time by Hybrid scheme = ' num2str(t_PCE) 's']);

figure,
for k = 1 : order
    subplot(5,2,k)
    hold on;
    for l = 1 : 10     % the first 10 trajectories samples
        plot(T,X(k,:,l))
        title(sprintf('[X_%s]',num2str(k-1)))
    end
    hold off;
end


%% PC_PDE method

B2 = .5*mu(2)^2*IP_Legendre2D(psi(1:end-1)',psi(1:end-1),x,y) + ...
  mu(2)*sigma(2)*IP_Legendre2D(psi(3)*psi(1:end-1)',psi(1:end-1),x,y) + ...
 .5*sigma(2)^2*IP_Legendre2D(psi(3)*psi(3)*psi(1:end-1)',psi(1:end-1),x,y);  
figure, imagesc(B2); colorbar;

A_bar = gallery('tridiag',NbCell-1,-1,0,1); % 1st derivative
B_bar = gallery('tridiag',NbCell-1,1,-2,1); % 2nd derivative   
D = diag(xx(2:end-1));

LHS = eye((NbCell-1)*order)/dt + kron(A1,D*A_bar)/(4*dx) ...
      - kron(B2,B_bar)/(2*dx^2);
RHS = eye((NbCell-1)*order)/dt - kron(A1,D*A_bar)/(4*dx) ...
      + kron(B2,B_bar)/(2*dx^2);
  
% ------------------------ boundary condition U ------------------------- %
psi = Legendre2D_bis(order);
alpha = @(x,y) mu(1)+sigma(1)*psi{2}(x,y);
beta = @(x,y) mu(2)+sigma(2)*psi{3}(x,y);

coefUM = zeros(order,NbStep+1);
tic;
for i = 1 : NbStep+1
    tmp1 = @(x,y) exp(-alpha(x,y)*T(i));
    for k = 1 : order
        tmp2 = @(x,y) psi{k}(x,y);
        fun = @(x,y) tmp1(x,y).*tmp2(x,y);
        coefUM(k,i) = integral2(fun,0,1,0,1);
    end
end
toc;
fM = (-A1*xx(end-1)/(4*dx) + B2/(2*dx^2))*coefUM;

tic;
fM_bis = zeros(order,NbStep+1);
tmp1 = @(x,y) (-xx(end-1)/(4*dx))*alpha(x,y)+(beta(x,y).^2)/(4*dx^2);
for i = 1 : NbStep+1
    tmp2 = @(x,y) exp(-alpha(x,y)*T(i));
    for k = 1 : order
        tmp3 = @(x,y) psi{k}(x,y);
        fun = @(x,y) tmp1(x,y).*tmp2(x,y).*tmp3(x,y);
        fM_bis(k,i) = integral2(fun,0,1,0,1);
    end
end
toc;
figure, imagesc(fM-fM_bis); colorbar; xlabel('t_i'); ylabel('Psi_k');

F = zeros(order*(NbCell-1),NbStep+1);
index = [1:1:order];
F(index*(NbCell-1),:) = fM(index,:);
% ------------------------------ fin ------------------------------------ %

U = zeros((NbCell-1)*order,NbStep+1);
U0 = zeros((NbCell-1)*order,1);
U0(1:NbCell-1) = xx(2:end-1);
U(:,1) = U0;

tic;
for i = 1 : NbStep
    U(:,i+1) = LHS\(RHS*U(:,i)+F(:,i+1)+F(:,i));
end
toc;
%save('UPDE','U');

Var12_PDE = zeros(NbCell-1,NbStep+1);   % xi
Var1_PDE = zeros(NbCell-1,NbStep+1);  % xi_1
Var2_PDE = zeros(NbCell-1,NbStep+1);  % xi_2
for i = 1 : NbStep+1
    Utmp = reshape(U(:,i),[NbCell-1 order]);
    Var12_PDE(:,i) = sum(Utmp(:,2:end).^2,2);
    Var1_PDE(:,i) = Utmp(:,2).^2 + Utmp(:,4).^2 + Utmp(:,7).^2;
    Var2_PDE(:,i) = Utmp(:,3).^2 + Utmp(:,6).^2 + Utmp(:,10).^2;
end
%save('V1PDE','Var1_PDE'); save('V2PDE','Var2_PDE');
%save('V12PDE','Var12_PDE');

u = zeros(NbCell-1,order,NbStep+1);
for i = 1 : NbStep+1
    u(:,:,i) = reshape(U(:,i),[NbCell-1 order]);
end

figure,
for k = 1 : order
    ue = reshape(u(:,k,:),[NbCell-1 NbStep+1]);
    subplot(4,3,k)
    mesh(T,xx(2:end-1),ue); 
    colorbar; xlabel('t'); ylabel('x'); 
    title(sprintf('u_%s',num2str(k-1)))
end

% ------------------------ boundary condition V ------------------------- %
coefV0 = zeros(order,NbStep+1); coefVM_tmp = zeros(order,NbStep+1); 
tic;
for i = 1 : NbStep+1
    tmp1 = @(x,y) exp(-2*alpha(x,y)*T(i));
    tmp2 = @(x,y) (beta(x,y).^2).*((1-tmp1(x,y))./(2*alpha(x,y)));
    for k = 1 : order
        tmp3 = @(x,y) psi{k}(x,y);
        fun1 = @(x,y) tmp2(x,y).*tmp3(x,y);
        fun2 = @(x,y) tmp1(x,y).*tmp3(x,y);
        coefV0(k,i) = integral2(fun1,0,1,0,1);
        coefVM_tmp(k,i) = integral2(fun2,0,1,0,1);
    end
end
coefVM = coefV0 + coefVM_tmp;
toc;
g0 = (A1*xx(2)/(4*dx) + B2/(2*dx^2))*coefV0;
gM = (-A1*xx(end-1)/(4*dx) + B2/(2*dx^2))*coefVM;
G = zeros(order*(NbCell-1),NbStep+1);
G((index-1)*(NbCell-1)+1,:) = g0(index,:);
G(index*(NbCell-1),:) = gM(index,:);
% --------------------------------- fin --------------------------------- %

V = zeros((NbCell-1)*order,NbStep+1);
V0 = zeros((NbCell-1)*order,1);
V0(1:NbCell-1) = xx(2:end-1).^2;
V(:,1) = V0;

tic;
for i = 1 : NbStep
    V(:,i+1) = LHS\(RHS*V(:,i)+G(:,i+1)+G(:,i));
end
toc;

% total variance
Var_PDE = zeros(NbCell-1,NbStep+1);   
for i = 1 : NbStep+1
    Utmp = reshape(U(:,i),[NbCell-1 order]);
    Vtmp = reshape(V(:,i),[NbCell-1 order]);
    Var_PDE(:,i) = Vtmp(:,1) - Utmp(:,1).^2;
end
%save('VPDE','Var_PDE');


%% Exact solutions

w = sqrt(3)*sigma(1)*T;
E1 = X0*exp(-mu(1)*T).*sinh(w)./w;
E2 = X0^2*exp(-2*mu(1)*T).*sinh(2*w)./(2*w);
E3 = mu(2)^2+sigma(2)^2;
inf = mu(1)-sqrt(3)*sigma(1); sup = mu(1)+sqrt(3)*sigma(1);
E4 = log(sup/inf)/(4*sqrt(3)*sigma(1));
fun = @(x)exp(-2*x*T)/(2*x);
q = integral(fun,inf,sup,'ArrayValued',true);
E5 = q/(2*sqrt(3)*sigma(1));

Var12_EX = E2 - E1.^2;  % save('V12EX','Var12_EX');
Var_EX = E3.*(E4-E5) + Var12_EX;  % save('VEX','Var_EX');


%% Plots

% load('./data/V1MC'); load('./data/V2MC'); load('./data/V3MC'); 
% load('./data/V12MC'); load('./data/VMC');
% load('./data/V1PCE'); load('./data/V2PCE'); load('./data/V3PCE'); 
% load('./data/V12PCE'); load('./data/V123PCE'); load('./data/VPCE'); 
% load('./data/V12EX'); load('./data/VEX'); load('./data/V1PDE'); 
% load('./data/V2PDE'); load('./data/V12PDE'); load('./data/VPDE'); 

figure, hold on;
plot(T,Var_PCE,'Linewidth',2);
plot(T,Var12_PCE,'Linewidth',2);
plot(T,Var3_PCE,'Linewidth',2);
plot(T,Var123_PCE,'Linewidth',2);
hold off; ylim([0 .01]); xlabel('t'); 
setleg(legend('total','param','noise','mix'));

figure, hold on;
plot(T,Var12_EX,'Linewidth',2);
plot(T,Var1,'Linewidth',2);
plot(T,Var12_PCE,'Linewidth',2);
plot(T,Var12_PDE(end,:),'Linewidth',2);
hold off; xlabel('t'); setleg(legend('EX','MC','PCE','PDE'));

figure, hold on;
plot(T(2:end),Var12_PCE(2:end)./Var_PCE(2:end),'Linewidth',2);
plot(T(2:end),Var3_PCE(2:end)./Var_PCE(2:end),'Linewidth',2);
plot(T(2:end),Var123_PCE(2:end)./Var_PCE(2:end),'Linewidth',2);
hold off; ylim([0 1]); xlabel('t'); setleg(legend('param','noise','mix')); 

figure, % compare methods
subplot(2,2,1)
hold on;
plot(T,Var12_MC,'Linewidth',2);
plot(T,Var1_PCE,'Linewidth',2);
hold off; ylim([-.0001 .0005]);
xlabel('t'); title('Var ( E ( X_t | xi_1 ) )'); legend('MC','PCE');
subplot(2,2,2)
hold on; ylim([-.0002 .0005]);
plot(T,Var2_MC,'Linewidth',2);
plot(T,Var2_PCE,'Linewidth',2);
hold off;
xlabel('t'); title('Var ( E ( X_t | xi_2 ) )'); legend('MC','PCE');
subplot(2,2,3)
hold on;
plot(T,Var3_MC,'Linewidth',2);
plot(T,Var3_PCE,'Linewidth',2);
hold off; ylim([0 .007]);
xlabel('t'); title('Var ( E ( X_t | W ) )'); legend('MC','PCE');
subplot(2,2,4)
hold on;
plot(T,Var_MC,'Linewidth',2);
plot(T,Var_PCE,'Linewidth',2);
hold off;
xlabel('t'); title('Var ( X_t )'); legend('MC','PCE');

figure,
subplot(1,2,1)
mesh(T,xx(2:end-1),Var12_PDE); colorbar;
xlabel('t'); ylabel('x'); zlabel('Var(y1)'); title('Total variance evolution of a parametrized Ornstein-Uhlenbeck process');
subplot(1,2,2)
mesh(T,z(2:end-1),Var_PDE); colorbar;
xlabel('t'); ylabel('x'); title('Total variance');

figure, hold on;
plot(T,Var12_EX,'Linewidth',2);
plot(T,Var12_MC,'Linewidth',2);
plot(T,Var12_PCE,'Linewidth',2);
plot(T,Var12_PDE(end,:),'Linewidth',2);
hold off; xlabel('t'); setleg(legend('Exact','MC SDE','PC SDE','PC PDE'));
ylabel('Var(y1)');
%print('fig1','-depsc');
figure, hold on;
plot(T,Var_EX,'Linewidth',2);
plot(T,Var_MC,'Linewidth',2);
plot(T,Var_PCE,'Linewidth',2);
plot(T,Var_PDE(end,:),'Linewidth',2);
hold off; xlabel('t'); setleg(legend('EX','MC','PCE','PDE'));
%print('fig2','-depsc');
figure, hold on;
%plot(T,Var1_EX./Var_EX,'Linewidth',2);
%plot(T(2:end),Var1_MC(2:end)./Var12_MC(2:end),'Linewidth',2);
%plot(T,Var1_PCE./Var12_PCE,'Linewidth',2);
plot(T,Var1_PDE(end,:)./Var12_PDE(end,:),'Linewidth',2);
plot(T,Var2_PDE(end,:)./Var12_PDE(end,:),'Linewidth',2);
hold off; xlabel('t'); %setleg(legend('EX','MC','PCE','PDE'));
legend('xi1','xi2');
%print('fig3','-depsc');

err_MC = Var12_MC./Var_MC - Var12_EX./Var_EX;
err_PCE = Var12_PCE./Var_PCE - Var12_EX./Var_EX;
err_PDE = Var12_PDE(end,:)./Var_PDE(end,:) - Var12_EX./Var_EX;
figure, hold on;
plot(T,err_PCE,'Linewidth',2);
plot(T,err_PDE,'Linewidth',2);
hold off; xlabel('t'); setleg(legend('PCE','PDE')); print('fig6','-depsc');


%% elliptic 

alpha = mu(1)+sqrt(3)*sigma(1)*(2*rand-1); 
beta = mu(2)+sqrt(3)*sigma(2)*(2*rand-1); 
sol = oned_linear_FEM(NbCell,alpha,beta);
plot(xx,sol,'Linewidth',2);

% MC_FEM
xi_A = rand(NbTraj,2);  xi_B = rand(NbTraj,2); 
xi_C1 = [xi_B(:,1),xi_A(:,2)];  xi_C2 = [xi_A(:,1),xi_B(:,2)]; 

uA = PDE_MC(mu,sigma,xi_A,NbCell,NbTraj);
uB = PDE_MC(mu,sigma,xi_B,NbCell,NbTraj);
uC1 = PDE_MC(mu,sigma,xi_C1,NbCell,NbTraj);
uC2 = PDE_MC(mu,sigma,xi_C2,NbCell,NbTraj);

Var1_tauD = mean(uB.*uC1,1) - mean(.5*(uB+uC1),1).^2; 
Var2_tauD = mean(uB.*uC2,1) - mean(.5*(uB+uC2),1).^2; 
VarT1_tauD = mean(.5*(uB.^2+uC1.^2),1) - mean(.5*(uB+uC1),1).^2;
VarT2_tauD = mean(.5*(uB.^2+uC2.^2),1) - mean(.5*(uB+uC2),1).^2;

% stochastic Galerkin method
x = [(dx/2):dx:(1-dx/2)]';
[Kks,Mks,bks] = get_elt_arrays(dx,ones(NbCell,1),x,ones(NbCell,1),NbCell);

nvtx = length(xx);
Q = NbCell -1; elt2vert = [1:(Q+1); 2:(Q+2)]';

K = sparse(nvtx,nvtx);
M = sparse(nvtx,nvtx);
b0 = zeros(nvtx,1);
for i = 1 : 2
    nrow = elt2vert(:,i);
    for j = 1 : 2
        ncol = elt2vert(:,j);
        K = K + sparse(nrow,ncol,Kks(:,i,j),nvtx,nvtx);
        M = M + sparse(nrow,ncol,Mks(:,i,j),nvtx,nvtx);
    end
    b0 = b0 + sparse(nrow,1,bks(:,i),nvtx,1);
end
K([1,end],:) = []; K(:,[1,end]) = [];
M([1,end],:) = []; M(:,[1,end]) = [];
b0(1) = []; b0(end) = [];

L = kron(B2,K) + kron(A1,M);
%b = zeros((NbCell-1)*order,1); b(1:NbCell-1) = ones(NbCell-1,1); 
b = zeros((NbCell-1)*order,1); b(1:NbCell-1) = b0; 
u_tmp = L\b;
uu = reshape(u_tmp,NbCell-1,order);
uh = [zeros(1,order);uu;zeros(1,order)];

toto = (uh(:,2:end).^2)*ones(order-1,1);
toto1 = uh(:,2).^2 + uh(:,4).^2 + uh(:,7).^2;
toto2 = uh(:,3).^2 + uh(:,6).^2 + uh(:,10).^2;

figure, hold on;
plot(xx,toto1./toto,'b','Linewidth',2);
plot(xx,Var1_tauD./VarT1_tauD,'b--','Linewidth',2);
plot(xx,toto2./toto,'r','Linewidth',2);
plot(xx,Var2_tauD./VarT2_tauD,'r--','Linewidth',2);
plot(xx,(toto-toto1-toto2)./toto,'y','Linewidth',2);
plot(xx,1-(Var1_tauD./VarT1_tauD+Var2_tauD./VarT2_tauD),'y--','Linewidth',2);
hold off; xlabel('x'); 
setleg(legend('S1 PC','S1 MC','S2 PC','S2 MC','Smix PC','Smix MC'));
title('Sobol indices for the mean exit time of the parametrized OU process from a bounded domain');

