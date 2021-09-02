function v = MC_PC_EF(mu,sigma)

NbTraj = 10000; NbStep = 1000; Tfin = 10; dt = Tfin/NbStep; T = [0:dt:Tfin];       
NbCell = 100;  xx = linspace(0,1,NbCell+1);  dx = xx(2)-xx(1);
syms x y;
order = 10; 
psi = Legendre2D(x,y,order);
A1 = mu(1)*IP_Legendre2D(psi(1:end-1)',psi(1:end-1),x,y) + ...
     sigma(1)*IP_Legendre2D(psi(2)*psi(1:end-1)',psi(1:end-1),x,y); 
B2 = .5*mu(2)^2*IP_Legendre2D(psi(1:end-1)',psi(1:end-1),x,y) + ...
  mu(2)*sigma(2)*IP_Legendre2D(psi(3)*psi(1:end-1)',psi(1:end-1),x,y) + ...
 .5*sigma(2)^2*IP_Legendre2D(psi(3)*psi(3)*psi(1:end-1)',psi(1:end-1),x,y);  

% MC
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

% PCE
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

end

