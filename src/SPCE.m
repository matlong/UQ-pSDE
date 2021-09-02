function X = SPCE(order)

NbTraj = 10000;  mu = [1,.1];  sigma = [.05,.05];  X0 = .99;  
NbStep = 1000;  Tfin = 10;  dt = Tfin/NbStep;  T = [0:dt:Tfin];    
syms x y;
psi = Legendre2D(x,y,order);

A1 = mu(1)*IP_Legendre2D(psi(1:end-1)',psi(1:end-1),x,y) + ...
     sigma(1)*IP_Legendre2D(psi(2)*psi(1:end-1)',psi(1:end-1),x,y); 
B1 = mu(2)*IP_Legendre2D(psi(1),psi(1:end-1)',x,y) + ...
     sigma(2)*IP_Legendre2D(psi(3),psi(1:end-1)',x,y);
 
X = TrajsPCE(NbTraj,order,NbStep,Tfin,X0,A1,B1);

end
