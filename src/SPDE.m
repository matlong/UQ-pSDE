function u = SPDE(order)

NbTraj = 10000;  mu = [1,.1];  sigma = [.05,.05];  X0 = .99;  
NbStep = 1000;  Tfin = 10;  dt = Tfin/NbStep;  T = [0:dt:Tfin];       
NbCell = 100;  xx = linspace(0,1,NbCell+1);  dx = xx(2)-xx(1);

syms x y;
psi = Legendre2D(x,y,order);

A1 = mu(1)*IP_Legendre2D(psi(1:end-1)',psi(1:end-1),x,y) + ...
     sigma(1)*IP_Legendre2D(psi(2)*psi(1:end-1)',psi(1:end-1),x,y); 
B2 = .5*mu(2)^2*IP_Legendre2D(psi(1:end-1)',psi(1:end-1),x,y) + ...
  mu(2)*sigma(2)*IP_Legendre2D(psi(3)*psi(1:end-1)',psi(1:end-1),x,y) + ...
 .5*sigma(2)^2*IP_Legendre2D(psi(3)*psi(3)*psi(1:end-1)',psi(1:end-1),x,y);  
A_bar = gallery('tridiag',NbCell-1,-1,0,1); 
B_bar = gallery('tridiag',NbCell-1,1,-2,1); 
D = diag(xx(2:end-1));
LHS = eye((NbCell-1)*order)/dt + kron(A1,D*A_bar)/(4*dx) ...
      - kron(B2,B_bar)/(2*dx^2);
RHS = eye((NbCell-1)*order)/dt - kron(A1,D*A_bar)/(4*dx) ...
      + kron(B2,B_bar)/(2*dx^2);
  
P = Legendre2D_bis(order);
alpha = @(x,y) mu(1)+sigma(1)*P{2}(x,y);
beta = @(x,y) mu(2)+sigma(2)*P{3}(x,y);
coefUM = zeros(order,NbStep+1);
for i = 1 : NbStep+1
    tmp1 = @(x,y) exp(-alpha(x,y)*T(i));
    for k = 1 : order
        tmp2 = @(x,y) P{k}(x,y);
        fun = @(x,y) tmp1(x,y).*tmp2(x,y);
        coefUM(k,i) = integral2(fun,0,1,0,1);
    end
end
fM = (-A1*xx(end-1)/(4*dx) + B2/(2*dx^2))*coefUM;
F = zeros(order*(NbCell-1),NbStep+1);
index = [1:1:order];
F(index*(NbCell-1),:) = fM(index,:);

U = zeros((NbCell-1)*order,NbStep+1);
U0 = zeros((NbCell-1)*order,1);
U0(1:NbCell-1) = xx(2:end-1);
U(:,1) = U0;
for i = 1 : NbStep
    U(:,i+1) = LHS\(RHS*U(:,i)+F(:,i+1)+F(:,i));
end

u = U;

end
