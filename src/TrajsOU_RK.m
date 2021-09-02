function X = toto_OU(NbTraj,NbStep,Tfin,X0,xi,mu,sigma,omega)

DeltaT = Tfin/NbStep;  
alpha = mu(1)+sqrt(3)*sigma(1)*(2*xi(:,1)-1); 
beta = mu(2)+sqrt(3)*sigma(2)*(2*xi(:,2)-1);   

X = zeros(NbTraj,NbStep+1);
X(:,1) = X0;    
for i = 1 : NbStep
      X(:,i+1) = X(:,i).*exp(-alpha*DeltaT) + beta.*sqrt((1-exp ...
                 (-2*alpha*DeltaT))./(2*alpha)).*omega(:,i); 
end

end


