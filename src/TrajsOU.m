%% Trajectories of Ornstein-Uhlenbeck process %%

function X = TrajsOU(NbTraj,NbStep,Tfin,X0,param,mu,sigma)

% input('Choose method : 1-Euler scheme , 2-Transition density : ');
method = 2; 

DeltaT = Tfin/NbStep;  % StepSize
alpha = mu(1)+sqrt(3)*sigma(1)*(2*param(:,1)-1); % drift coef 
beta = mu(2)+sqrt(3)*sigma(2)*(2*param(:,2)-1); % diffusion coef 
Z = param(:,3:end);   
DeltaW = sqrt(DeltaT)*Z;  % increasement of Brownian motion

X = zeros(NbTraj,NbStep+1);
X(:,1) = X0;    % initial condition
for i = 1 : NbStep
    if (method == 1)  % Euler
      X(:,i+1) = X(:,i) - alpha.*X(:,i)*DeltaT + beta.*DeltaW(:,i);
    elseif (method == 2)  % Exact
      X(:,i+1) = X(:,i).*exp(-alpha*DeltaT) + ...
                 beta.*sqrt((1-exp(-2*alpha*DeltaT))./(2*alpha)).*Z(:,i); 
    end
end
end