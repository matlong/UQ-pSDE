%% Hybrid scheme (with only one loop for NbStep) %%

function X = TrajsPCE(NbTraj,order,NbStep,Tfin,X0,A,B)

DeltaT = Tfin/NbStep;
DeltaW = sqrt(DeltaT)*randn(NbStep,NbTraj);

X = zeros(order,NbStep+1,NbTraj);
X(1,1,:) = X0;    % only [X_0]=X0
Xtmp = zeros(order,NbTraj);  % reduce dimension
for i = 1 : NbStep
    Xtmp(:,:) = X(:,i,:);  % for the instant i
    X(:,i+1,:) = Xtmp - A*Xtmp*DeltaT + B*DeltaW(i,:);
end
end