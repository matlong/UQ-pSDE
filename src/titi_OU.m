function Y = titi_OU(NbTraj,NbStep,Tfin,X0,mu,sigma,xi,n)

omega = randn(NbTraj,NbStep,n);
X = zeros(NbTraj,NbStep+1,n);
for i = 1 : n
    X(:,:,i) = toto_OU(NbTraj,NbStep,Tfin,X0,xi,mu,sigma,omega(:,:,i));
end

Y = mean(X,3);

end