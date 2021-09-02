function u = PDE_MC(mu,sigma,xi,ne,NbTraj)

alpha = mu(1)+sqrt(3)*sigma(1)*(2*xi(:,1)-1); 
beta = mu(2)+sqrt(3)*sigma(2)*(2*xi(:,2)-1);

u = zeros(NbTraj,ne+1);
for i = 1 : NbTraj
    u(i,:) = oned_linear_FEM(ne,alpha(i),beta(i));
end

end