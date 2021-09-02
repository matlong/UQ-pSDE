%% Return values of 1D-Legendre polynomials %%

function Px = evalLegendre1D(x,order)
% x : vector
% degree = order+1 
Px = zeros(length(x),order);
Px(:,1) = 1;
Px(:,2) = 2*x-1;
for n = 1 : (order-1)
    Px(:,n+2) = ((2*n+1).*(2*x-1).*Px(:,n+1)-n.* Px(:,n))/(n+1);
end
for k = 1 : order
    Px(:,k) = sqrt(2*(k-1)+1)*Px(:,k);
end
end
