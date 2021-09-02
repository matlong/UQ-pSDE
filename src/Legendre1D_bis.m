function r = Legendre1D_bis(order)

psi = cell(1,order);
psi{1} = @(x) 1;
if (order > 1)
    psi{2} = @(x) 2*x-1;
    for n = 2 : order-1
        psi{n+1} = @(x) (1/n)*((2*(n-1)+1)*(2*x-1).*psi{n}(x) - ...
                           (n-1)*psi{n-1}(x));
    end
end
for k = 1 : order
    r{k} = @(x) sqrt(2*(k-1)+1)*psi{k}(x); 
end
end