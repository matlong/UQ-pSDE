function pc = Legendre2D(order2)
%
%
%
order1 = floor((1/2)*(1+sqrt(8*order2-7)));

syms x y;
psi_x = Legendre1D(x,order1);
psi_y = Legendre1D(y,order1);

for i = 1 : order1
    for j = 1 : order1
        k = ((i+j-1)*(i+j))/2-i+1;
        psi_tmp(k) = simplify(expand(psi_x(i)*psi_y(j)));
    end
end

psi = psi_tmp(1:order2+1);

pc = psi;

end