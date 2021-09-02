function pc = Legendre1D(x,order)
%
%
%
psi(1) = sym(1);
if order > 1
    psi(2) = 2*x-1;
    for n = 2:order-1
        psi(n+1) = simplify(expand( ...
            (1/(n))*((2*(n-1)+1)*(2*x-1)*psi(n)-(n-1)*psi(n-1)) ...
            ));
    end
end
for k = 1:order
    pc(k) = sqrt(2*(k-1)+1)*psi(k); 
end

end