function psi = Legendre2D_bis(order)

order1D = floor((1/2)*(1+sqrt(8*order-7)));

px = cell(1,order1D); py = cell(1,order1D);
px{1} = @(x) 1; py{1} = @(y) 1;
if (order1D > 1)
    px{2} = @(x) 2*x-1; py{2} = @(y) 2*y-1;
    for n = 2 : order1D-1
        px{n+1} = @(x) (1/n)*((2*(n-1)+1)*(2*x-1).*px{n}(x) - ...
                        (n-1)*px{n-1}(x));
        py{n+1} = @(y) (1/n)*((2*(n-1)+1)*(2*y-1).*py{n}(y) - ...
                        (n-1)*py{n-1}(y));
    end
end

for i = 1 : order1D
    psi_x{i} = @(x) sqrt(2*(i-1)+1)*px{i}(x);
    for j = 1 : order1D
        psi_y{j} = @(y) sqrt(2*(j-1)+1)*py{j}(y);
        k = ((i+j-1)*(i+j))/2-i+1;
        fun{k} = @(x,y) psi_x{i}(x).*psi_y{j}(y);
    end
end

for k = 1 : order+1
    psi{k} = @(x,y) fun{k}(x,y);
end

end