function r = IP_Legendre1D(f,g,x)

h = collect(simplify(f*g),x);

r = double(int(h,x,0,1));