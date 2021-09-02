function r = IP_Legendre2D(f,g,x,y)

h = expand(simplify(f*g));

r = double(int(int(h,x,0,1),y,0,1));
        