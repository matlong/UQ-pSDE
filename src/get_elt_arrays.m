function [Kks,Mks,bks] = get_elt_arrays(h,p,q,f,ne)

Kks = zeros(ne,2,2);
Kks(:,1,1) = p./h;
Kks(:,1,2) = -p./h;
Kks(:,2,1) = -p./h;
Kks(:,2,2) = p./h;

Mks = zeros(ne,2,2);
Mks(:,1,1) = -q./2;
Mks(:,1,2) = -q;
Mks(:,2,1) = q./2;
Mks(:,2,2) = q;

bks = zeros(ne,2);
bks(:,1) = f.*h./2;
bks(:,2) = f.*h./2;

end