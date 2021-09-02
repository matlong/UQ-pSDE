function [uh,A,b,K,M] = oned_linear_FEM(ne,alpha,beta)

% 1D FE mesh
h = 1/ne; xx = [0:h:1]'; nvtx = length(xx);
Q = ne -1; elt2vert = [1:(Q+1); 2:(Q+2)]';
p = (beta^2/2)*ones(ne,1); f = ones(ne,1);
x = [(h/2):h:(1-h/2)]'; q = alpha*x;

% global matrices & vectors
K = sparse(nvtx,nvtx);
M = sparse(nvtx,nvtx);
b = zeros(nvtx,1);
[Kks,Mks,bks] = get_elt_arrays(h,p,q,f,ne);
for i = 1 : 2
    nrow = elt2vert(:,i);
    for j = 1 : 2
        ncol = elt2vert(:,j);
        K = K + sparse(nrow,ncol,Kks(:,i,j),nvtx,nvtx);
        M = M + sparse(nrow,ncol,Mks(:,i,j),nvtx,nvtx);
    end
    b = b + sparse(nrow,1,bks(:,i),nvtx,1);
end

% boundary conditions
K([1,end],:) = []; K(:,[1,end]) = [];
M([1,end],:) = []; M(:,[1,end]) = [];
b(1) = []; b(end) = [];

% solve linear system
A = K + M;
u_tmp = A\b;
uh = [0;u_tmp;0];

%plot(xx,uh,'-k');

end