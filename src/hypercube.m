function rep = hypercube(N,k)
%Function using the latin hypercube method
%N: Number of points to generate
%k: Dimension of the sample

%Generation of the random numbers matrix
Matrix = rand(N,k);

%Loop to generate the permutations and mixing the terms
for j=1:k
    per = randperm(N);
    Matrix(:,j) = (per'-1+Matrix(:,j))/N;
end
rep = Matrix;

end
