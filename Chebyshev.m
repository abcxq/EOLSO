function X=Chebyshev(N,dim,ub,lb)
X(1, :) = rands(1, dim);
for j = 1:dim
    for i = 1:N-1
        X(i+1, j) = cos(i*acos(X(i, j)));
    end
end
X = 0.5*(X+1).*(ub-lb)+lb;
