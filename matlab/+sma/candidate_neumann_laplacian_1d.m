function [Dxx,G,audit]=candidate_neumann_laplacian_1d(n,dx)
%SMA.CANDIDATE_NEUMANN_LAPLACIAN_1D Symmetric cell-centred no-flux operator.
% Dxx is negative semidefinite and satisfies Dxx*ones(n,1)=0.  G stores
% interface differences, so -Dxx=G'*G for the uniform cell grid.

assert(isscalar(n) && n==floor(n) && n>=2 && isfinite(dx) && dx>0, ...
    'A grid with at least two cells and positive spacing is required.');
G=spdiags([-ones(n-1,1),ones(n-1,1)],[0,1],n-1,n)/dx;
Dxx=-(G.'*G);
one=ones(n,1);
symmetry=norm(Dxx-Dxx.','fro');
constant_error=norm(Dxx*one,inf);
row_sum_error=norm(sum(Dxx,2),inf);
lambda_max=max(eig(full((Dxx+Dxx.')/2)));
audit=struct('n',n,'dx',dx,'symmetry_error',symmetry, ...
    'constant_field_error',constant_error, ...
    'no_flux_row_sum_error',row_sum_error, ...
    'largest_eigenvalue',lambda_max,'negative_semidefinite',lambda_max<=1e-12, ...
    'production_connected',false);
end
