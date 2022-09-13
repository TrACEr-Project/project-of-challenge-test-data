function A = spca_vec_to_fac(x, Z, R)
%SPCA_VEC_TO_FAC Converts a vector to a cell array of factor matrices.
%
%   A = SPCA_VEC_TO_FAC(X,Z) converts the vector X into a cell array
%   of factor matrices consistent with the size of the tensor Z.
%
%   See also SPCA_FAC_TO_VEC, SPCA_FUN, SPCA_OPT, SPCA_FG.
%

%% Set-up
N = ndims(Z);
sz = size(Z);

%% Create A
A = cell(N+1,1);
for n = 1:N
    idx1 = sum(sz(1:n-1))*R + 1;
    idx2 = sum(sz(1:n))*R;
    A{n} = reshape(x(idx1:idx2),sz(n),R);
end
A{N+1} = x(idx2+1:end);