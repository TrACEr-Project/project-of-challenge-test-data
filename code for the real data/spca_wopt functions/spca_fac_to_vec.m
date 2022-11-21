function x = spca_fac_to_vec(A)
%SPCA_FAC_TO_VEC Converts a set of factor matrices to a vector
%
%  X = SPCA_FAC_TO_VEC(A) converts a cell array of factor matrices A to a vector
%  by vectorizing each matrix and stacking them.
%
%   See also SPCA_VEC_TO_FAC, SPCA_FUN, SPCA_OPT, SPCA_FG.
%

%% Set-up
N = length(A);

%% Get sizes
sz = zeros(N-1,1);
for n = 1:N-1
    sz(n) = size(A{n},1);
end
R = size(A{1},2);
P = sum(sz)*R;

%% Create x
x = zeros(P,1);
for n = 1:N-1
    idx1 = sum(sz(1:n-1))*R + 1;
    idx2 = sum(sz(1:n))*R;
    x(idx1:idx2) = reshape(A{n},sz(n)*R,1);    
end
x = [x;A{N}];
