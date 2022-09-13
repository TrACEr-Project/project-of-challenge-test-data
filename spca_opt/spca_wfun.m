function [f,g] = spca_wfun(Zdata,W,x,R,normZsqr, beta)
%SPCA_WFUN Computes function and gradient for matrix factorization for the
%scalings are taken out
%
%   [F,G] = SPCA_WFUN(Z,x,R,normZsqr) calculates the function and gradient
%   for the function ||Z - A{1}SigmaA{2}'||_F^2, where R is the number of components 
%
%   See also SPCA_WOPT, SPCA_WFG, SPCA_VEC_TO_FAC, SPCA_FAC_TO_VEC
%

%% Convert x to factor matrices (i.e., a cell array).
A = spca_vec_to_fac(x,Zdata,R);

%% Compute the function and gradient
if isa(Zdata,'tensor')
    if ~exist('normZsqr','var')
        normZsqr = norm(Zdata)^2;
    end
    [f,G] = spca_wfg(Zdata,W, A,normZsqr, beta);
end

%% Convert gradient to a vector
g = spca_fac_to_vec(G);


