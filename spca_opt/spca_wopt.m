function [P, P0, output] = spca_wopt(Z,W,R,varargin)
%SPCA_WOPT solves the following optimization problem:
%     min_{A,\Sigma,B} || W*(Z - A*diag(\Sigma)*B')|| +\beta*|\Sigma|
%
% using gradient-based optimization algorithms, where W is an indicator
% matrix containing zeros whenever data is missing (otherwise ones) and 
% R is the number of components. 
%
%   P = SPCA_WOPT(Z, W, R) 
% 
%   P = SPCA_WOPT(Z, W, R,'param', value,...) specifies additional
%   parameters for the method. Specifically...
%
%   'alg' - Specfies optimization algorithm (default: 'ncg')
%      'ncg'   Nonlinear Conjugate Gradient Method
%      'lbfgs' Limited-Memory BFGS Method
%      'tn'    Truncated Newton
%
%   'init' - Initialization for factor matrices. (default:'random'). 
%   This can be a cell array with the initial matrices or one of the 
%   following strings:
%      'random' Randomly generated via randn function
%      'nvecs'  Selected as leading left singular vectors of X(n)
%
%   'alg_options' - Parameter settings for selected optimization
%   algorithm. For example, type OPTIONS = NCG('defaults') to get
%   the NCG algorithm options which can then me modified as passed
%   through this function to NCG.
%
%    'beta': sparsity penalty parameter on the weights of rank-one components.
%   
%   See also SPCA_WFUN, SPCA_WFG, SPCA_VEC_TO_FAC, SPCA_FAC_TO_VEC

%% Check for POBLANO
if ~exist('poblano_params','file')
    error(['SPCA_WOPT requires Poblano Toolbox for Matlab. This can be ' ...
           'downloaded at http://software.sandia.gov/trac/poblano.']);
end

%% Set parameters
params = inputParser;
params.addParamValue('alg','ncg', @(x) ismember(x,{'ncg','tn','lbfgs'}));
params.addParamValue('init','random', @(x) (iscell(x) || ismember(x,{'random','nvecs'})));
params.addParamValue('alg_options', '', @isstruct);
params.addParamValue('beta', 0.001, @(x) x >= 0);
params.parse(varargin{:});

%% Set up optimization algorithm
switch (params.Results.alg)
    case 'ncg'
        opthandle = @ncg;
    case 'tn'
        opthandle = @tn;
    case 'lbfgs'
        opthandle = @lbfgs;
end

%% Set up optimization algorithm options
if isempty(params.Results.alg_options)
    options = feval(opthandle, 'defaults');
else
    options = params.Results.alg_options;
end

%% Set up function handle
if ~isa(Z,'tensor')
    Z=tensor(Z);
end

normZsqr = norm(Z)^2;
funhandle = @(x) spca_wfun(Z,W,x,R,normZsqr, params.Results.beta);
    
%% Initial guess
sz = size(Z);
N = length(sz);

if iscell(params.Results.init)
    P0 = params.Results.init;
elseif strcmpi(params.Results.init,'random')
    P0 = cell(N+1,1);
    for n=1:N
        P0{n} = randn(sz(n),R);
        for j=1:R
            P0{n}(:,j) = P0{n}(:,j) / norm(P0{n}(:,j));
        end
    end
    P0{N+1} = ones(R,1);
elseif strcmpi(params.Results.init,'nvecs')
    P0 = cell(N+1,1);
    for n=1:N
        P0{n} = nvecs(Z,n,R);
    end
    [~,S,~]=svds(Z.data,R);
    P0{N+1} = diag(S);    
else
    error('Initialization type not supported')
end

%% Fit SPCA_WOPT
out = feval(opthandle, funhandle, spca_fac_to_vec(P0), options);

T  = spca_vec_to_fac(out.X,Z,R);
P  = ktensor(T{end}, T(1:end-1)); 
if nargout > 1
    output.ExitFlag  = out.ExitFlag;
    output.FuncEvals = out.FuncEvals;
    output.f = out.F;
    output.G = spca_vec_to_fac(out.G,Z,R);
    output.OptOut = out;
end


