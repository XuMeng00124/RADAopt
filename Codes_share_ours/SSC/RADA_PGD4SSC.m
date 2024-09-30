function [U_best, out]= RADA_PGD4SSC(U, opts, L, Y_now)
%% The RADA-PGD algorithm for sparse SC problem
%
% min_{UUt} max_Y - < L, UUt > + <Y,UUt>
% s.t.  UUt=U*U'; U'*U = I; ||Y||_\infinity<=mu.
%
% Input:
%       X      --- data matrices
%       (U,Y)  --- initial guess of the variables
%
%       opts --- option structure with fields:
%            epsilon: stopping tolerance
%            mxitr: maximum number of iteration
%            R: max_Y \|Y\|, s.t. \|Y\|_\infinity<=mu
%            lambda: regularization paramter
%            nupdateU: number of Riemannian gradient descent steps at each
%            iteration
%            beta0: regularization paramter
%            rho: decreasing factor for parameter beta
%
% Output:
%       U_best   ---  solution
%       out      ---  output information
%
% Reference:
% Meng Xu, Bo Jiang, Ya-Feng Liu, and Anthony Man-Cho So, A Riemannian
% Alternating Descent Ascent Algorithmic Framework for Nonconvex-Linear
% Minimax Problems on Riemannian Manifolds, 2024
%
%
% Author: Meng Xu, Bo Jiang, Ya-Feng Liu, and Anthony Man-Cho So
%   Version 0.0 ... 2024/10
%
%  Contact information: xumeng22@mails.ucas.ac.cn, jiangbo@njnu.edu.cn, yafliu@lsec.cc.ac.cn
%-------------------------------------------------------------------------

if isempty(U)
    error('input U is an empty matrix');
else
    [d, r] = size(U);
    U_best = zeros(d,r);
end

if ~isfield(opts, 'epsilon');      opts.epsilon = 1e-8; end
if ~isfield(opts, 'adbeta');      opts.adbeta = 0; end
if ~isfield(opts, 'mxitr');     opts.mxitr  = 10000; end
if ~isfield(opts, 'R');      opts.R = d; end
if ~isfield(opts, 'rho');      opts.rho = 1.5; end
if ~isfield(opts, 'lambda');      opts.lambda = 1e-8; end
if ~isfield(opts, 'beta0');      opts.beta0 = 1; end

%-------------------------------------------------------------------------------

sparse = opts.lasso_constant;
lambda = opts.lambda;
epsilon = opts.epsilon;
rho = opts.rho;
beta0 = opts.beta0;
n = opts.n;
beta0 = beta0*n*sqrt(r);
beta = beta0;

F_hist = zeros(opts.mxitr,1);
Phi_hist = zeros(opts.mxitr,1);
nls_hist = zeros(opts.mxitr,1);
beta_hist = zeros(opts.mxitr,1);
residual = 1e10;
F_histbest = inf;

G = U;

grad = (L+Y_now) * U;
GU = grad'*U;
dtU = grad - U*GU;     nrmG  = norm(dtU, 'fro');

%% Main Iteration

for itr = 1 : opts.mxitr

    npU=1;
    for nupdateU = 1:npU

        U = G;
        UUt = U*U';
        f = sum(L.*UUt,"all") + sparse * sum(sum(abs(UUt)));
        Y = max(-sparse,min((1/(lambda+beta)*UUt + beta/(lambda+beta)*Y_now),sparse));
        F = sum((L+Y).*UUt,"all") - lambda/2*norm(Y,'fro')^2 - beta/2*norm(Y-Y_now,'fro')^2;

        Phi_hist(itr) = F;
        F_hist(itr) = f;
        if itr>=2 && abs(F_hist(itr)-F_hist(itr-1))<1e-8
            breakindex = breakindex+1;
        else
            breakindex = 0;
        end
        nls_hist(itr) = 1; %nls;
        if F_hist(itr) < F_histbest

            U_best = U;
            F_histbest = F_hist(itr);
        end

        if nupdateU == npU

            if opts.adbeta

                tau1 = 0.999;
                residualnext = max(max(abs(lambda*Y+beta*(Y_now - Y))));
                if (residualnext >= 1e-3 && residualnext >= tau1 * residual)
                    beta0 = 0.9*beta0;
                end

                beta = max(beta0/((itr)^rho),1e-9);
                residual = residualnext;
            else
                beta = beta0/((itr+1)^rho);
            end
            Y_now = Y;
            beta_hist(itr) = beta;
        end

        Y = max(-sparse,min((1/(lambda+beta)*UUt + beta/(lambda+beta)*Y_now),sparse));
        [G,~] = eigs(-(lambda+beta) * (L + Y) + UUt, r, 'largestreal');

        G_temp = L+Y_now;
        GUUt_temp = G_temp*U*U';
        dtU_temp = GUUt_temp + GUUt_temp'-2*U*(U'*GUUt_temp);
        nrmG = norm(dtU_temp,'fro');

    end

    if nrmG <= epsilon && norm(Y - max(-sparse,min((UUt + Y),sparse)),'fro') <= epsilon
        break;
    end

end

out.nrmG = nrmG;
out.fval = F;
out.itr = itr;
out.obj = F_hist(1:itr);
out.Phi_hist = Phi_hist(1:itr);
out.beta_hist = beta_hist(1:itr);
out.F_best = F_histbest;

end