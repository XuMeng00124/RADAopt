function [U_best, out]= RADA_RGD4SSCGr(U, opts, L, Y)
%% The RADA-RGD algorithm for sparse SC problem
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
%            zeta: stepsize for Riemannian gradient descent
%            c1: parameter for line search
%            eta: parameter for line search
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
if ~isfield(opts, 'zeta');       opts.zeta  = 1e-3; end
if ~isfield(opts, 'c1');     opts.c1  = 1e-4; end
if ~isfield(opts, 'eta');       opts.eta  = 0.1; end
if ~isfield(opts, 'mxitr');     opts.mxitr  = 10000; end
if ~isfield(opts, 'R');      opts.R = d; end
if ~isfield(opts, 'nupdateU'); opts.nupdateU = 1; end

%-------------------------------------------------------------------------------
% copy parameters
c1 = opts.c1;
eta     = opts.eta;
sparseparameter = opts.lasso_constant;
npU = opts.nupdateU;
epsilon = opts.epsilon;
R = opts.R;
lambda = opts.lambda;
beta0 = opts.beta0;
n = opts.n;
beta0 = beta0*n*sqrt(r);
beta = beta0;
rho = opts.rho;

F_hist = zeros(opts.mxitr,1);
nls_hist = zeros(opts.mxitr,1);
beta_hist = zeros(opts.mxitr,1);
F_histbest = inf;

residual = 1e10;

Y_now = Y;
UUt = U*U';
Y = max(-sparseparameter,min((1/(lambda+beta)*UUt + beta/(lambda+beta)*Y_now),sparseparameter));
F = sum((L+Y).*UUt,"all") - lambda/2*norm(Y,'fro')^2 - beta/2*norm(Y-Y_now,'fro')^2;
G = L+Y;
GUUt = (G*U)*U';
dtU = GUUt + GUUt' - 2*U*(U'*GUUt);    nrmG  = norm(dtU, 'fro');

Q = 1; Cval = F;  zeta = opts.zeta;
%% Main Iteration

for itr = 1 : opts.mxitr

    for nupdateU = 1:npU

        UP = UUt;     % FP = F;   GP = G;
        dtUP = dtU;

        nls = 1; deriv = c1*nrmG^2;

        %% Line Search
        while 1

            [V,~] = eigs(UP - zeta*dtU, r, 'largestreal');
            U = V;

            UUt = U*U';

            Y = max(-sparseparameter,min((1/(lambda+beta)*UUt + beta/(lambda+beta)*Y_now),sparseparameter));

            F = sum((L+Y).*UUt,"all") - lambda/2*norm(Y,'fro')^2 - beta/2*norm(Y-Y_now,'fro')^2;


            if F <= Cval - zeta*deriv || nls >= 5

                f = sum(L.*UUt,"all") + sparseparameter * sum(sum(abs(UUt)));

                F_hist(itr) = f;
                nls_hist(itr) = nls;

                if itr>=2 && F_histbest-F_hist(itr) < 1e-8
                    breakindex = breakindex+1;
                else
                    breakindex = 0;
                end

                if F_hist(itr) < F_histbest - 1e-8

                    U_best = U;
                    Y_best = Y;
                    F_histbest = F_hist(itr);
                end

                break;
            end
            zeta = eta*zeta;          nls = nls+1;
        end

        if nupdateU == npU

            tau1 = 0.999;

            if opts.adbeta
                residualnext = max(max(abs(lambda*Y+beta*(Y_now - Y))));
                if (residualnext >= 1e-3 && residualnext >= tau1 * residual)
                    beta0 = 0.9*beta0;
                end

                beta = max(beta0/((itr)^rho),1e-9);
                residual = residualnext;
            else
                beta = beta0/((itr+1)^rho);
            end

            beta_hist(itr) = beta;
            Y_now = Y;

            Y = max(-sparseparameter,min((1/(lambda+beta)*UUt + beta/(lambda+beta)*Y_now),sparseparameter));

        end

        G = L+Y;
        GUUt = (G*U)*U';
        dtU = GUUt + GUUt' - 2*U*(U'*GUUt);     nrmG  = norm(dtU, 'fro');
        S = UUt - UP;
        y = dtU - dtUP;     SY = abs(iprod(S,y));
        if mod(nupdateU,2)==0
            zeta = (norm(S,'fro')^2)/SY;
        else
            zeta  = SY/(norm(y,'fro')^2);
        end

        zeta = max(min(zeta, 1e20), 1e-20);

        Cval = F + 2 * R^2 * beta;
    end

    if nrmG <= epsilon && norm(Y - max(-sparseparameter,min((UUt + Y),sparseparameter)),'fro') <= epsilon
        break;
    end

end

out.nrmG = nrmG;
out.fval = F;
out.itr = itr;
out.obj = F_hist(1:itr);
out.nls_hist = nls_hist(1:itr);
out.beta_hist = beta_hist(1:itr);
out.F_best = F_histbest;
out.Y_best = Y_best;
out.Y_end = Y;
out.U_end = U;

end

function a = iprod(x,y)
%a = real(sum(sum(x.*y)));
a = real(sum(sum(conj(x).*y)));
end