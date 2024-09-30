function [U_best, out]= RADA_RGD4FPCA(U, opts, X, Y)
%% The RADA-RGD algorithm for fair PCA problem
%
% min_{U} max_y -\sum_{i=1}^{m} y_i * < A_i*A_i', U*U' > 
% s.t.  U'*U = I; \sum_{i=1}^{m} y_i = 1; y_i>=0, i=1,2,...,m.
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
%            R: max_y \|y\|, s.t. \sum_{i=1}^{m} y_i = 1; y_i>=0, i=1,2,...,m
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

if ~isfield(opts, 'epsilon');      opts.epsilon = 1e-4; end
if ~isfield(opts, 'zeta');       opts.zeta  = 1e-3; end
if ~isfield(opts, 'rhols');     opts.rhols  = 1e-4; end
if ~isfield(opts, 'eta');       opts.eta  = 0.1; end
if ~isfield(opts, 'mxitr');     opts.mxitr  = 10000; end
if ~isfield(opts, 'R');      opts.R = sqrt(2); end
if ~isfield(opts, 'warmstart');      opts.warmstart = 0; end
if ~isfield(opts, 'nupdateU'); opts.nupdateU = 10; end

%-------------------------------------------------------------------------------
% copy parameters
epsilon = opts.epsilon;
beta0 = opts.beta0;
rhols = opts.rhols;
eta = opts.eta;
npU = opts.nupdateU;
n = opts.n;
R = opts.R;
rho = opts.rho;
lambda = opts.lambda;
Y_now = Y;
beta0 = beta0*n^2*sqrt(r);
beta = beta0;

F_hist = zeros(opts.mxitr,1);
nls_hist = zeros(opts.mxitr,1);
F_histbest = inf;
residual = 1e10;

XtU = X'*U;
f = sum(XtU.*XtU,2);
ytemp = simplexproj(Y_now - 1/(lambda+beta)*f - lambda/(lambda+beta)*Y_now);
F = -f'*ytemp - lambda/2*(ytemp'*ytemp) - beta/2*((ytemp-Y_now)'*(ytemp-Y_now));
yXU = ytemp.*XtU;
G =  X * (-2 *yXU);
GU = G'*U;
dtU = G - U*GU;     nrmG  = norm(dtU, 'fro');

Cval = F;  zeta = opts.zeta;

%% Main Iteration

for itr = 1 : opts.mxitr

    for nupdateU = 1:npU % npU RGD steps to update U

        UP = U;     %FP = F;   GP = G;
        dtUP = dtU;

        nls = 1; deriv = rhols*nrmG^2;

        %% Line Search
        while 1

            [U, ~] = myQR(UP - zeta*dtU, r);

            XtU = X'*U;%FPCA
            f = sum(XtU.*XtU,2);%FPCA
            Y = simplexproj(Y_now - 1/(lambda+beta)*f - lambda/(lambda+beta)*Y_now);%FPCA
            F = -f'*Y - lambda/2*(Y'*Y) - beta/2*((Y-Y_now)'*(Y-Y_now));%FPCA

            if F <= Cval - zeta*deriv || nls >= 5

                F_hist(itr) = -min(f);
                nls_hist(itr) = nls;

                if F_hist(itr) < F_histbest - 1e-8

                    % breakindex = 0;
                    U_best = U;
                    Y_best = Y;
                    F_histbest = F_hist(itr);
                end

                break;
            end
            zeta = eta*zeta;          nls = nls+1;
        end

        % Update Y_k, beta_k
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
            Y = simplexproj(Y_now - 1/(lambda+beta)*f - lambda/(lambda+beta)*Y_now);%FPCA

        end

        yXU = Y.*XtU;G =  X * (-2 *yXU);%FPCA

        GU = G'*U;

        dtU = G - U*GU;     nrmG  = norm(dtU, 'fro');
        S = U - UP;

        y = dtU - dtUP;     SY = abs(iprod(S,y));
        if mod(nupdateU,2)==0
            zeta = (norm(S,'fro')^2)/SY;
        else
            zeta  = SY/(norm(y,'fro')^2);
        end

        zeta = max(min(zeta, 1e20), 1e-20);
        Cval = F + 2 * R^2 * beta;
        % end

    end

    if nrmG <= epsilon && norm(Y - simplexproj(Y - f),'fro') <= epsilon
        break;
    end

end

out.nrmG = nrmG;
out.fval = F;
out.itr = itr;
out.obj = F_hist(1:itr);
out.nls_hist = nls_hist(1:itr);
out.F_best = F_histbest;
% out.lambda = lambda2;
out.Y_best = Y_best;

out.Y_end = Y;

end

function a = iprod(x,y)
%a = real(sum(sum(x.*y)));
a = real(sum(sum(conj(x).*y)));
end

function [Q, RR] = myQR(XX,k)
[Q, RR] = qr(XX, 0);
diagRR = sign(diag(RR)); ndr = diagRR < 0;
if nnz(ndr) > 0
    Q = Q*spdiags(diagRR,0,k,k);
    %Q(:,ndr) = Q(:,ndr)*(-1);
end

end