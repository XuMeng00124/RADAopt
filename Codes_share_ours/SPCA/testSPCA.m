%% Test the RADA-RGD Algorithm on Sparse Principal Component Analysis
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
clear;

nrun = 20;

d = 1000;
Result11 = RunRep(d,  50, 10, 0.5, 1, nrun);
Result12 = RunRep(d,  50, 10, 0.75, 1, nrun);
Result13 = RunRep(d,  50, 10, 1.00, 1, nrun);
Result14 = RunRep(d,  50, 10, 1.25, 1, nrun);
Result15 = RunRep(d,  50, 10, 1.5, 1, nrun);

%save('Result1.mat', 'Result11',"Result12","Result13","Result14","Result15");

d = 2000; mu = 3;
Result21 = RunRep(d, 50, 4, mu, 1, nrun);
Result22 = RunRep(d, 50, 6, mu, 1, nrun);
Result23 = RunRep(d, 50, 8, mu, 1, nrun);
Result24 = RunRep(d, 50, 10, mu, 1, nrun);
Result25 = RunRep(d, 50, 12, mu, 1, nrun);

%save('Result2.mat', 'Result21',"Result22","Result23","Result24","Result25");

function Result = RunRep(n, m, r, mu, start, reps)

T1 = [];
S1 = [];
F1 = [];
I1 = [];
V1 = [];

for num = start:reps

    fprintf(1, '\n   --------- #%d/%d, n = %d, r = %d, mu = %.2f ---------\n', num, reps, n, r, mu);

    A = randn(m, n);
    A = A - repmat(mean(A, 1), m, 1);
    A = A ./ repmat(sqrt(sum(A .* A)), m, 1);

    [U, ~, V] = svd(A, 'econ');
    % D = diag(S(1:r, 1:r));
    D = sort(abs(randn([m 1])) .^ 4) + 1.0e-5;

    A = U * diag(D) * V';
    A = A - repmat(mean(A, 1), m, 1);
    A = A ./ repmat(sqrt(sum(A .* A)), m, 1);

    [V, D] = eig(A'*A);

    eigVals = diag(D);
    [~, idx] = sort(eigVals, 'descend');
    Xpca = V(:, idx(1:r));
    Var0 = sum(Xpca.*(A'*(A*Xpca)),"all");

    [X_init,~] = svd(randn(n,r),0);  % random intialization

    % ================== RADA-RGD ====================

    Y0 = ones(size(X_init));
    opts.n = n*r;
    opts.R = norm(mu*ones(n,r),'fro');
    opts.epsilon = 1e-8;
    lambda_RADA = opts.epsilon/(2*opts.R);

    opts.mxitr = 50000;
    opts.lasso_constant = mu;

    opts.adbeta = 1;
    opts.nupdateU = 10;
    opts.beta0 = 1;
    opts.rho = 1.5;
    t_RADA = 0;
    t = tic;

    [U_RADA, out_RADA]= RADA_RGD4SPCA(X_init, opts, A, lambda_RADA, Y0);
    t_RADA =t_RADA +toc(t);

    fprintf('RADA-RGD & itr=%4d & f=%.8f & time=%fs & sparse=%.4f \\\\ \n',out_RADA.itr, out_RADA.F_best,...
        t_RADA, out_RADA.sparse);
    T1 = [T1; t_RADA];
    F1 = [F1; out_RADA.F_best];
    I1 = [I1; out_RADA.itr];
    S1 = [S1; out_RADA.sparse * 100];
    V1 = [V1; sum(U_RADA.*(A'*(A*U_RADA)),"all")/Var0];

end

fprintf(1, '=========== Summary: n = %d, r = %d, mu = %.3f ==========\n', n, r, mu);
fprintf(1, 'RADA-RGD:   time = %.3fs, sparsity = %.2f, loss = %.4f, iter = %.4f\n', mean(T1), mean(S1), mean(F1),mean(I1));

Result = [mean(F1);
    mean(V1);
    mean(I1);
    mean(T1);
    mean(S1)];

end