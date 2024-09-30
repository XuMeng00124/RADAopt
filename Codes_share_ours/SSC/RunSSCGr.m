function aver_result = RunSSCGr(d, r, lasso_constant, nrun, b0)

% nrunresults.ARPGDA = [];
nrunresults.RADAPGD = [];
nrunresults.RADARGD = [];
% nrunresults.NADMM = [];
% nrunresults.DSGM = [];

rng(20240000);

for nrunindex = 1:nrun

fprintf('========= r=%d, sparse=%4e, run=%d/%d ==========\n', r, lasso_constant, nrunindex, nrun);

opts.lasso_constant = lasso_constant;
%% Generate L
X = rand(d,d);
W = abs(X'*X);
I = eye(d);
D = diag(sum(W,1));
L = I - sqrt(inv(D))*W*sqrt(inv(D));

%% Generate inital U
d = size(L,1);
[nEigVec,~] = eigs(-L, r, 'largestreal');
U = nEigVec;

opts.epsilon = 1e-4;
opts.R = norm(lasso_constant*ones(d,d),'fro');
opts.n = d*d;
opts.mxitr = 1e4;
opts.nbreakindex = 50;
opts.lambda = opts.epsilon/(2*opts.R);

%% ADMM Result
% opts.lambda0 = 0.01;
% opts.sigma0 = 0.1;
% 
% t_NADMM = 0;
% t = tic;
% [~, ~, out_NADMM] = NADMM(U, L, lasso_constant, opts);
% t_NADMM = t_NADMM + toc(t);
% 
% fprintf(' NADMM & itr=%4d & f=%.8f & time=%fs \n', out_NADMM.itr, out_NADMM.F_best, t_NADMM);
% 
% nrunresults.NADMM=[nrunresults.NADMM; out_NADMM.F_best, t_NADMM, out_NADMM.itr];
%% ARPGDA

% Y0 = ones(d,d);%
% Y0 = Y0 + Y0';
% Y0 = max(-lasso_constant,min(Y0,lasso_constant));
% opts.beta0 = b0;
% opts.rho = 1.5;
% 
% t_ARPGDA = 0;
% t = tic;
% [~, out_ARPGDA]= ARPGDA4SSCGr(U, opts, L, Y0);
% t_ARPGDA =t_ARPGDA +toc(t);
% 
% fprintf('ARPGDA & itr=%4d & f=%.8f & time=%fs \n',out_ARPGDA.itr, out_ARPGDA.F_best, t_ARPGDA);
% nrunresults.ARPGDA = [nrunresults.ARPGDA; out_ARPGDA.F_best, t_ARPGDA, out_ARPGDA.itr];

%% RADA-RGD
Y0 = ones(d,d);%
Y0 = Y0 + Y0';
Y0 = max(-lasso_constant,min(Y0,lasso_constant));

t_RADARGD = 0;
opts.adbeta = 1;
opts.nupdateU = 3;
opts.beta0 = 1;
opts.rho = 1.5;
t = tic;
[~, out_RADARGD]= RADA_RGD4SSCGr(U, opts, L, Y0);
t_RADARGD =t_RADARGD +toc(t);

fprintf('RADA_RGD & itr=%4d & f=%.8f & time=%fs \n',out_RADARGD.itr, out_RADARGD.F_best, t_RADARGD);
nrunresults.RADARGD=[nrunresults.RADARGD; out_RADARGD.F_best, t_RADARGD, out_RADARGD.itr];

%% RADA-PGD

opts.adbeta = 1;
opts.beta0 = 1;
opts.rho = 1.5;
t_RADAPGD = 0;
t = tic;
[~, out_RADAPGD]= RADA_PGD4SSC(U, opts, L, Y0);
t_RADAPGD = t_RADAPGD + toc(t);
fprintf('RADAPGD & itr=%4d & f=%.8f & time=%fs  \n',out_RADAPGD.itr, out_RADAPGD.F_best, t_RADAPGD);
nrunresults.RADAPGD=[nrunresults.RADAPGD; out_RADAPGD.F_best, t_RADAPGD, out_RADAPGD.itr];

%% DSGM 

% t_DSGM = 0;
% opts.lambda0 = 0.01;
% t = tic;
% [~, out_DSGM]= DSGM4SSCGr(U, opts, L, Y0);
% t_DSGM =t_DSGM +toc(t);
% fprintf('DSGM & itr=%4d & f=%.8f & time=%fs \n',out_DSGM.itr, out_DSGM.F_best, t_DSGM);
% 
% nrunresults.DSGM = [nrunresults.DSGM; out_DSGM.F_best, t_DSGM, out_DSGM.itr];

end

% aver_result.NADMM = mean(nrunresults.NADMM, 1); 
% aver_result.DSGM = mean(nrunresults.DSGM, 1); 
% aver_result.ARPGDA = mean(nrunresults.ARPGDA, 1); 
aver_result.RADARGD = mean(nrunresults.RADARGD, 1); 
aver_result.RADAPGD = mean(nrunresults.RADAPGD, 1); 

end