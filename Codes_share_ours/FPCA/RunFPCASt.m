function aver_result = RunFPCASt(d,n,r,nrun,ADMMb0)

% nrunresults.ARPGDA = [];
nrunresults.RADARGD = [];
% nrunresults.RADMM = [];
% nrunresults.DSGM = [];

X = normrnd(0,1,[d,n]);

for nrunindex = 1:nrun

fprintf(1, '\n   --------- #%d/%d, d = %d, N = %d, r = %d ---------\n', nrunindex, nrun, d, n, r);

U0=rand(d,r); %draw U0 randomly
U0=orth(U0);

y0 = rand(n,1);%draw y0 randomly
y0 = simplexproj(y0);

opts=struct();  
opts.mxitr = 20000;
opts.epsilon = 1e-8;
epsilon = opts.epsilon;
R = 1;
opts.lambda = epsilon/(2*R);
opts.n = n;
opts.R = R;
opts.nbreakindex = 1000;

%% ARPGDA

% opts.beta0 = 1e3;opts.rho = 1.5;
% tARPGDA = 0;
% t = tic;
% [~, out_ARPGDA]= ARPGDA4FPCA(U0, opts, X, y0);
% tARPGDA =tARPGDA +toc(t);
% 
% fprintf('ARPGDA: itr=%d, f=%f, time=%fs \n',out_ARPGDA.itr, out_ARPGDA.F_best, tARPGDA);
% nrunresults.ARPGDA=[nrunresults.ARPGDA;out_ARPGDA.F_best,tARPGDA,out_ARPGDA.itr];

%% RADA-RGD

opts.nupdateU = 5;
opts.adbeta = 1;
opts.beta0 = 1e4;opts.rho = 1.5;
tRADARGD = 0;
t = tic;
[~, out_RADARGD]= RADA_RGD4FPCA(U0, opts, X, y0);
tRADARGD =tRADARGD +toc(t);

fprintf('RADA-RGD: itr=%d, f=%f, time=%fs \n',out_RADARGD.itr, out_RADARGD.F_best, tRADARGD);
nrunresults.RADARGD=[nrunresults.RADARGD;out_RADARGD.F_best,tRADARGD,out_RADARGD.itr];
%% RADMM

% t_RADMM = 0;
% opts.lambda0 = 10;
% opts.sigma = ADMMb0;
% t = tic;
% [~, out_RADMM]= RADMM4FPCA(U0, opts, X, y0);
% t_RADMM =t_RADMM +toc(t);
% 
% fprintf('RADMM & itr=%4d & f=%.8f & time=%fs \n',out_RADMM.itr, out_RADMM.F_best, t_RADMM);
% nrunresults.RADMM=[nrunresults.RADMM;out_RADMM.F_best,t_RADMM,out_RADMM.itr];

%% DSGM 

% opts.lambda0 = 10;
% tDSGM = 0;
% t = tic;
% [~, out_DSGM]= DSGM4FPCA(U0, opts, X, y0);
% tDSGM =tDSGM +toc(t);
% fprintf('DSGM: itr=%d, f=%f, time=%fs \n',out_DSGM.itr, out_DSGM.F_best, tDSGM);
% nrunresults.DSGM=[nrunresults.DSGM;out_DSGM.F_best,tDSGM,out_DSGM.itr];

end

% aver_result.RADMM = mean(nrunresults.RADMM, 1); 
% aver_result.DSGM = mean(nrunresults.DSGM, 1); 
% aver_result.ARPGDA = mean(nrunresults.ARPGDA, 1); 
aver_result.RADARGD = mean(nrunresults.RADARGD, 1); 


end