%% Test RADA on Sparse Spectral Clustering
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

close all;

nrun = 20;

d = 200; sparse = 5e-3;
aver11 = RunSSCGr(d, 2, sparse, nrun, 1e-7);
aver12 = RunSSCGr(d, 3, sparse, nrun, 1e-8);
aver13 = RunSSCGr(d, 4, sparse, nrun, 1e-8);
aver14 = RunSSCGr(d, 5, sparse, nrun, 1e-8);
aver15 = RunSSCGr(d, 6, sparse, nrun, 1e-8);
% save('Result1tk3.mat', 'aver11',"aver12","aver13","aver14","aver15");

d = 500; r = 3;
aver21 = RunSSCGr(d, r, 1e-3, nrun, 1e-9);%
aver22 = RunSSCGr(d, r, 2e-3, nrun, 1e-9);%
aver23 = RunSSCGr(d, r, 5e-3, nrun, 1e-9);%
aver24 = RunSSCGr(d, r, 1e-2, nrun, 1e-10); %
aver25 = RunSSCGr(d, r, 2e-2, nrun, 1e-10);%
%save('Result2tk3.mat', 'aver21',"aver22","aver23","aver24","aver25");
