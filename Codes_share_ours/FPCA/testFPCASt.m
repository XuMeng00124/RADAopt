%% Test RADA-RGD on Fair Principal Component Analysis
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

d = 1000;

r = 2;

ADMMb0 = 1e-7;
aver11 = RunFPCASt(d,20,r,nrun,ADMMb0); 
aver12 = RunFPCASt(d,25,r,nrun,ADMMb0); 
aver13 = RunFPCASt(d,30,r,nrun,ADMMb0);
aver14 = RunFPCASt(d,35,r,nrun,ADMMb0);
aver15 = RunFPCASt(d,40,r,nrun,ADMMb0);

% save('Result1.mat', 'aver11',"aver12","aver13","aver14","aver15");

n = 20;
r = 4;

aver21 = RunFPCASt(200,n,r,nrun,ADMMb0); 
aver22 = RunFPCASt(400,n,r,nrun,ADMMb0); 
aver23 = RunFPCASt(600,n,r,nrun,ADMMb0);
aver24 = RunFPCASt(800,n,r,nrun,ADMMb0);
aver25 = RunFPCASt(1000,n,r,nrun,ADMMb0);
