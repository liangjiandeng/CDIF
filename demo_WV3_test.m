%% ================================================================
% This is the demo code for 
% "A New Context-Aware Details Injection Fidelity with 
%  Adaptive Coefficients Estimation for Variational Pansharpening"
%  by J.-L. Xiao, T.-Z. Huang, L.-J. Deng, Z.-C. Wu and G. Vivone.

% If you use this code, please cite the following paper:

% J.-L. Xiao, T.-Z. Huang, L.-J. Deng, Z.-C. Wu and G. Vivone, 
% A New Context-Aware Details Injection Fidelity with Adaptive Coefficients
% Estimation for Variational Pansharpening,
% IEEE Trans. Geosci. Remote Sens., doi:10.1109/TGRS.2022.3154480.

%==========================================================================

clear;
clc;
close all;
addpath(genpath(pwd));
%% Load data which type is double precision and the range is [0 1]
load 'WV3_test.mat';

%% Initialization
maxit = 100;             
lambda = 0.00005;
lambda2 = 0.0000001;
lambda3 = 0.1;
eta_1 = 0.00007;
eta_2 = 0.0000000001;
eta_3 = 0.0000000001;
eta_4 = 0.000000009;
ratio = 4;
nclusters =5;%The number of clusters
sf      = 4;
sensor = 'WV3';
[~,~,L]  = size(lrms);
sz       = size(pan);
Nways    = [sz, L];
opts.sf     = sf;
opts.Nways  = Nways;
opts.sensor = sensor;
opts.nclusters = nclusters;
opts.tol    = 2*1e-5;
opts.lambda = lambda;
opts.lambda2 = lambda2;
opts.lambda3 = lambda3;
opts.eta_1 = eta_1;
opts.eta_2 = eta_2;
opts.eta_3 = eta_3;
opts.eta_4 = eta_4;
opts.maxit = maxit;
opts.sz = sz;

%% Adaptive Coefficients Estimation
[G1,G2] = G_Estimate(lrms,pan ,opts);

%% The Fusion Algorithm
I_guide = CDIF_fusion(lrms, pan,G1,G2, opts);

%% Plotting
close all
location = [140 160 100 120];
showRGB8(gt, gt, location);title('Orginal');
showRGB8(gt, I_guide, location);title('Fusion by CDIF');
%% Show Metrics
Eva_Xfin   = Quality_assess(gt, I_guide, sf);
fprintf('CDIF_test      PSNR: %.4f   SSIM: %.4f   SAM: %.4f   SCC: %.4f   ERGAS: %.4f   Q8: %.4f\n',...
     Eva_Xfin.PSNR,Eva_Xfin.SSIM,Eva_Xfin.SAM,Eva_Xfin.SCC,Eva_Xfin.ERGAS,Eva_Xfin.Q8)