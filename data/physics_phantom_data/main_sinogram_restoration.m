% main test codes for sinogram restoration
% i.e., physical phantom data 40mAs

%% clear
% clear all
% close all
clc
%% load data
load physphantom_sinogram.mat 
load physphantom_sigma.mat
sino_raw =sino_40mas;
sigma = sigma_40mas;
kappa = 2294.5; % scaling factor in MP2012.
GrayWin = [4, 55];
%% imaging geoms
sg = sino_geom('fan','nb', 672, 'na',1160, 'ds', 1.85, ...
    'dsd', 1361.2, 'dod',615.18 ,  'orbit_start',-90,...
    'source_offset',0.0,'channel_offset',-1.25,'orbit',360, 'down', 1);%0.909976
ig = image_geom('nx',512, 'ny', 512,'dx',1.2,'offset_x',0,'down', 1);
% G = Gtomo2_dscmex(sg, ig,'mask',ig.mask); % system object
tmp = fbp2(sg, ig);
%% FBP ramp recon
% xfbp_ramp = fbp2(sino_raw, tmp,'window','');
% figure,imshow(xfbp_ramp,GrayWin);
%% PWLS restoration
sinit=max(sino_raw,0);
niter = 20;
beta_pwls = 30;
tic
sino_pwls= eml_pwls_qm(sinit, sinit, niter,beta_pwls,sigma,kappa);
toc
xfbp_pwls = fbp2(sino_pwls, tmp,'window','');
figure,imshow(xfbp_pwls,GrayWin);
minI = min(xfbp_true(:));
maxI = max(xfbp_true(:));
[psnrPwls, ssimPwls, fsimPwls, ergasPwls, msamPwls] = MSIQA(normalized(xfbp_pwls)*255, normalized(xfbp_true)*255);
disp(['Pwls·¨µÄ PSNR Îª£º  ' num2str(psnrPwls)])

%% END