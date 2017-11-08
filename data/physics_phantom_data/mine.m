% main test codes for sinogram restoration
% i.e., physical phantom data 40mAs

%% clear
% clear all
% close all
clc
%% load data
load physphantom_sinogram.mat 
load physphantom_sigma.mat
sino_raw = exp(-0.00015*max(sino_17mas,0));
sigma = sigma_40mas;
kappa = 2294.5; % scaling factor in MP2012.
GrayWin = [4, 55];
%% imaging geoms
sg = sino_geom('fan','nb', 672, 'na',1160, 'ds', 1.85, ...
    'dsd', 1361.2, 'dod',615.18 ,  'orbit_start',-90,...
    'source_offset',0.0,'channel_offset',-1.25,'orbit',360, 'down', 1);%0.909976
ig = image_geom('nx',512, 'ny', 512,'dx',1.2,'offset_x',0,'down', 1);
G = Gtomo2_dscmex(sg, ig,'mask',ig.mask); % system object
tmp = fbp2(sg, ig);
%% FBP ramp recon
xfbp_ramp = fbp2(sino_raw, tmp,'window','');
Imin = min(xfbp_ramp(:));


imtool(-10000*log(5*(xfbp_ramp-Imin)+0.05))
figure,imshow(xfbp_ramp,GrayWin);title('xfbp ramp');
%% PWLS restoration
sinit=max(sino_raw,0);
niter = 20;
beta_pwls = 30;
sino_pwls= eml_pwls_qm(sinit, sinit, niter,beta_pwls,sigma,kappa);
xfbp_pwls = fbp2(sino_pwls, tmp,'window','');
figure,imshow(xfbp_pwls,GrayWin);title('xfbp pwls');
%% END