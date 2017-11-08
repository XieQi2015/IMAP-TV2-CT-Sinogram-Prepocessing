% physical phantom FBP recons 17 ~ 100 mAs
%%
close all
clear all
clc
%%
%% reading data
load physphantom_sinogram
%% 17 mAs
%% imaging geoms
sg = sino_geom('fan','nb', 672, 'na',1160, 'ds', 1.85, ...
    'dsd', 1361.2, 'dod',615.18 ,  'orbit_start',-90,...
    'source_offset',0.0,'channel_offset',-1.25,'orbit',360, 'down', 1);
ig = image_geom('nx',512, 'ny', 512,'dx',1.2,'offset_x',0,'down', 1);
% % % % system matrix
G = Gtomo2_dscmex(sg, ig); % linux system
% % % % display window
GrayWin = [4, 55];
%% FBP recons with ramp window
xfbp_17mas = fbp_fan_arc(sino_17mas, G,'');
%% figures show
figure,imshow(xfbp_17mas,GrayWin);
%%
%% 40 mAs
%% imaging geoms
sg = sino_geom('fan','nb', 672, 'na',1160, 'ds', 1.85, ...
    'dsd', 1361.2, 'dod',615.18 ,  'orbit_start',-90,...
    'source_offset',0.0,'channel_offset',-1.25,'orbit',360, 'down', 1);
ig = image_geom('nx',512, 'ny', 512,'dx',1.2,'offset_x',0,'down', 1);
% % % % system matrix
G = Gtomo2_dscmex(sg, ig); % linux system
% % % % display window
GrayWin = [4, 55];
%% FBP recons with ramp window
xfbp_40mas = fbp_fan_arc(sino_40mas, G,'');
%% figures show
figure,imshow(xfbp_40mas,GrayWin);
%%
%% 60 mAs
%% imaging geoms
sg = sino_geom('fan','nb', 672, 'na',1160, 'ds', 1.85, ...
    'dsd', 1361.2, 'dod',615.18 ,  'orbit_start',90,...
    'source_offset',0.0,'channel_offset',-1.25,'orbit',360, 'down', 1);
ig = image_geom('nx',512, 'ny', 512,'dx',1.2,'offset_x',0,'down', 1);
% % % % system matrix
G = Gtomo2_dscmex(sg, ig); % linux system
% % % % display window
GrayWin = [4, 55];
%% FBP recons with ramp window
xfbp_60mas = fbp_fan_arc(sino_60mas, G,'');
%% figures show
figure,imshow(xfbp_60mas,GrayWin);
%%
%% 100 mAs
%% imaging geoms
sg = sino_geom('fan','nb', 672, 'na',1160, 'ds', 1.85, ...
    'dsd', 1361.2, 'dod',615.18 ,  'orbit_start',-90+360/1160*2,...
    'source_offset',0.0,'channel_offset',-1.25,'orbit',360, 'down', 1);
ig = image_geom('nx',512, 'ny', 512,'dx',1.2,'offset_x',0,'down', 1);
% % % % system matrix
G = Gtomo2_dscmex(sg, ig); % linux system
% % % % display window
GrayWin = [4, 55];
%% FBP recons with ramp window
xfbp_100mas = fbp_fan_arc(sino_100mas, G,'');
%% figures show
figure,imshow(xfbp_100mas,GrayWin);
%%
%% gold standard
%% imaging geoms
sg = sino_geom('fan','nb', 672, 'na',1160, 'ds', 1.85, ...
    'dsd', 1361.2, 'dod',615.18 ,  'orbit_start',90,...
    'source_offset',0.0,'channel_offset',-1.25,'orbit',360, 'down', 1);
ig = image_geom('nx',512, 'ny', 512,'dx',1.2,'offset_x',0,'down', 1);
% % % % system matrix
G = Gtomo2_dscmex(sg, ig); % linux system
% % % % display window
GrayWin = [4, 55];
%% FBP recons with ramp window
xfbp_true = fbp_fan_arc(sino_100mas_mean, G,'');
%% figures show
figure,imshow(xfbp_true,GrayWin);
%% figures show
close all
figure,imshow(xfbp_true,GrayWin);
figure,imshow(xfbp_17mas,GrayWin);
figure,imshow(xfbp_40mas,GrayWin);
figure,imshow(xfbp_60mas,GrayWin);
figure,imshow(xfbp_100mas,GrayWin);
%% END
save physphantom_fbp_recon xfbp_true xfbp_17mas xfbp_40mas xfbp_60mas xfbp_100mas GrayWin