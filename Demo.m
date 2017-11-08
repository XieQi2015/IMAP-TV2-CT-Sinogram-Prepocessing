%==========================================================================
% This script compares CT Sinogram Preprocessing methods
% listed as follows:
%   1. PL
%   2. PWLS
%   3. POCS-TV
%   4. PWLS-TV2
%   5. IMAP-TV
%   6. IMAP-TV2 (the proposed method)
%
% Three quality assessment (QA) indices -- PSNR, FSIM, NMSE
% -- are calculated for each methods after Preprocessing.
%
% You can:
%       1. Type 'Demo' to to run various methods and see the pre-computed results.
%       2. Select competing methods by turn on/off the "enList" in Demo.m
%       3. Select competing case by simple change "NmasList" in Demo.m
% See also 
% IMAP_TV2_CPU, IMAP_TV2_GPU, IMAP_TV, PWLS_TV2, pwls_tv, eml_pwls_qm, possion_MLE.
%
% more detail can be found in [1]
%
% [1] Qi Xie, Dong Zeng, Qian Zhao, Deyu Meng*, Zongben Xu, Zhengrong Liang, and Jianhua Ma*
%     Robust Low-dose CT Sinogram Prepocessing via Exploiting Noise-generating Mechanism
%     IEEE Transactions on Medical Imaging, 2017
%
% by Qi Xie, 2017.
%==========================================================================

clc;
clear;close all;
%% data loading
addpath(genpath('Lib'));
load('data\physics_phantom_data\physphantom_sinogram.mat');
load('data\physics_phantom_data\physphantom_sigma.mat');
load('data\physics_phantom_data\physphantom_fbp_recon.mat');
load('data\physics_phantom_data\PWLS_rec.mat');
SaveName = 'Result_1';
%% input

NameMethod  = {   'PL',      'PWLS',    'POCS_TV',      'PWLS-TV2',     'IMAP-TV',      'IMAP-TV2'    }   ;
enList      = [    1,           1,            0,               1,             1,               1]  >.5    ; % set to 0 for turning off
NmasList    = [40];   %case to use, can use a subset of [17,40,60]
UseGPU      = 1;

%% initialization
X_true      = xfbp_true;
PSNR        = zeros(length(NmasList),sum(enList));
SSIM        = zeros(length(NmasList),sum(enList));
FSIM        = zeros(length(NmasList),sum(enList));
NMSE        = zeros(length(NmasList),sum(enList));
IndList     = 1:length(NameMethod);
IndList     = IndList(enList);
 
% imaging geoms
GrayWin = [-12, 90];
sg = sino_geom('fan','nb', 672, 'na',1160, 'ds', 1.85, ...
    'dsd', 1361.2, 'dod',615.18 ,  'orbit_start',-90,...
    'source_offset',0.0,'channel_offset',-1.25,'orbit',360, 'down', 1);%0.909976
ig = image_geom('nx',512, 'ny', 512,'dx',1.2,'offset_x',0,'down', 1);
G = Gtomo2_dscmex(sg, ig,'mask',ig.mask); % system object
tmp = fbp2(sg, ig);

%% main
if UseGPU
    getGPU = 1;
    disp('Setting up GPU...')
    tic
    gpuArray(getGPU);
    GPUtime = toc;
    disp(['Take ' num2str(toc) 's to set up GPU']);
end

for i  = 1:length(NmasList)
    Nmas   = NmasList(i);
    disp('==============================================='  )
    disp(['===========  Case of ' num2str(Nmas),  ' mAs  ============' ] )
    
    signY    = eval(['sino_' num2str(Nmas) 'mas']);
    sigma    = eval(['sigma_' num2str(Nmas) 'mas']);
    Xnois{i} = eval(['xfbp_' num2str(Nmas) 'mas']);
    sizeY    = size(signY);
    sigma    = mean(sigma,2);
    temp     = polyfit((1:sizeY(1))',(sigma+sigma(end:-1:1))/2,20);
    tempSi   = polyval(temp,(1:sizeY(1))');
    sigma    = repmat(tempSi,1,sizeY(2));
    I0       = 1./sigma;
    signY    = min(max(signY/ 2294.5, 0),log(I0));
    signP    = exp(-signY)./sigma;
     
    nM = 0;
    
    %% Using PL method
    nM = nM+1;
    if enList(nM)
        disp('...')
        disp(['performing ',NameMethod{nM}, ' ... ']);
        if UseGPU
            tic
            %         Y_MLE = possion_MLE(signP,I0,signY,11,0.0005);% 17 mAs
            Y_MLE = possion_MLE(signP,I0,signY,11,0.0008);% 40 mAs
            %             Y_MLE = possion_MLE(signP,I0,signY,11,0.0012); % 60 mAs
            Time = toc;
        else
            tic
            %         Y_MLE = possion_MLE_CPU(signP,I0,signY,11,0.0005);% 17 mAs
            Y_MLE = possion_MLE_CPU(signP,I0,signY,11,0.0008);% 40 mAs
            %             Y_MLE = possion_MLE_CPU(signP,I0,signY,11,0.0012); % 60 mAs
            Time = toc;
        end

        disp([NameMethod{nM}, ' done in ' num2str(Time), ' s.'])
        X{i,nM} = fbp2(Y_MLE*2294.5, tmp,'window','');
    end
    
    %% Using PWLS method
    nM = nM+1;
    if enList(nM)
        disp('...')
        disp(['performing ',NameMethod{nM}, ' ... ']);
        niter = 20;
        beta_pwls = 250;%
        tic
        Y_PWLS = eml_pwls_qm(signY*2294.5, signY*2294.5, niter,beta_pwls,sigma,2294);
        Time = toc;
        disp([NameMethod{nM}, ' done in ' num2str(Time), ' s.'])
        X{i,nM}   = fbp2(Y_PWLS, tmp,'window','');
    end
    
    %% Using POCS-TV method
    nM = nM+1;
    if enList(nM)
        disp('...')
        disp(['performing ',NameMethod{nM}, ' ... ']);
        X_its = pwls_tv(X_noise, signY*kappa, G, 1e-2, 1e4, a);
        X{i,nM} = X_its;
    end
    
    %% Using PWLS-TV^2 method
    nM = nM+1;
    if enList(nM)
        disp('...')
        disp(['performing ',NameMethod{nM}, ' ... ']);
        if UseGPU
            tic
            [Par,beta] = ParSet_compare(Nmas);
            [Y_PWLS_FC, Var,Z]= PWLS_TV2(signP,I0,beta,Par);
            Time = toc;
        else
            tic
            [Par,beta] = ParSet_compare(Nmas);
            [Y_PWLS_FC, Var,Z]= PWLS_TV2_CPU(signP,I0,beta,Par);
            Time = toc;
        end
        disp([NameMethod{nM}, ' done in ' num2str(Time), ' s.'])
        X{i,nM} = fbp2(Y_PWLS_FC*2294.5, tmp,'window','');
    end
    
    %% Using IMAP-TV method
    nM = nM+1;
    if enList(nM)
        disp('...')
        disp(['performing ',NameMethod{nM}, ' ... ']);
        if UseGPU
            tic
            Par = ParSetIMAP_TV(Nmas,1);
            [Y_diff1, I_diff1]= IMAP_TV(signP,I0,Par);
            Time = toc;
        else
            tic
            Par = ParSetIMAP_TV(Nmas,1);
            [Y_diff1, I_diff1]= IMAP_TV_CPU(signP,I0,Par);
            Time = toc;
        end
        disp([NameMethod{nM}, ' done in ' num2str(Time), ' s.'])
        X{i,nM} = fbp2(Y_diff1*2294.5, tmp,'window','');
    end
    
    %% Using IMAP-TV^2 method
    nM = nM+1;
    if enList(nM)
        disp('...')
        disp(['performing ',NameMethod{nM}, ' ... ']);
        if UseGPU
            tic
            [Par,b] = ParSet(Nmas);
            Y = IMAP_TV2_GPU(signP,I0,b,Par);
            Time = toc;
        else
            tic
            [Par,b] = ParSet(Nmas);
            Y = IMAP_TV2_CPU(signP,I0,b,Par);
            Time = toc;
        end
        disp([NameMethod{nM}, ' done in ' num2str(Time), ' s.'])
        X{i,nM} = fbp2(Y*2294.5, tmp,'window','');
    end
    
    %% caculating performance measurements
    minI = min(X_true(:));
    maxI = max(X_true(:));
    for j = IndList(end:-1:1)
        if Nmas == 60
            X{i,j} =  X{i,j}(end:-1:1,end:-1:1);
        end
        [PSNR(i,j),SSIM(i,j),FSIM(i,j)] = MSIQA((X{i,j}-minI)/(maxI-minI)*255, normalized(X_true)*255);
        NMSE(i,j)   = norm(X{i,j}(:)-X_true(:))/norm(X_true(:));
    end
    
    
end
saveroad0  = ['result/', SaveName];
mkdir(saveroad0);
save([saveroad0, '/Result'],'X','PSNR','SSIM','FSIM');
for i = 1: length(NmasList)
fprintf('\n');
disp(['=============== Result of', num2str(Nmas),  ' mAs ===============']);
fprintf(' %9.9s    %5.4s       %5.4s    %5.5s    \n','method','PSNR', 'FSIM', ' NMSE');
    for j = 1:length(NameMethod)
        if enList(j)
            fprintf(' %9.9s    %5.3f       %5.3f    %5.3f    \n',...
                NameMethod{j}, PSNR(i,j), FSIM(i,j), NMSE(i,j));
        end
    end
end
fprintf('================== Result =====================\n');
close all
GetPic_Anthropomorphic (X, X_true, Xnois, NameMethod, NmasList, enList)


