# IMAP-TV2-CT-Sinogram-Prepocessing
Outline:

    Folder structure
    Usage
    Citation
    URL of third-party toolboxes and functions


Folder structure:

    data\           : example data
        |-- physics_phantom_data\      : physics_phantom_data
    lib\            : library of code
        |-- KBRreg\            : helper functions of the ITSReg method
        |-- MyTool\            : a toolbox for the perposed method 
        |-- quality_assess\    : functions of quality assessment indices
        |       |-- ssim_index.m            : SSIM [1]
        |       |-- FeatureSIM.m            : FSIM [2]
        |       |-- ErrRelGlobAdimSyn.m     : ERGAS
        |       |-- SpectAngMapper.m        : SAM
        |       |-- MSIQA.m                 : interface for calculating PSNR and the four indices above
        |-- CompetingMethod : competing methods (down loaded or implemented based on reference papers)
        |       |-- eml_pwls_qm.m           : script for PWLS method [3]
        |       |-- possion_MLE.m           : script for PL method, GPU version [4]
        |       |-- possion_MLE_CPU.m       : script for PL method, CPU version [4]
        |       |-- pwls_tv.m               : script for POCS-TV method [5]
        |       |-- PWLS_TV2.m              : script for PWLS-TV2 method, GPU version 
        |       |-- PWLS_TV2_CPU.m          : script for PWLS-TV2 method, CPU version
        |       |-- IMAP_TV.m               : script for IMAP-TV method, GPU version
        |       |-- IMAP_TV_CPU.m           : script for IMAP-TV method, CPU version  
        |       |-- ParSet_compare.m        : script for setting parameters in PWLS-TV2 method
        |       |-- ParSetIMAP_TV.m         : script for setting parameters in IMAP-TV method
        |-- GetPic_Anthropomorphic          : function for showing result of the example
    Demo.m          : scripts that applies the methods and calculates the QA indices
    IMAP_TV2_GPU.m  : core function of the proposed IMAP-TV2 sinogram preprocessing, GPU version
    IMAP_TV2_CPU.m  : core function of the proposed IMAP-TV2 sinogram preprocessing, CPU version

Usage:

    For CT sinogram, you can simply follow these steps:
        1.If you have install the Michigan Image Reconstruction Toolbox (MIRT), please skip the first step.
          Else, download the MIRT form the following link:
          'http://www.eecs.umich.edu/~fessler/',
          then install the MIRT.
        2.Prepare the the corrupted projection data,the incident X-ray intensity;
        3.Add 'Lib' into path;
        4.Use the function IMAP_TV2_GPU (or IMAP_TV2_CPU) as follows:
            [ clean_img ] = IMAP_TV2_GPU(P,I0)
        5.If you want to change the parameters in the method, use
            [ clean_img ] = IMAP_TV2_GPU(P,I0,b,Par)
          to set a initial b, and other parameters.
    Please type 'help IMAP_TV2_GPU' to get more information.

    You may find example codes in file Demo.m

    Also, you can use the demo to see some comparison. You can:
      1. Type 'Demo' to to run various methods and see the pre-computed results.
      2. Use 'help Demo' for more information.
      3. Select competing methods by turn on/off the enable-bits in Demo.m

Citation:

     Qi Xie, Dong Zeng, Qian Zhao, Deyu Meng*, Zongben Xu, Zhengrong Liang, and Jianhua Ma*
     Robust Low-dose CT Sinogram Preprocessing  via Exploiting Noise-generating Mechanism
     IEEE Transactions on Medical Imaging, 2017

URL of the toolboxes and functions and citation of the competing method:

    [1]  ssim_index.m       https://ece.uwaterloo.ca/~z70wang/research/ssim/
    [2]  FeatureSIM.m       http://www4.comp.polyu.edu.hk/~cslzhang/IQA/FSIM/FSIM.htm
    [3]  P. J. La Rivi`ere, J. Bian, and P. A. Vargas. Penalized-likelihood sinogram restoration for computed tomography. 
         IEEE Transactions on Medical Imaging, 25(8):1022–1036, 2006.
    [4]  J. Wang, T. Li, H. Lu, and Z. Liang. Penalized weighted least-squares approach to sinogram noise reduction and 
         image reconstruction for low-dose x-ray computed tomography. IEEE Transactions on Medical Imaging, 25(10):1272–1283, 2006.
    [5]  E. Y. Sidky and X. Pan. Image reconstruction in circular conebeam computed tomography by constrained, total-variation 
         minimization. Physics in Medicine & Biology, 53(17):4777–4807, 2008.
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Resources for research purpose only, shall not be used for commercial purposes! All copyrights belong to the original anthors. 
The technology has applied for patents. If you want to purchase the patents for commercial purposes, please contact the 
corresponding author: Deyu Meng, dymeng@mail.xjtu.edu.cn. Thank you!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

