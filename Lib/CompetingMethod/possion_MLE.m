function Y = possion_MLE(P,I0,Y0,sigE,beta)
% This script is for the PL method,
% more detail can be found in [1]
%
% [1] Qi Xie, Dong Zeng, Qian Zhao, Deyu Meng*, Zongben Xu, Zhengrong Liang, and Jianhua Ma*
%     Robust Low-dose CT Sinogram Prepocessing via Exploiting Noise-generating Mechanism
%     IEEE Transactions on Medical Imaging, 2017
%
% by Qi Xie, 2017.
sigma = (1./I0);
Y     = (max(log(I0./P),0))*2294.5;
niter = 20;
beta_pwls = 80;
Y = eml_pwls_qm( Y ,  Y , niter,beta_pwls,sigma,2294.5);
Y = gpuArray(Y/2294.5);
P     = gpuArray(P);
I0    = gpuArray(I0);
beta  = gpuArray(beta);

sigE  = gpuArray(sigE);
r = sigE/1;
sizeY = size(Y0);
P  = P+r;
V   = max(abs(I0(:).*P(:).*exp(-Y(:))))*100;
beta = V*beta;
fftnD            = zeros(sizeY,'gpuArray');
ndim = 2;
for i = 1:ndim
    fftnD     = fftnD +permute(( abs(psf2otf([+1; -1], sizeY([i:end,1:i-1])))).^2,[ndim-i+2:ndim,1:ndim-i+1]);
end

for n = 1:100  
    Y = max(real(ifft2(fft2(V.*Y - (I0.*exp(-Y)+r)   +   I0.*exp(Y)./(exp(Y)+r).*P)./(V+beta*fftnD))),0);
end
Y = gather(Y);
end
