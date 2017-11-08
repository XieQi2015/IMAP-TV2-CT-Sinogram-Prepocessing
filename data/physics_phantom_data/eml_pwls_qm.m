function gi= eml_pwls_qm(gi,yi,niter,beta,sigma,kappa)
% gi: sinogram to estimate; 
% yi: measured sinogram
% sigma: 1/Iio;
sigma_elect_noise = 11;
for iter = 1:niter
    [gi_expect,gi_wsum] = pwls_wj(gi);
%     [~,gi_wsum] = pwls_wj(gi);
%     gi_expect = gi;
    sigma_tmp = sigma.*exp(gi_expect/kappa);
    sigma_est = sigma_tmp.*(1+sigma_tmp*(sigma_elect_noise-1.25)); % JH Ma, 2012 Med Phys.
%     sigma_est = sigma.*exp(gi_expect/kappa); % old in J Wang
    gi=(yi+beta*sigma_est.*gi_wsum)./(1+beta*sigma_est); % J Wang,2006 TMI.
    gi=max(gi,0);
end

function [ye,yw] = pwls_wj(yi)
[nx,ny] = size(yi);
w = 1; % 3x3 window
weight = [0 1 0;0.25 0 0.25;0 1 0]/2.5; %  J Wang, PWLS: 2006 TMI.
yi_ext = symmetric_extension(yi,w);
p = nx * ny;
[Y,X] = meshgrid(w+1:w+ny, w+1:w+nx);
X = X(:);
Y = Y(:);
[dY,dX] = meshgrid(-w:w,-w:w);
dim = 2*w+1;
Xp = repmat( reshape(X,[1,1,p]) ,[dim dim 1]) + repmat(dX,[1 1 p]);
Yp = repmat( reshape(Y,[1,1,p]) ,[dim dim 1]) + repmat(dY,[1 1 p]);
I = sub2ind(size(yi_ext), Xp,Yp);
Yi = yi_ext(I);
Yi = reshape(Yi,[dim^2,p]);
ye = mean(Yi,1);
ye = reshape(ye,[nx,ny]);
yw = Yi' * weight(:);
yw = reshape(yw,[nx,ny]);


function xi_padded = symmetric_extension(xi,k)
% symmetric_extension - perform a symmetric extension of the signal.
n1 = size(xi,1);
n2 = size(xi,2);
xi_padded = zeros(n1+2*k,n2+2*k);
xi_padded(k+1:end-k,k+1:end-k) = xi;
% extension
xi_padded(1:k,:) = xi_padded(2*k:-1:k+1,:);
xi_padded(end-k+1:end,:) = xi_padded(end-k:-1:end-2*k+1,:);
xi_padded(:,1:k) = xi_padded(:,2*k:-1:k+1);
xi_padded(:,end-k+1:end) = xi_padded(:,end-k:-1:end-2*k+1);