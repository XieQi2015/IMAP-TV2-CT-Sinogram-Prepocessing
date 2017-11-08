function  [Y,var,Z] = PWLS_TV2_CPU(SignoP,I0,beta,Par)
% This script is for the PWLS-TV2 method,
% more detail can be found in [1]
%
% [1] Qi Xie, Dong Zeng, Qian Zhao, Deyu Meng*, Zongben Xu, Zhengrong Liang, and Jianhua Ma*
%     Robust Low-dose CT Sinogram Prepocessing via Exploiting Noise-generating Mechanism
%     IEEE Transactions on Medical Imaging, 2017
%
% by Qi Xie, 2017.
sizeS  = size(SignoP);
N      = prod(sizeS);
ndim   = length(sizeS);
if nargin <3
    beta  = 1;
end
if nargin<4
    MaxIter        = 1000;
    tol            = 1e-7;
    mu             = 1e-5;
    rho            = 1.05; 
    initIterNum    = 40;
    var            = 11;
    wTV            = [1,1,1,1];
    
else
    if isfield(Par,'maxIter')         ; MaxIter = Par.maxIter;          else  MaxIter = 1000;          end;
    if isfield(Par,'tol')             ; tol = Par.tol;                  else  tol = 1e-7;              end;
    if isfield(Par,'mu')              ; mu = Par.mu;                    else  mu = 1;                  end;
    if isfield(Par,'rho')             ; rho = Par.rho;                  else  rho = 1.05;              end;
    if isfield(Par,'initIterNum')     ; initIterNum = Par.initIterNum;  else  initIterNum = 40;        end;
    if isfield(Par,'wTV')             ; wTV = Par.wTV;                  else  wTV = [1,1,1,1];         end;
    if isfield(Par,'var')             ; var = Par.var;                  else  var = 11;                end;
end

sigma = gather(1./I0);
Y0     = gather(max(log(I0./SignoP),0))*Par.kappa;
niter = 20;
beta_pwls = Par.pwlsBeta;
[Y, sigma_est]  = eml_pwls_qm( Y0 ,  Y0 , niter,beta_pwls,sigma,Par.kappa);
Y = (Y/Par.kappa);
L        = 1*max(max(sigma_est(:), eps).^(-2))/2;
Dy               = diff2_ND(Y,sizeS);

Q                = zeros(size(Dy));
Z                = Q;
wbeta            = Z;
for k = 1:4
    wbeta(N*(k-1)+1:N*k) = beta*wTV(k);
end
normD            = sum(Dy.^2);
normD            = sqrt(normD);
fftnD            = zeros(sizeS);
for i = 1:ndim
    fftnD     = fftnD + permute((abs(psf2otf([+1, -1 ; -1, +1], sizeS([i:end,1:i-1])))).^2,[ndim-i+2:ndim,1:ndim-i+1]);
    fftnD     = fftnD + permute((abs(psf2otf([+1; -2 ; +1], sizeS([i:end,1:i-1])))).^2,[ndim-i+2:ndim,1:ndim-i+1]);
end
subtol = 0.3e-2;
theEps = 1e-6;
%% Main loop
for i = 1:MaxIter
    
    %% UPdate Z
    if i<=initIterNum
        Z  = (mu*Dy+Q)./(wbeta+mu);                        
    else
        Z  = ClosedWL1(Dy+Q/mu,wbeta/mu,theEps);
    end
    %% Update Y
    Y       = UpdateY(Y,Y0,sigma_est,L,mu,Z,Q,fftnD,sizeS);
    Dy      = diff2_ND(Y,sizeS);
    %% Updtae beta
    if  i>40 
        for k = 1:4
            wbeta(N*(k-1)+1:N*k) = N/sum(log(theEps+abs(Dy(N*(k-1)+1:N*k))) -log(theEps));
        end
    end
    
   
    temp    = Dy-Z;
    Q    = Q+mu*(temp);
    reErr   =  sum(temp(:).^2)/normD;
    mu  = mu*rho;
  
    if reErr<tol
        disp('convergence');
        break
    end
end
Y   = gather(Y);
var = gather(var);
Z   = gather(Z);
if i==MaxIter
    disp(['not convergence，relative error is ', num2str(reErr)])
end
end

function Y   =UpdateY(Y_ini,Y0,sigma_est,L,mu,Z,Q,fftnD,sizeS)
%更新y，用FSTA法

t        = 1;
stoption = 0;
difZ     = diff2T_ND(Z-Q/mu,sizeS);
Y        = Y_ini;
X        = Y;
% normY    = norm(Y(:));
i        = 1;
while ~stoption
    perX     = X;
    perT     = t;
    tempX    = 1/mu*(L*Y-sigma_est.*(Y-Y0))+difZ;
    X        = ifft2(fft2(tempX) ./(fftnD + L/mu) );
    X        = max(0,real(X));
    t        = ( 1+sqrt(1+ 4 * perT^2 ) )/2;
    Y        = X + (perT/t).*(X-perX);
    %      reErr    = sqrt(sum((X(:)-perX(:)).^2))/normY;
    i        = i+1;
    stoption = i>30;%&&(reErr<tol);
end
end

