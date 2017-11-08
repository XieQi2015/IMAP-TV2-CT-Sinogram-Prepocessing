function  [Y,I,var,Z] = IMAP_TV_CPU(SignoP,I0,Par)
% This script is for the IMAP-TV method,
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
if nargin<3
    beta           = 1;
    MaxIter        = 1000;
    tol            = 1e-7;
    mu             = 1e-5;
    rho            = 1.05; 
    Y                = max(log(I0./SignoP),0);
    Dy               = diff_ND(Y,sizeS);
else
    if isfield(Par,'beta')            ; beta = Par.beta;                else  beta = 1;                end;
    if isfield(Par,'maxIter')         ; MaxIter = Par.maxIter;          else  MaxIter = 1000;          end;
    if isfield(Par,'tol')             ; tol = Par.tol;                  else  tol = 1e-7;              end;
    if isfield(Par,'mu')              ; mu = Par.mu;                    else  mu = 1;                  end;
    if isfield(Par,'rho')             ; rho = Par.rho;                  else  rho = 1.05;              end;
    if isfield(Par,'var')             ; var = Par.var;                  else  var = 11;                end;
    if isfield(Par,'PwlsInit')&&Par.PwlsInit==1;
        sigma = gather(1./I0);
        Y     = gather(max(log(I0./SignoP),0))*Par.kappa;
        niter = 20;
        beta_pwls = Par.pwlsBeta;
        Y = eml_pwls_qm( Y ,  Y , niter,beta_pwls,sigma,Par.kappa);
        Y = (Y/Par.kappa);
        Dy               = diff_ND(Y,2);
    else
        Y                = max(log(I0./SignoP),0);
        Dy               = diff_ND(Y,2);
    end;
end

I                = max(round(SignoP),1);
normD            = sum(Dy(:).^2);
Q                = zeros(size(Dy));
Z                = Q;
% wbeta            = ones(size(Dy))*beta;
normD            = sqrt(normD);
var              = 11;%sum((SignoP(:)-I(:)).^2)/N*100 ;
fftnD            = zeros(sizeS);
for i = 1:ndim
    fftnD     = fftnD +permute(( abs(psf2otf([+1; -1], sizeS([i:end,1:i-1])))).^2,[ndim-i+2:ndim,1:ndim-i+1]);
end
subtol = 0.3e-2;

%% Main loop
for i = 1:MaxIter
    
    %% Update I
    I   = UpdateI(I,Y,I0,SignoP,var,sizeS);
    %% Update Z
    Z  = (mu*Dy+Q)./(beta+mu);                  
    %% Update Y
    Y       = UpdateY(Y,I,I0,mu,Z,Q,fftnD,sizeS,subtol);
    Dy      = diff_ND(Y,2);
    %% Update beta
    beta = N/sum(Dy(:).^2);
    %% Update mu & Q
    temp    = Dy-Z;
    Q    = Q+mu*(temp);
    reErr   =  sum(temp(:).^2)/normD;
    mu  = mu*rho;

    subtol  = reErr/10;
    if reErr<tol
        disp('convergence');
        break
    end
end
Y   = gather(Y);
I   = gather(I);
var = gather(var);
Z   = gather(Z);
if i==MaxIter
    disp(['not convergence£¬relative error is ', num2str(reErr)])
end


end

function I   = UpdateI(I,Y,I0,SignoP,var,sizeS)
Jugefun  = @(X,Ind) -Y(Ind)+log(I0(Ind))-log(X(Ind))+(2*SignoP(Ind)-2*X(Ind)+1)/2/var;
stoption = 0;
BadInd   = ones(sizeS)>0;
BadInd(BadInd)   = (Jugefun(I,BadInd)<0);
while ~stoption
    I(BadInd)         = I(BadInd)-1;
    BadInd(BadInd)        = (Jugefun(I,BadInd)<0);
    stoption              = (sum(BadInd(:))==0);
end

Jugefun  = @(X,Ind) Y(Ind)-log(I0(Ind))+log(X(Ind)+1)-(2*SignoP(Ind)-2*X(Ind)-1)/2/var;
stoption = 0;
BadInd   = ones(sizeS)>0;
BadInd(BadInd)   = (Jugefun(I,BadInd)<0);
while ~stoption
    I(BadInd)         = I(BadInd)+1;
    BadInd(BadInd)        = (Jugefun(I,BadInd)<0);
    stoption              = (sum(BadInd(:))==0);
end
end

function Y   = UpdateY(Y0,I,I0,mu,Z,Q,fftnD,sizeS,tol)
t        = 1;
stoption = 0;
L        = 1*max(I0(:));
difZ     = diffT_ND(Z-Q/mu,sizeS,2);
Y        = Y0;
X        = Y;
i        = 1;
while ~stoption
    perX     = X;
    perT     = t;
    tempX = L/mu*(Y-( I-I0.*exp(-Y) )/L)+difZ;
    X        = ifft2(fft2(tempX) ./(fftnD + L/mu) );
    X        = max(0,real(X));
    t        = ( 1+sqrt(1+ 4 * perT^2 ) )/2;
    Y        = X + (perT/t).*(X-perX);
    i        = i+1;
    stoption = i>6;
    if i>=300
        break
    end
end
end
