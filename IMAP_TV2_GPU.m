function  Y = IMAP_TV2_GPU(P,I0,b,Par)
% sloving following CT Sinogram Prepocessing problem.
%   min {    sum_i  {  ln(sigma)+(P_i -Q_i)^2/(2*sigma)-Q_i(ln(I0))+Q_iY_i+ln(Q_i!)+I0e^(-Y_i)  } 
%            +2/b|| D_2Y ||^{1/2}_{1/2} - 2Mln(b)     }
%
% Input arguments:
%   P    ...... the corrupted projection data
%   I0   ...... denotes the incident X-ray intensity
%   b    ...... inital scale parameter of the prior distribution
%   Par  ...... an option structure whose fields are as follows:
%      lambda    ... the compromise parameter in ITS, usually setted in [0.1,10];
%      mu        ... initial mu in ADMM algorithm;
%      rho       ... parameter control the increasing speed of mu
%      MaxIter   ... max iteration step number
%      tol       ... termination threshold
%      initIter  ... iteration step number of initialization      
%      sig       ... variance of electronic noise background
%      PwlsInit  ... wehter use PWLS to initialize
%      pwlsBeta  ... parameter of the initial PWLS method
%      kappa     ... parameter of the initial PWLS method
%      Just      ... control parameter for avoiding over smooth
%      wTV       ... weight of different partial derivatives
% Output arguments:
%   Y     ......  the reconstructed sinogram data
%
% more detail can be found in [1]
%
% [1] Qi Xie, Dong Zeng, Qian Zhao, Deyu Meng*, Zongben Xu, Zhengrong Liang, and Jianhua Ma*
%     Robust Low-dose CT Sinogram Preprocessing via Exploiting Noise-generating Mechanism
%     IEEE Transactions on Medical Imaging, 2017
%
% by Qi Xie, 2017.
%==========================================================================

%% initialze parameter
P      = gpuArray(P);
I0     = gpuArray(I0);
b      = gpuArray(b);
sizeS  = size(P);
N      = prod(sizeS);
ndim   = length(sizeS);
if nargin <3
    b  = 1/1.3;
end
if nargin<4
    MaxIter        = 1000; % max iteration step number
    tol            = 1e-7; % termination threshold
    mu             = 1e-5; % initial mu in ADMM algorithm;
    rho            = 1.05; % parameter control the increasing speed of mu
    initIter       = 40;   % iteration step number of initialization   
    sig            = 11;   % variance of electronic noise background
    wTV            = [1,1,1,1]; % weight of different partial derivatives
    Y              = max(log(I0./P),0); % the corrupted sinogram data
    Dy             = diff2_ND_GPU2(Y,sizeS); % the 2 order TV of sinogram data
    Just           = 8.5;  % control parameter for avoiding over smooth
else
    if isfield(Par,'maxIter')         ; MaxIter = Par.maxIter;          else  MaxIter = 1000;          end;
    if isfield(Par,'tol')             ; tol = Par.tol;                  else  tol = 1e-7;              end;
    if isfield(Par,'mu')              ; mu = Par.mu;                    else  mu = 1;                  end;
    if isfield(Par,'rho')             ; rho = Par.rho;                  else  rho = 1.05;              end;
    if isfield(Par,'initIter')        ; initIter = Par.initIter;        else  initIter = 40;           end;
    if isfield(Par,'wTV')             ; wTV = Par.wTV;                  else  wTV = [1,1,1,1];         end;
    if isfield(Par,'sig')             ; sig = Par.sig;                  else  sig = 11;                end;
    if isfield(Par,'Just')            ; Just = Par.Just;                else  Just = 8.5;              end;
    if isfield(Par,'PwlsInit')
        sigma = gather(1./I0);
        Y     = gather(max(log(I0./P),0))*Par.kappa;
        niter = 20;
        beta_pwls = Par.pwlsBeta;
        Y = eml_pwls_qm( Y ,  Y , niter,beta_pwls,sigma,Par.kappa);
        Y = gpuArray(Y/Par.kappa);
        Dy               = diff2_ND_GPU2(Y,sizeS);
    else
        Y                = max(log(I0./P),0);
        Dy               = diff2_ND_GPU2(Y,sizeS);
    end;
end
tol   = gpuArray(tol);
mu    = gpuArray(mu);
rho   = gpuArray(rho);

Q     = max(round(P),1); % the quanta of rays
normD = sum(Dy.^2);
Lam   = zeros(size(Dy),'gpuArray'); % the Lagrange multiplier in ADMM
Z     = Lam;
wb    = Z;
for k = 1:4
    wb(N*(k-1)+1:N*k) = b*wTV(k);
end
normD            = sqrt(normD);
fftnD            = zeros(sizeS,'gpuArray'); % will be used in uwhen pdating Y
for i = 1:ndim
    fftnD     = fftnD + permute((abs(psf2otf([+1, -1 ; -1, +1], sizeS([i:end,1:i-1])))).^2,[ndim-i+2:ndim,1:ndim-i+1]);
    fftnD     = fftnD + permute((abs(psf2otf([+1; -2 ; +1], sizeS([i:end,1:i-1])))).^2,[ndim-i+2:ndim,1:ndim-i+1]);
end

%% main loop
for i = 1:MaxIter
    
    %% Update Q
    Q   = UpdateQ(Q,Y,I0,P,sig,sizeS);
    
    %% Update Z
    if i<=initIter
        Z  = (mu*Dy+Lam)./(70./wb+mu);  % for first several steps, use convex thresholding for initialization
    else
        Z = Halfreg(Dy+Lam/mu,1./(wb*mu));
    end
    
    %% Update Y
    Y       = UpdateY(Y,Q,I0,mu,Z,Lam,fftnD,sizeS);
    Dy      = diff2_ND_GPU2(Y,sizeS);
    
    %% Update b
    for k = 1:4
        wb(N*(k-1)+1:N*k) = sum(abs(Dy(N*(k-1)+1:N*k)).^.5)/N ;
        wb(N*(k-1)+1:N*k) = wb(N*(k-1)+1:N*k)*Just; % avoiding over smooth
    end
    
    %% Update parameter in ADMM algorithm
    temp    = Dy-Z;
    Lam    = Lam+mu*(temp);
    reErr   =  sum(temp(:).^2)/normD;
    mu  = mu*rho;

    subtol  = reErr/10;
    if reErr<tol
        disp('convergence');
        break
    end
end
Y   = gather(Y);
if i==MaxIter
    disp(['not convergence£¬relative error is ', num2str(reErr)])
end
end

%% script for updating Q
function Q   = UpdateQ(Q,Y,I0,P,sig,sizeS)
% solving following problem
% min_Q {sum_i{ (P_i-Q_i)^2/(2*sig)+Q_i(Y_i-ln(I_0i))+ln(Q_i!) }}

%% Decline
Jugefun  = @(X,Ind) -Y(Ind)+log(I0(Ind))-log(X(Ind))+(2*P(Ind)-2*X(Ind)+1)/2/sig;
stoption = 0;
BadInd   = ones(sizeS,'gpuArray')>0;
BadInd(BadInd)   = (Jugefun(Q,BadInd)<0);
while ~stoption
    Q(BadInd)             = Q(BadInd)-1;
    BadInd(BadInd)        = (Jugefun(Q,BadInd)<0);
    stoption              = (sum(BadInd(:))==0);
end

%% Rise
Jugefun  = @(X,Ind) Y(Ind)-log(I0(Ind))+log(X(Ind)+1)-(2*P(Ind)-2*X(Ind)-1)/2/sig;
stoption = 0;
BadInd   = ones(sizeS,'gpuArray')>0;
BadInd(BadInd)   = (Jugefun(Q,BadInd)<0);
while ~stoption
    Q(BadInd)         = Q(BadInd)+1;
    BadInd(BadInd)        = (Jugefun(Q,BadInd)<0);
    stoption              = (sum(BadInd(:))==0);
end
end

%% script for updating Y
function Y   = UpdateY(Y0,I,I0,mu,Z,Lam,fftnD,sizeS)
% solving following problem framework of the accelerated first order method
% min_Y {mu/2|| DY-Z+Lam/mu||^2_2 + sum_i{ Q_iY_i+I_0i(exp(-Y_i))}}

t        = 1;
stoption = 0;
L        = 1*max(I0(:));
difZ     = diff2T_ND_GPU2(Z-Lam/mu,sizeS);
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
end
end
