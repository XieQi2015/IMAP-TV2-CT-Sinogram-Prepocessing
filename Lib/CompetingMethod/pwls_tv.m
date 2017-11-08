function J=pwls_tv(I,y,G,beta,iter,dt,a)
%function J=ndiTV(I,P,y,alpha,G,beta,iter,dt,a) %written by YWZhang

if ~exist('alpha')
   alpha=0;% CS based on TV without prior image
end

if ~exist('P')
   P=0;
end

if ~exist('iter')
   iter=10;
end

if ~exist('a')
   a=10.^(-8);
end

[ny,nx]=size(I); 


for i=1:iter,  %% do iterations
    g=G'*(G*I-y);
    Gg=(G*g);
    tau=sum(g(:).^2)/sum(Gg(:).^2);
    I_f=I(:,[2:nx nx])-I;
    I_b=I-I(:,[1 1:nx-1]);
    I_u=I-I([1 1:ny-1],:);
    I_d=I([2:ny ny],:)-I;
    I_x_1=I(:,[2:nx nx])-I([1 1:ny-1],[2:nx nx]);
    I_y_1=I([2:ny ny],:)-I([2:ny ny],[1 1:ny-1]);
    den1=(a+I_b.^2+I_u.^2).^(0.5);
    den2=(a+I_f.^2+I_x_1.^2).^(0.5);
    den3=(a+I_d.^2+I_y_1.^2).^(0.5);

  if ~exist('P')
    v=(I_b+I_u)./den1-I_f./den2-I_d./den3;
    norm=(sum(sum(v.^2))).^(0.5);
    v=v./norm;
    u=tau*g+beta*dt*v;
    I=I-u;
  else
    Q=I-P;
    Q_f=Q(:,[2:nx nx])-Q;
    Q_b=Q-Q(:,[1 1:nx-1]);
    Q_u=Q-Q([1 1:ny-1],:);
    Q_d=Q([2:ny ny],:)-Q;
    Q_x_1=Q(:,[2:nx nx])-Q([1 1:ny-1],[2:nx nx]);
    Q_y_1=Q([2:ny ny],:)-Q([2:ny ny],[1 1:ny-1]);
    den4=(a+Q_b.^2+Q_u.^2).^(0.5);
    den5=(a+Q_f.^2+Q_x_1.^2).^(0.5);
    den6=(a+Q_d.^2+Q_y_1.^2).^(0.5);
    v1=(I_b+I_u)./den1-I_f./den2-I_d./den3;
    v2=(Q_b+Q_u)./den4-Q_f./den5-Q_d./den6;
    v=alpha*(v2)+(1-alpha)*(v1);
    norm=(sum(sum(v.^2))).^(0.5);
    v=v./norm;
    u=tau*g+beta*dt*v;
    I=I-u;
  end
  
end    
  
J=I;
