myfun=@(x,lam) lam.^x./(prod(1:x)).*exp(-lam);
% lam = 0:0.1:20;
% y = myfun(5,lam);
% plot(lam,y);grid;

lama = 10;
clear x yy
for i = 1:50
   x(i) = 0+i*1;
   yy(i) =  myfun(x(i),lama);   
end

guassionfun = @(x,lam) 1/(sqrt(2*pi*lam)).*exp(-(x-lam).^2/lam);

yyy = guassionfun(x,lama);
 plot(x,yy)
hold on
plot(x,yyy,'--r')



plot(x,yy);grid;




 a = (sino_17mas*2e-4-2);
 b = exp(-a+log(1239));
 min(b(:))
 max(b(:))
 
 

 
 
 a = (sino_17mas/ 2294.5);
 b = exp(-a)./sigma_17mas;
 min(b(:))
 max(b(:))