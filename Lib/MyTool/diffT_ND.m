function diffT_F = diffT_ND(ten_F, sizeD, ndim)
N       = prod(sizeD); 
diffT_F = zeros(sizeD);
for i = 1:ndim
    tenX   = reshape(ten_F(1+(i-1)*N:i*N),sizeD);
    dfX1   = diff(tenX,1,i);
    tempX  = permute(tenX,[i:ndim,1:i-1]);
    dfX    = -cat(1,tempX(1,:,:)-tempX(end,:,:),permute(dfX1,[i:ndim,1:i-1])) ;
    dfX    = permute(dfX,[ndim-i+2:ndim,1:ndim-i+1]);
    diffT_F = diffT_F+dfX;
end
end