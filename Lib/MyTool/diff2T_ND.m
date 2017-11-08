function diff2T_F = diff2T_ND(tenF, sizeD)
ndim     = 2;
diff2T_F = zeros(sizeD);
N        = prod(sizeD);
for i = 1:ndim
diff1T_F = zeros(sizeD);
    for j = 1:ndim
        k = ((i-1)*ndim+j);
        tenX   = reshape(tenF(N*(k-1)+1:N*k),sizeD);
        dfX1   = diff(tenX,1,j);
        tempX  = permute(tenX,[j:ndim,1:j-1]);
        dfX    = -cat(1,tempX(1,:,:)-tempX(end,:,:),permute(dfX1,[j:ndim,1:j-1]));
        dfX    = permute(dfX,[ndim-j+2:ndim,1:ndim-j+1]);
        diff1T_F = diff1T_F+dfX;
    end
    dfX1   = diff(diff1T_F,1,i);
    tempX  = permute(diff1T_F,[i:ndim,1:i-1]);
    dfX    = -cat(1,tempX(1,:,:)-tempX(end,:,:),permute(dfX1,[i:ndim,1:i-1]));
    dfX    = permute(dfX,[ndim-i+2:ndim,1:ndim-i+1]);
    diff2T_F = diff2T_F+dfX;
end
end