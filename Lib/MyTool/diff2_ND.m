function [diff2_X ]= diff2_ND(tenX,sizeX)
N        = prod(sizeX);
ndim     = 2;
diff2_X  = zeros(4*N,1);
for i = 1:ndim
    dfX1   = diff(tenX,1,i);
    tempX  = permute(tenX,[i:ndim,1:i-1]);
    dfX    = cat(1,permute(dfX1,[i:ndim,1:i-1]),tempX(1,:,:)-tempX(end,:,:)) ;
    dfX    = permute(dfX,[ndim-i+2:ndim,1:ndim-i+1]);
    for j = 1:ndim
        dfX1    = diff(dfX,1,j);
        tempX   = permute(dfX,[j:ndim,1:j-1]);
        dfX2    = cat(1,permute(dfX1,[j:ndim,1:j-1]),tempX(1,:,:)-tempX(end,:,:)) ;
        dfX2    = permute(dfX2,[ndim-j+2:ndim,1:ndim-j+1]);        
        k = ((i-1)*ndim+j);
        diff2_X(N*(k-1)+1:N*k)=dfX2(:);
    end
end
end