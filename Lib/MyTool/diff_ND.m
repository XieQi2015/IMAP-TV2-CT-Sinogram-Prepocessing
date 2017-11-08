function diff_X = diff_ND(tenX, ndim)
diff_X = [];
for i = 1:ndim
    dfX1   = diff(tenX,1,i);
    tempX  = permute(tenX,[i:ndim,1:i-1]);
    dfX    = cat(1,permute(dfX1,[i:ndim,1:i-1]),tempX(1,:,:)-tempX(end,:,:)) ;
    diff_X = cat(ndim+1,diff_X,permute(dfX,[ndim-i+2:ndim,1:ndim-i+1]));    
end
end