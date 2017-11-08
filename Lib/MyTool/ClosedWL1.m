function [X,svp]=ClosedWL1(Y,C,oureps)
%加权1范数最小化问题的解析解
%这里的目标函数如下
%         sum(w*|Y_i|)+1/2*||Y-X||_F^2
% 其中w_i =C/(sigmaX_i+oureps),oureps是一个足够小的常数
absY      = abs(Y);
signY     = sign(Y);
temp      = (absY-oureps).^2-4*(C-oureps*absY);
ind       = temp>0;
svp       = sum(ind(:));
absY      = absY.*ind;
absY(ind) = max(absY(ind)-oureps+sqrt(temp(ind)),0)/2;
X         = absY.*signY;
end