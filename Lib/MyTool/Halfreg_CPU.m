% Half regularization function
function W1 = Halfreg_CPU(W,lambdamu)
    lambdamu = lambdamu*2;
    W1=zeros(size(W));
    Nzb=find(abs(W)-3*2^(1/3)/4.*((lambdamu).^(2/3))>0);% new non zeros b
    phi = acos(lambdamu(Nzb)./8.*((abs(W(Nzb))/3).^(-1.5)));
    W1(Nzb)=real(2/3*W(Nzb).*(1 + cos(2/3*(pi - phi))));

%    W1 = sign(W).*max(abs(W)-lambdamu/2,0); %% Soft


%% Hard
%    W1=zeros(size(W));
%    Nzb=find(abs(W)-sqrt(lambdamu)>0);% new non zeros b
%    W1(Nzb)=W(Nzb);
end