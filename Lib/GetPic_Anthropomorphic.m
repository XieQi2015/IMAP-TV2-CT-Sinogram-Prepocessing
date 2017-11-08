function GetPic_Anthropomorphic (Y,Y_true,Y_Noise, tempNameMethod, NmasList, enList)
%% get PIC
Par.lineW  = 3;
Par.times  = 2.8;
Par.outX   = 0.00;
local1     = [0.42,0.67];
local2     =  [0.42,0.64];

numMethod = length(tempNameMethod);
lineNum   = floor(numMethod/4)+1;

for i = 1:length(NmasList)
    figure('units','normalized','position',[0.05+(i-1)*0.01,0.482-0.29*lineNum/2-(i-1)*0.02,0.9,0.29*lineNum],'name',...
        ['Recovery Result of ', num2str(NmasList(i)), ' mAs']);
    X{1} = Y_true;
    X{2} = Y_Noise{i};
    k    = 2;
    NameMethod{1} = 'Clean';
    NameMethod{2} = 'Noisy';
    for j = 1:numMethod
        if enList(j)
            k = k+1;
            X{k} = Y{i,j};
            NameMethod{k} = tempNameMethod{j};
        end
    end
    
    for j = 1:length(NameMethod)
        k = k+1;
        Par.outX   = 0.025;
        Par.color  = [1,0.2,0.2];
        Par.Osz    = [0.15,0.17];
        Z         = WindowBig(normalized(max(min(X{j}(130:end,:),55),4)),local1,Par);
        Par.Osz    = [0.15,0.17]/1.05;
        Par.outX   = 0.0;
        Par.color  = [0.2,1,0.2];
        Z         = WindowBig(Z(:,end:-1:1,:),local2,Par);
        Z         = Z(:,end:-1:1,:);
        subplot(lineNum,4,j);
        imshow(Z); title(NameMethod{j}); 
%         imwrite(Z,[saveroad, '\', NameMethod{j},'.png']);
    end
    
end

end









