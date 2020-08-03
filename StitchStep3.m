clear
clc
close all
conductorInfoTotal=[];

rootFolder = 'C:\Users\ctc\Documents\GitHub\Chapter3\';
folderStart = 'C:\Users\ctc\Documents\GitHub\Chapter3\Step3_ConvergenceProof\';
% rootFolder = '/Volumes/THESIS/Github/Chapter3/';
% folderStart = '/Volumes/THESIS/Github/Chapter3/Step3_ConvergenceProof/';

myFiles = dir(fullfile(folderStart,'*.mat'));
for k = 1:length(myFiles)
    load(strcat(folderStart,'',myFiles(k).name),'cdata')
    conductorInfoTotal = [conductorInfoTotal;cdata];
    disp(k/length(myFiles))
end

conductorInfoTotal=sortrows(conductorInfoTotal,'Index');
conductorInfo = conductorInfoTotal;
% plot(conductorInfoTotal.Index,conductorInfoTotal.convergeCurrent)
% types = unique(conductorInfoTotal.Type);
% for i = 1:size(types,1)
%     figure
%     subplot(3,1,1)
%     scatter(conductorInfoTotal(conductorInfoTotal.Type==types(i),:).DiamCompleteCable,conductorInfoTotal(conductorInfoTotal.Type==types(i),:).AllowableAmpacity)
%     title(string(types(i)))
%     subplot(3,1,2)
%     scatter(conductorInfoTotal(conductorInfoTotal.Type==types(i),:).Index,conductorInfoTotal(conductorInfoTotal.Type==types(i),:).convergeCurrent)%.*conductorInfoTotal(conductorInfoTotal.Type==types(i),:).AllowableAmpacity
%     subplot(3,1,3)
%     scatter(conductorInfoTotal(conductorInfoTotal.Type==types(i),:).DiamCompleteCable,conductorInfoTotal(conductorInfoTotal.Type==types(i),:).convergeCurrent)%.*conductorInfoTotal(conductorInfoTotal.Type==types(i),:).AllowableAmpacity
% end

save(strcat(rootFolder,'conductorInfoStep3.mat'),'conductorInfo')