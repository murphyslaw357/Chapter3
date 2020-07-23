clear
clc
close all
conductorInfoTotal=[];

rootFolder = 'C:\Users\ctc\Documents\GitHub\Chapter3\';
folderStart = 'C:\Users\ctc\Documents\GitHub\Chapter3\Step1_PolyTrain\';
myFiles = dir(fullfile(folderStart,'*.mat'));
for k = 1:length(myFiles)
    load(strcat(folderStart,'\',myFiles(k).name),'conductorInfo')
    conductorInfoTotal = [conductorInfoTotal;conductorInfo];
    disp(k/length(myFiles))
end

conductorInfoTotal=sortrows(conductorInfoTotal,'Index');
conductorInfo = conductorInfoTotal;
save(strcat(rootFolder,'conductorInfoPoly.mat'),'conductorInfo')