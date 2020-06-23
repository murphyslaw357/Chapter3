clear
clc
close all
conductorInfoTotal=[];

rootFolder = 'C:\Users\ctc\Documents\GitHub\Chapter3\';
folderStart = 'C:\Users\ctc\Documents\GitHub\Chapter3\Step2_ConvergenceProof\';
myFiles = dir(fullfile(folderStart,'*.mat'));
for k = 1:length(myFiles)
    load(strcat(folderStart,'\',myFiles(k).name),'conductorInfo')
    conductorInfoTotal = [conductorInfoTotal;conductorInfo];
    disp(k/length(myFiles))
%   baseFileName = myFiles(k).name;
%   fullFileName = fullfile(myFolder, baseFileName);
%   fprintf(1, 'Now reading %s\n', fullFileName);
%   [wavData, Fs] = wavread(fullFileName);
%   % all of your actions for filtering and plotting go here
end

conductorInfoTotal=sortrows(conductorInfoTotal,'Index');
conductorInfo = conductorInfoTotal;
% save(strcat(rootFolder,'conductorInfoPoly.mat'),'conductorInfo')