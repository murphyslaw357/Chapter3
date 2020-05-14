clear
clc
close all
conductorDataTotal=[];

folderStart = 'C:\Users\ctc\Documents\GitHub\Chapter3\May2020';
myFiles = dir(fullfile(folderStart,'*.mat'));
for k = 1:length(myFiles)
    load(strcat(folderStart,'\',myFiles(k).name),'conductorData')
    conductorDataTotal = [conductorDataTotal;conductorData];
    disp(k/length(myFiles))
%   baseFileName = myFiles(k).name;
%   fullFileName = fullfile(myFolder, baseFileName);
%   fprintf(1, 'Now reading %s\n', fullFileName);
%   [wavData, Fs] = wavread(fullFileName);
%   % all of your actions for filtering and plotting go here
end

conductorDataTotal=sortrows(conductorDataTotal,'Index');
writetable(conductorDataTotal,strcat(folderStart,'\conductorDataTotal.csv'));