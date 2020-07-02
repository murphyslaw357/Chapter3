clear
clc
close all

if(ispc==1)
    foldersource='C:\Users\ctc\Documents\GitHub\Chapter3\';
elseif(ismac==1)
    foldersource='/Users/Shaun/Documents/GitHub/Chapter3/';
elseif(isunix==1)
    foldersource='/mnt/HA/groups/nieburGrp/Shaun/Chapter3/';
end

% conductorDataTotal=[];
% 
% folderStart = 'C:\Users\ctc\Documents\GitHub\Chapter3\May2020';
% myFiles = dir(fullfile(folderStart,'*.mat'));
% for k = 1:length(myFiles)
%     load(strcat(folderStart,'\',myFiles(k).name),'conductorData')
%     conductorDataTotal = [conductorDataTotal;conductorData];
%     disp(k/length(myFiles))
% %   baseFileName = myFiles(k).name;
% %   fullFileName = fullfile(myFolder, baseFileName);
% %   fprintf(1, 'Now reading %s\n', fullFileName);
% %   [wavData, Fs] = wavread(fullFileName);
% %   % all of your actions for filtering and plotting go here
% end
% 
% conductorDataTotal=sortrows(conductorDataTotal,'Index');
load(strcat(foldersource,'conductorInfoStep2.mat'))
[conductorCount,~] = size(conductorInfo);

convergeLimit=conductorInfo.convergeCurrent;
figure('Renderer', 'painters', 'Position', [10 10 700 400]);
plot(convergeLimit.*100,'linewidth',1)
ylim([0 2])
hold on
yL = get(gca,'YLim');
line([68 68],yL,'LineStyle','--','Color','r');
line([132 132],yL,'LineStyle','--','Color','r');
line([161 161],yL,'LineStyle','--','Color','r');
line([175 175],yL,'LineStyle','--','Color','r');
line([230 230],yL,'LineStyle','--','Color','r');
line([302 302],yL,'LineStyle','--','Color','r');
line([339 339],yL,'LineStyle','--','Color','r');

ax = gca;
ax.YGrid = 'on';
ax.GridLineStyle = '-';
xlabel('Conductor Index')
ylabel('Minimum Convergence Current (%)')

set(gcf, 'Color', 'w');
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)