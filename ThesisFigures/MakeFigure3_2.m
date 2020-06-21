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

[conductorCount,~] = size(conductorDataTotal);

figure('Renderer', 'painters', 'Position', [10 10 700 400]);
plot(conductorDataTotal.lowestRise,'linewidth',1)

yL = get(gca,'YLim');
line([68 68],yL,'LineStyle','--','Color','r');
line([132 132],yL,'LineStyle','--','Color','r');
line([161 161],yL,'LineStyle','--','Color','r');
line([175 175],yL,'LineStyle','--','Color','r');
line([230 230],yL,'LineStyle','--','Color','r');
line([302 302],yL,'LineStyle','--','Color','r');
line([339 339],yL,'LineStyle','--','Color','r');
xlabel('Conductor Index')
ylabel('Minimum Conductor Temperature Rise (°C)')

ax = gca;
ax.YGrid = 'on';
ax.GridLineStyle = '-';

set(gcf, 'Color', 'w');
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)

% if(ispc)
%     export_fig C:\Users\ctc\Documents\GitHub\Chapter3\Figure3_1.png -m3
% elseif(ismac)
%     export_fig /Volumes/THESIS/Github/Chapter3/Figure3_1.png -m3
% end

%fig = gcf;
%fig.PaperPositionMode = 'auto';
%print(gcf,strcat(foldersource,'Figure3_1.jpg'),'-r1200','-djpeg')

% saveas(gcf,strcat(foldersource,'Figure3_1.jpg'))