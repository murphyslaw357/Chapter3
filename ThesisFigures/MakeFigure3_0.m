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

% conductorData=importfile46(strcat(foldersource,'conductorDataResults.csv'));
load(strcat(foldersource,'conductorInfoStep2_HighCurrentLowResolution.mat'))
[conductorCount,~] = size(conductorInfo);

figure('Renderer', 'painters', 'Position', [10 10 700 400]);
plot(conductorInfo.guessErr_min)
hold on
plot(conductorInfo.guessErr_mean)
plot(conductorInfo.guessErr_max)
plot(conductorInfo.guessErr_std)
yL = get(gca,'YLim');
line([68 68],yL,'LineStyle','--','Color','r');
line([132 132],yL,'LineStyle','--','Color','r');
line([161 161],yL,'LineStyle','--','Color','r');
line([175 175],yL,'LineStyle','--','Color','r');
line([230 230],yL,'LineStyle','--','Color','r');
line([302 302],yL,'LineStyle','--','Color','r');
line([339 339],yL,'LineStyle','--','Color','r');
xlabel('Conductor Index')
ylabel('Error [°C]')
legend('Min','Mean','Max','Std')
set(gcf, 'Color', 'w');
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)

% if(ispc)
%     export_fig C:\Users\ctc\Documents\GitHub\Chapter3\Figure3_3.png -m3
% elseif(ismac)
%     export_fig /Volumes/THESIS/Github/Chapter3/Figure3_3.png -m3
% end