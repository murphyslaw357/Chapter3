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

load(strcat(foldersource,'conductorInfoStep3.mat'))
[conductorCount,~] = size(conductorInfo);

figure('Renderer', 'painters', 'Position', [10 10 700 400]);


%p = plot([conductorInfo.guessMIN conductorInfo.guessMAX conductorInfo.guessMEAN conductorInfo.guessSTD],'linewidth',1)
% y = [conductorInfo.guessMIN conductorInfo.guessMAX conductorInfo.guessMEAN conductorInfo.guessSTD];
scatter(conductorInfo.Index,conductorInfo.guessMIN,50,'.')
hold on
scatter(conductorInfo.Index,conductorInfo.guessMAX,50,'.')
scatter(conductorInfo.Index,conductorInfo.guessMEAN,50,'.')
scatter(conductorInfo.Index,conductorInfo.guessSTD,50,'.')

yL = get(gca,'YLim');
line([68 68],yL,'LineStyle','--','Color','r');
line([132 132],yL,'LineStyle','--','Color','r');
line([161 161],yL,'LineStyle','--','Color','r');
line([175 175],yL,'LineStyle','--','Color','r');
line([230 230],yL,'LineStyle','--','Color','r');
line([302 302],yL,'LineStyle','--','Color','r');
line([339 339],yL,'LineStyle','--','Color','r');
legend('Min.','Max.','Mean','Std.')
xlabel('Conductor Index')
ylabel('Starting Guess Error [°C]')
% yticks([0 0.025 0.05 0.075 0.1])
ax = gca;
ax.YGrid = 'on';
ax.GridLineStyle = '-';

set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)
set(gcf, 'Color', 'w');