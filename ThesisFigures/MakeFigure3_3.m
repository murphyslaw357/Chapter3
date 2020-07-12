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
load(strcat(foldersource,'conductorInfoStep2.mat'))
[conductorCount,~] = size(conductorInfo);

figure('Renderer', 'painters', 'Position', [10 10 700 400]);

cMax=conductorInfo.Cmax;
cMin=conductorInfo.Cmin;
subplot(2,1,1);
plot(cMax,'linewidth',1)
% ylim([0 1])
yL = get(gca,'YLim');
line([68 68],yL,'LineStyle','--','Color','r');
line([132 132],yL,'LineStyle','--','Color','r');
line([161 161],yL,'LineStyle','--','Color','r');
line([175 175],yL,'LineStyle','--','Color','r');
line([230 230],yL,'LineStyle','--','Color','r');
line([302 302],yL,'LineStyle','--','Color','r');
line([339 339],yL,'LineStyle','--','Color','r');

legend('Maximum c')
xlabel('Conductor Index')
ylabel('Convergence Coefficient')
% yticks([0 0.25 0.5 0.75 1])
ax = gca;
ax.YGrid = 'on';
ax.GridLineStyle = '-';

subplot(2,1,2); 
plot(cMin,'linewidth',1)
yL = get(gca,'YLim');
line([68 68],yL,'LineStyle','--','Color','r');
line([132 132],yL,'LineStyle','--','Color','r');
line([161 161],yL,'LineStyle','--','Color','r');
line([175 175],yL,'LineStyle','--','Color','r');
line([230 230],yL,'LineStyle','--','Color','r');
line([302 302],yL,'LineStyle','--','Color','r');
line([339 339],yL,'LineStyle','--','Color','r');
legend('Minimum c')
xlabel('Conductor Index')
ylabel('Convergence Coefficient')
yticks([0 0.33e-3 0.66e-3 1e-3])
ax = gca;
ax.YGrid = 'on';
ax.GridLineStyle = '-';

set(gcf, 'Color', 'w');
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)

%savefig(strcat(foldersource,'Figure3_3.fig'))
set(gcf, 'Color', 'w');

% if(ispc)
%     export_fig C:\Users\ctc\Documents\GitHub\Chapter3\Figure3_3.png -m3
% elseif(ismac)
%     export_fig /Volumes/THESIS/Github/Chapter3/Figure3_3.png -m3
% end