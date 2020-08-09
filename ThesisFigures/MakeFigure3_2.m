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
load(strcat(foldersource,'conductorInfoStep3.mat'))
[conductorCount,~] = size(conductorInfo);

figure('Renderer', 'painters', 'Position', [10 10 700 400]);

cMax=conductorInfo.Cmax;
cMin=conductorInfo.Cmin;
subplot(2,1,1);
scatter(conductorInfo.Index,cMax,50,'.')
ylim([0 0.4])
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
yticks([0 0.1 0.2 0.3 0.4])
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)
ax = gca;
ax.YGrid = 'on';
ax.GridLineStyle = '-';

subplot(2,1,2); 
scatter(conductorInfo.Index,cMin,50,'.')
ylim([0.5e-8 2e-8])
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

% yticks([0 0.75e-8 1.5e-8 2.25e-8 3e-8])
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