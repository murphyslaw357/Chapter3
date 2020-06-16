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

conductorData=importfile46(strcat(foldersource,'conductorDataResults.csv'));
[conductorCount,~] = size(conductorData);

figure('Renderer', 'painters', 'Position', [10 10 500 600]);

cMax=conductorData.Cmax;
cMin=conductorData.Cmin;
subplot(2,1,1);
plot(cMax,'linewidth',0.9)
ylim([0 1])
yL = get(gca,'YLim');
line([68 68],yL,'LineStyle','--','Color','r');
line([132 132],yL,'LineStyle','--','Color','r');
line([161 161],yL,'LineStyle','--','Color','r');
line([175 175],yL,'LineStyle','--','Color','r');
line([230 230],yL,'LineStyle','--','Color','r');
line([302 302],yL,'LineStyle','--','Color','r');
line([339 339],yL,'LineStyle','--','Color','r');

legend('Maximum C')
xlabel('Conductor Index')
ylabel('Convergence Coefficient')
subplot(2,1,2); 
plot(cMin,'linewidth',0.9)
yL = get(gca,'YLim');
line([68 68],yL,'LineStyle','--','Color','r');
line([132 132],yL,'LineStyle','--','Color','r');
line([161 161],yL,'LineStyle','--','Color','r');
line([175 175],yL,'LineStyle','--','Color','r');
line([230 230],yL,'LineStyle','--','Color','r');
line([302 302],yL,'LineStyle','--','Color','r');
line([339 339],yL,'LineStyle','--','Color','r');
legend('Minimum C')
xlabel('Conductor Index')
ylabel('Convergence Coefficient')
set(gca,'FontSize',10)

savefig(strcat(foldersource,'Figure3_2.fig'))
set(gcf, 'Color', 'w');

if(ispc)
    export_fig C:\Users\ctc\Documents\GitHub\Chapter3\Figure3_2.png -m3
elseif(ismac)
    export_fig /Volumes/THESIS/Github/Chapter3/Figure3_2.png -m3
end