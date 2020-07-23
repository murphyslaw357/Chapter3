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

load(strcat(foldersource,'tempSolutionMethodComparison.mat'))
%% Murphy
figure('Renderer', 'painters', 'Position', [10 10 700 275]);
scatter(morganTemps-ambtemps,murphyTemps-morganTemps)
title('Murphy Temp. Error')
xlabel('Conductor Temperature Rise [°C]')
ylabel('Solution Error [°C]')
xlim([0 250])
ylim([-5 5])
yticks([-5 -3.75 -2.5 -1.25 0 1.25 2.5 3.75 5])
ax = gca;
ax.YGrid = 'on';
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)
set(gcf, 'Color', 'w');
%% IEEE
figure('Renderer', 'painters', 'Position', [10 10 700 275]);
scatter(morganTemps-ambtemps,IEEE738Temps-morganTemps)
title('IEEE738 Temp. Error')
xlabel('Conductor Temperature Rise [°C]')
ylabel('Solution Error [°C]')
xlim([0 250])
yticks([-60 -50 -40 -30 -20 -10 0 10])
ax = gca;
ax.YGrid = 'on';
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)
set(gcf, 'Color', 'w');
%% Black
figure('Renderer', 'painters', 'Position', [10 10 700 275]);
scatter(morganTemps-ambtemps,blackTemps-morganTemps)
title('Black Temp. Error')
xlabel('Conductor Temperature Rise [°C]')
ylabel('Solution Error [°C]')
xlim([0 250])
yticks([0 250 500 750 1000 1250 1500 1750 2000])
ax = gca;
ax.YGrid = 'on';
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)
set(gcf, 'Color', 'w');
%% Santos
figure('Renderer', 'painters', 'Position', [10 10 700 275]);
scatter(morganTemps-ambtemps,santosTemps-morganTemps)
title('Santos Temp. Error')
xlabel('Conductor Temperature Rise [°C]')
ylabel('Solution Error [°C]')
xlim([0 250])
ax = gca;
ax.YGrid = 'on';
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)
set(gcf, 'Color', 'w');
%% House
figure('Renderer', 'painters', 'Position', [10 10 700 275]);
scatter(morganTemps-ambtemps,houseTemps-morganTemps)
title('House Temp. Error')
xlabel('Conductor Temperature Rise [°C]')
ylabel('Solution Error [°C]')
xlim([0 250])
yticks([-75 -50 -25 0 25 50 75 100])
ax = gca;
ax.YGrid = 'on';
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)
set(gcf, 'Color', 'w');
%% CIGRE
figure('Renderer', 'painters', 'Position', [10 10 700 275]);
scatter(morganTemps-ambtemps,CIGRETemps-morganTemps)
title('CIGRE Temp. Error')
xlabel('Conductor Temperature Rise [°C]')
ylabel('Solution Error [°C]')
xlim([0 250])
ax = gca;
ax.YGrid = 'on';
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)
set(gcf, 'Color', 'w');