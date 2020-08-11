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
criteria=morganTemps<=250*1.1;
%% Murphy
figure('Renderer', 'painters', 'Position', [10 10 700 275]);
subplot(1,2,1)
scatter(morganTemps(criteria)-ambtemps(criteria),murphyTemps(criteria)-morganTemps(criteria),'.')
xlabel('Conductor Temperature Rise [°C]')
ylabel('Solution Error [°C]')
xlim([0 350])
ylim([-5 5])
yticks([-5 -3.75 -2.5 -1.25 0 1.25 2.5 3.75 5])
ax = gca;
ax.YGrid = 'on';
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)
subplot(1,2,2)
scatter(winds(criteria),murphyTemps(criteria)-morganTemps(criteria),'.')
xlabel('Wind Speed [m/s]')
ylabel('Solution Error [°C]')
ylim([-5 5])
yticks([-5 -3.75 -2.5 -1.25 0 1.25 2.5 3.75 5])
ax = gca;
ax.YGrid = 'on';
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)
set(gcf, 'Color', 'w');
sgtitle('Ch. 3 Newton-Raphson Method - Solution Error')
%% IEEE
figure('Renderer', 'painters', 'Position', [10 10 700 275]);
subplot(1,2,1)
scatter(morganTemps(criteria)-ambtemps(criteria),IEEE738Temps(criteria)-morganTemps(criteria),'.')
xlabel('Conductor Temperature Rise [°C]')
ylabel('Solution Error [°C]')
xlim([0 350])
yticks([-60 -50 -40 -30 -20 -10 0 10])
ax = gca;
ax.YGrid = 'on';
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)
subplot(1,2,2)
scatter(winds(criteria),IEEE738Temps(criteria)-morganTemps(criteria),'.')
xlabel('Wind Speed [m/s]')
ylabel('Solution Error [°C]')
yticks([-60 -50 -40 -30 -20 -10 0 10])
ax = gca;
ax.YGrid = 'on';
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)
set(gcf, 'Color', 'w');
sgtitle('IEEE738 Method - Solution Error')
%% Black
figure('Renderer', 'painters', 'Position', [10 10 700 275]);
subplot(1,2,1)
%scatter(morganTemps(criteria==1)-ambtemps(criteria==1),blackTemps(criteria==1)-morganTemps(criteria==1),'.')
%scatter(morganTemps(criteria==1 & winds~=1)-ambtemps(criteria==1 & winds~=1),blackTemps(criteria==1 & winds~=1)-morganTemps(criteria==1 & winds~=1),'.')
criteria2 = criteria==1 & (winds<0.5 | winds>1.5);
scatter(morganTemps(criteria2)-ambtemps(criteria2),blackTemps(criteria2)-morganTemps(criteria2),'.')

xlabel('Conductor Temperature Rise [°C]')
ylabel('Solution Error [°C]')
% xlim([0 350])
yticks([0 250 500 750 1000 1250 1500 1750 2000])
ax = gca;
ax.YGrid = 'on';
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)
subplot(1,2,2)
%scatter(winds(criteria==1),blackTemps(criteria==1)-morganTemps(criteria==1),'.')
scatter(winds(criteria2),blackTemps(criteria2)-morganTemps(criteria2),'.')
xlabel('Wind Speed [m/s]')
ylabel('Solution Error [°C]')
yticks([0 250 500 750 1000 1250 1500 1750 2000])
ax = gca;
ax.YGrid = 'on';
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)
set(gcf, 'Color', 'w');
sgtitle('Black Method - Solution Error')
%% Santos
figure('Renderer', 'painters', 'Position', [10 10 700 275]);
subplot(1,2,1)
scatter(morganTemps(criteria)-ambtemps(criteria),santosTemps(criteria)-morganTemps(criteria),'.')
xlabel('Conductor Temperature Rise [°C]')
ylabel('Solution Error [°C]')
xlim([0 350])
ax = gca;
ax.YGrid = 'on';
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)
subplot(1,2,2)
scatter(winds(criteria),santosTemps(criteria)-morganTemps(criteria),'.')
xlabel('Wind Speed [m/s]')
ylabel('Solution Error [°C]')
ax = gca;
ax.YGrid = 'on';
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)
set(gcf, 'Color', 'w');
sgtitle('Santos Method - Solution Error')
%% House
figure('Renderer', 'painters', 'Position', [10 10 700 275]);
subplot(1,2,1)
scatter(morganTemps(criteria)-ambtemps(criteria),houseTemps(criteria)-morganTemps(criteria),'.')
xlabel('Conductor Temperature Rise [°C]')
ylabel('Solution Error [°C]')
xlim([0 350])
yticks([-75 -50 -25 0 25 50 75 100])
ax = gca;
ax.YGrid = 'on';
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)
subplot(1,2,2)
scatter(winds(criteria),houseTemps(criteria)-morganTemps(criteria),'.')
xlabel('Wind Speed [m/s]')
ylabel('Solution Error [°C]')
yticks([-75 -50 -25 0 25 50 75 100])
ax = gca;
ax.YGrid = 'on';
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)
set(gcf, 'Color', 'w');
sgtitle('House Method - Solution Error')
%% CIGRE
figure('Renderer', 'painters', 'Position', [10 10 700 275]);
subplot(1,2,1)
scatter(morganTemps(criteria)-ambtemps(criteria),CIGRETemps(criteria)-morganTemps(criteria),'.')
xlabel('Conductor Temperature Rise [°C]')
ylabel('Solution Error [°C]')
xlim([0 350])
ax = gca;
ax.YGrid = 'on';
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)
subplot(1,2,2)
scatter(winds(criteria),CIGRETemps(criteria)-morganTemps(criteria),'.')
xlabel('Wind Speed [m/s]')
ylabel('Solution Error [°C]')
ax = gca;
ax.YGrid = 'on';
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)
set(gcf, 'Color', 'w');
sgtitle('CIGRE Method - Solution Error')