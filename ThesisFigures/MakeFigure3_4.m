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

load(strcat(foldersource,'murphyVsMorgan_132.mat'))
simCount=1:weatherPermutationCount;
figure('Renderer', 'painters', 'Position', [10 10 700 450]);
subplot(2,1,1);
plot(morganTemps-ambtemps)
yticks([0 100 200 300 400])
ax = gca;
ax.YGrid = 'on';
ax.GridLineStyle = '-';
xlabel('Simulation Count')
ylabel('Conductor Temperature Rise [°C]')
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)
subplot(2,1,2); 
err=murphyTemps-morganTemps;
scatter(morganTemps(minRiseThresh==1)-ambtemps(minRiseThresh==1),err(minRiseThresh==1))
hold on
scatter(morganTemps(minRiseThresh==0 & abs(err)>1e-2)-ambtemps(minRiseThresh==0 & abs(err)>1e-2),err(minRiseThresh==0 & abs(err)>1e-2))
scatter(morganTemps(minRiseThresh==0 & abs(err)<1e-2)-ambtemps(minRiseThresh==0 & abs(err)<1e-2),err(minRiseThresh==0 & abs(err)<1e-2))
ylim([-5 1])
yticks([-5 -4 -3 -2 -1 0 1])
xlabel('Conductor Temperature Rise')
ylabel('Solution Error [°C]')
ax = gca;
ax.YGrid = 'on';
ax.GridLineStyle = '-';
legend('Error 1','Error 2','Error 3')
set(gca, 'FontName', 'Calibri')
set(gca,'fontsize', 11)
set(gcf, 'Color', 'w');