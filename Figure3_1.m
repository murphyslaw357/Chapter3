close all

if(ispc==1)
    foldersource='C:\Users\ctc\Documents\GitHub\Chapter3\';
elseif(ismac==1)
    foldersource='/Users/Shaun/Documents/GitHub/Chapter3/';
elseif(isunix==1)
    foldersource='/mnt/HA/groups/nieburGrp/Shaun/Chapter3/';
end

conductorData=importfileAB(strcat(foldersource,'conductorData.csv'));
[conductorCount,~] = size(conductorData);

diff=0.2/15;
convergeLimit=conductorData.convergeCurrent+diff;
plot(convergeLimit.*100)
ylim([0 7.5])
hold on
yL = get(gca,'YLim');
line([68 68],yL,'LineStyle','--','Color','r');
line([132 132],yL,'LineStyle','--','Color','r');
line([161 161],yL,'LineStyle','--','Color','r');
line([175 175],yL,'LineStyle','--','Color','r');
line([230 230],yL,'LineStyle','--','Color','r');
line([302 302],yL,'LineStyle','--','Color','r');
line([339 339],yL,'LineStyle','--','Color','r');
xlabel('Conductor Index')
ylabel('Minimum Convergence Current (%)')
fig = gcf;
fig.PaperPositionMode = 'auto';
print(gcf,strcat(foldersource,'Figure3_1.jpg'),'-r1200','-djpeg')

% saveas(gcf,strcat(foldersource,'Figure3_1.jpg'))