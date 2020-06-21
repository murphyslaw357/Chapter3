clear
clc
close all

folderStart = 'C:\Users\ctc\Documents\GitHub\Chapter3\';
conductorInfo=importfile1(strcat(folderStart,'ConductorInfo.csv'));
conductorCount=size(conductorInfo,1);
conductorInfo.polymodels = strings(conductorCount,1); 

for k = 1:2:conductorCount
    disp(k)
    fileName = strcat(folderStart,'June2020\',num2str(k),'matlab.mat');
    if ~isfile(fileName)
        warn(strcat('Missing file: ',fileName))
    else
        conductorDataTmp = load(fileName,'conductorData');
        if k<conductorCount
            conductorInfo(k:k+1,:).polymodels = conductorDataTmp.conductorData(k:k+1,:).polymodels;
        else
            conductorInfo(k:conductorCount,:).polymodels = conductorDataTmp.conductorData(k:conductorCount,:).polymodels;
        end
    end
    %load(strcat(folderStart,'\',myFiles(k).name),'conductorData')
end
save(strcat(folderStart,'conductorInfoPoly.mat'),'conductorInfo')