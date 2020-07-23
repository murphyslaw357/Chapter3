clear
clc
close all

if(ispc==1)
    startCond = 132
    foldersource='C:\Users\ctc\Documents\GitHub\Chapter3\';
    addpath(genpath(foldersource))
    addpath(genpath(strcat(foldersource,'Step1_PolyTrain\')))
elseif(ismac==1)
    foldersource='/Users/Shaun/Documents/GitHub/Chapter3/';
elseif(isunix==1)
    startCond = str2num(getenv('SGE_TASK_ID'))
    endCond = startCond
    foldersource='/mnt/HA/groups/nieburGrp/Shaun/Chapter3/';
end

%% Load conductor info
load(strcat(foldersource,'conductorInfoStep2.mat'))

[conductorCount,~]=size(conductorInfo);
% conductorInfo.ResistanceACLowdegc=conductorInfo.ResistanceDCLowdegc;
% conductorInfo.simulated=zeros(conductorCount,1);
% conductorInfo.Cmax=zeros(conductorCount,1);
% conductorInfo.Cmin=zeros(conductorCount,1);
% conductorInfo.convergeCurrent = zeros(conductorCount,1);
% conductorInfo.lowestRise = realmax.*ones(conductorCount,1);
% conductorInfo.lilBottomEnd = zeros(conductorCount,1);
% conductorInfo.lilTopEnd = zeros(conductorCount,1);
% conductorInfo.guessMAPE = zeros(conductorCount,1);
% conductorInfo.guessMIN = zeros(conductorCount,1);
% conductorInfo.guessMAX = zeros(conductorCount,1);
% conductorInfo.guessMEAN = zeros(conductorCount,1);
% conductorInfo.guessSTD = zeros(conductorCount,1);
%% Setup weather data
epsilons=0.9;
H=0;
phi=90*pi/180;
maxpsol=1050;
alphas=0.9;
spacer=10;

psols=0:maxpsol/spacer:maxpsol;
winds=0:10/spacer:10;
ambtemps=-33:98/spacer:65;
currents=[1.5:-0.01:0.02, 0.019:-0.001:0.002];
inputCombo = allcomb(currents,psols,winds,ambtemps);
currents=inputCombo(:,1);
psols=inputCombo(:,2);
winds=inputCombo(:,3);
ambtemps=inputCombo(:,4);

weatherPermutationCount = size(inputCombo,1);

%% Run conductor simulation 
c=startCond
cdata=conductorInfo(c,:);  
output = ones(weatherPermutationCount,5);
maxcurrent=ceil(cdata.AllowableAmpacity);
diam=cdata.DiamCompleteCable*0.0254;
beta=(cdata.ResistanceACHighdegcMeter-...
    cdata.ResistanceACLowdegcMeter)/(cdata.HighTemp-cdata.LowTemp);
alpha=cdata.ResistanceACHighdegcMeter-beta*cdata.HighTemp;  
polymodel=str2func(cdata.polymodels);
GuessTcs=GetGuessTemp(currents.*maxcurrent,ambtemps,diam,phi,winds,...
    alpha,beta,epsilons,alphas,psols,polymodel); 
output(:,1)=GuessTcs;
I2Rs = zeros(weatherPermutationCount,1);
Prads = zeros(weatherPermutationCount,1);
Pcons = zeros(weatherPermutationCount,1);
murphyTemps=zeros(weatherPermutationCount,1);
morganTemps=zeros(weatherPermutationCount,1);
minRiseThresh=zeros(weatherPermutationCount,1);
tic
for counter=1:weatherPermutationCount
    if(mod(counter,1000)==0)
        toc
        disp(strcat(num2str(counter/weatherPermutationCount),'_',num2str(counter)))
    end
    GuessTc=GuessTcs(counter);        
    GuessTcRise=GuessTc-ambtemps(counter);

    if(GuessTcRise<=conductorInfo(c,:).minGuessRise)
        roott = ambtemps(counter);
        minRiseThresh(counter)=1;
    else
        [roott,~,~,~,~,~,~] = GetTempNewton(currents(counter)*...
            maxcurrent,ambtemps(counter),H,diam,phi,winds(counter),...
            alpha,beta,epsilons,alphas,psols(counter),GuessTc);
    end
    murphyTemps(counter)=roott;
    [roott2,~,~,~,~,~,~,~,~,~,~,~,~] = GetTempMorgan(currents(counter)*maxcurrent,ambtemps(counter),H,diam,phi,winds(counter),alpha,beta,epsilons,alphas,psols(counter));
    morganTemps(counter)=roott2;
    %I2Rs(counter)=I2R;
    %Prads(counter)=Prad;
    %Pcons(counter)=Pcon;
end

save(strcat(foldersource,'murphyVsMorgan_',num2str(startCond),'.mat'))