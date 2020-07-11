clear
clc
close all

if(ispc==1)
    startCond = 132;
    foldersource='C:\Users\ctc\Documents\GitHub\Chapter3\';
    cd(foldersource)
    restoredefaultpath
    addpath(genpath(foldersource))
    addpath(genpath(strcat(foldersource,'Step2_ConvergenceProof/')))
elseif(ismac==1)
    foldersource='/Users/Shaun/Documents/GitHub/Chapter3/';
elseif(isunix==1)
%     startCond = str2num(getenv('SGE_TASK_ID'))
    startCond = 132;
    foldersource='/mnt/HA/groups/nieburGrp/Shaun/Chapter3/';
    cd(foldersource)
    restoredefaultpath
    addpath(genpath(foldersource))
    addpath(genpath(strcat(foldersource,'Step2_ConvergenceProof/')))
end

c=startCond
load(strcat(foldersource,'GrPrSpline.mat'))
load(strcat(foldersource,'ReNuSpline.mat'))
load(strcat(foldersource,'NuReSpline.mat'))
%% Load conductor info
load(strcat(foldersource,'conductorInfoStep2.mat'))
[conductorCount,~]=size(conductorInfo);
%% Setup weather data
epsilons=0.9;
H=0;
phi=90*pi/180;
maxpsol=1050;
alphas=0.9;
spacer=30;

psols=0:maxpsol/spacer:maxpsol;
winds=0:10/spacer:10;
ambtemps=-33:98/spacer:65;
currents=[1.5:-0.01:0.02, 0.019:-0.001:0.002];
currents = currents(currents>conductorInfo(c,:).convergeCurrent);
inputCombo = allcomb(currents,psols,winds,ambtemps);
currents=inputCombo(:,1);
psols=inputCombo(:,2);
winds=inputCombo(:,3);
ambtemps=inputCombo(:,4);
weatherPermutationCount = size(inputCombo,1);

%% Run conductor simulation 
cdata=conductorInfo(c,:);  
output = zeros(weatherPermutationCount,4);
maxcurrent=ceil(cdata.AllowableAmpacity);
diam=cdata.DiamCompleteCable*0.0254;
beta=(cdata.ResistanceACHighdegcMeter-...
    cdata.ResistanceACLowdegcMeter)/(cdata.HighTemp-cdata.LowTemp);
alpha=cdata.ResistanceACHighdegcMeter-beta*cdata.HighTemp;  
polymodel=str2func(cdata.polymodels);
GuessTcs=GetGuessTemp(currents.*maxcurrent,ambtemps,diam,phi,winds,...
    alpha,beta,epsilons,alphas,psols,polymodel); 
% if(any(c==[75,111,113,190,254,257,262,288]))
%     polymodel1=str2func(conductorInfo(c-1,:).polymodels);
%     GuessTcs1=GetGuessTemp(currents.*maxcurrent,ambtemps,diam,phi,winds,...
%         alpha,beta,epsilons,alphas,psols,polymodel1);
%     polymodel2=str2func(conductorInfo(c+1,:).polymodels);
%     GuessTcs2=GetGuessTemp(currents.*maxcurrent,ambtemps,diam,phi,winds,...
%         alpha,beta,epsilons,alphas,psols,polymodel2);
%     GuessTcs=(GuessTcs1+GuessTcs2)./2;
%     disp('avg_fn')
% end
output(:,1)=GuessTcs;
I2Rs = zeros(weatherPermutationCount,1);
Prads = zeros(weatherPermutationCount,1);
Pcons = zeros(weatherPermutationCount,1);

tic
for bigcounter=1:10000:weatherPermutationCount
    if(bigcounter+10000<weatherPermutationCount)
        endcounter=10000;
        toc
        disp(bigcounter/weatherPermutationCount)
    else
        endcounter=weatherPermutationCount-bigcounter;
        toc
        disp(bigcounter/weatherPermutationCount)
    end
    for counter=bigcounter:bigcounter+endcounter
        GuessTc=GuessTcs(counter);

        [roott,I2R,~,Prad,~,Pcon,~] = GetTempNewton(currents(counter)*...
            maxcurrent,ambtemps(counter),H,diam,phi,winds(counter),...
            alpha,beta,epsilons,alphas,psols(counter),f,ff,ffinv,GuessTc);
        I2Rs(counter)=I2R;
        Prads(counter)=Prad;
        Pcons(counter)=Pcon;
        output(counter,2)=roott;
    end
end

% Pcons=Pcons(currents>conductorInfo(c,:).convergeCurrent);
% Prads=Prads(currents>conductorInfo(c,:).convergeCurrent);
% I2Rs=I2Rs(currents>conductorInfo(c,:).convergeCurrent);
% psols=psols(currents>conductorInfo(c,:).convergeCurrent);
% winds=winds(currents>conductorInfo(c,:).convergeCurrent);
% ambtemps=ambtemps(currents>conductorInfo(c,:).convergeCurrent);
% output=output(currents>conductorInfo(c,:).convergeCurrent,:);
% currents=currents(currents>conductorInfo(c,:).convergeCurrent,:);
save(strcat(foldersource,'matlab792020.mat'))
%%Remove tag before flight
% load('C:\Users\ctc\Documents\GitHub\Chapter3\ThermalResistanceComparison\matlab792020.mat')
% %

% 
% Kc=Pcons./(output(:,2)-ambtemps);
% Kr=Prads./(output(:,2)-ambtemps);
% figure
% scatter(winds,Kc)
% figure
% scatter(sqrt(winds),Kc)
% figure
% scatter3(sqrt(winds),ambtemps,Kc)
% figure
% plot(ambtemps,Kr)
% figure
% scatter(sqrt(winds),Kr)
% figure
% scatter(psols,Kr)
% figure
% scatter(currents,Kr)