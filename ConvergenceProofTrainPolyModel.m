clear
clc
close all

% numWorkers=feature('numcores')
% parpool('local',numWorkers)

if(ispc==1)
    foldersource='C:\Users\ctc\Documents\GitHub\Chapter3\';
elseif(ismac==1)
    foldersource='/Users/Shaun/Documents/GitHub/Chapter3/';
elseif(isunix==1)
    foldersource='/mnt/HA/groups/nieburGrp/Shaun/Chapter3June_2020/';
    addpath(genpath(foldersource));
end

load(strcat(foldersource,'GrPrSpline.mat'))
load(strcat(foldersource,'ReNuSpline.mat'))
load(strcat(foldersource,'NuReSpline.mat'))

conductorData=importfile1(strcat(foldersource,'ConductorInfo.csv'));
[conductorCount,~]=size(conductorData);

conductorData.ResistanceACLowdegc=conductorData.ResistanceDCLowdegc;
conductorData.ResistanceACLowdegcMeter=conductorData.ResistanceACLowdegc./conductorData.MetersperResistanceInterval;
conductorData.ResistanceACHighdegcMeter=conductorData.ResistanceACHighdegc./conductorData.MetersperResistanceInterval;
conductorData.simulated=zeros(conductorCount,1);
conductorData.polymodels = strings(conductorCount,1); 

Tref=25;
epsilons=0.9;
H=0;
phi=pi/2;
sigmab=5.6697e-8;
maxpsol=1050;
alphas=0.9;
spacer=15;
searchIncrement=0.1;
weatherPermutationCount=(spacer+1)^4;

psols=zeros(weatherPermutationCount,1);
winds=zeros(weatherPermutationCount,1);
ambtemps=zeros(weatherPermutationCount,1);
currents=zeros(weatherPermutationCount,1);
hs=zeros(weatherPermutationCount,1);
hprimes=zeros(weatherPermutationCount,1);
pconprimeprimes=zeros(weatherPermutationCount,1);
pradprimeprimes=zeros(weatherPermutationCount,1);
initguess=zeros(weatherPermutationCount,1);
pconas=zeros(weatherPermutationCount,1);
pconbs=zeros(weatherPermutationCount,1);
pconrs=zeros(weatherPermutationCount,1);
pradas=zeros(weatherPermutationCount,1);
pradbs=zeros(weatherPermutationCount,1);
pradrs=zeros(weatherPermutationCount,1);

bestpoly=zeros(conductorCount,1);
counter=0;

for imagnitude=0.01:(0.5)/spacer:0.51
    for psol=0:maxpsol/spacer:maxpsol
        for ambtemp=-33:98/spacer:65
            for Vw=0:10/spacer:10
                counter=counter+1;
                psols(counter)=psol*alphas;
                winds(counter)=Vw;
                ambtemps(counter)=ambtemp;
                currents(counter)=imagnitude;
            end
        end
    end
end

convergeCurrents=zeros(conductorCount,1);
deltainfo=zeros(weatherPermutationCount,conductorCount);
delta1info=zeros(weatherPermutationCount,conductorCount);
rootinfo=zeros(weatherPermutationCount,conductorCount);
cinfo=zeros(weatherPermutationCount,conductorCount);
stepinfo=zeros(weatherPermutationCount,conductorCount);

startCond=231;
endCond=232;
for c=startCond:endCond
    disp(c)
    cdata=conductorData(c,:);

    root=realmax.*ones(weatherPermutationCount,1);
    delta=zeros(weatherPermutationCount,1);
    delta1=zeros(weatherPermutationCount,1);
    prads=zeros(weatherPermutationCount,1);
    pcons=zeros(weatherPermutationCount,1);
    pir2s=zeros(weatherPermutationCount,1);
    cs=zeros(weatherPermutationCount,1);

    maxcurrent=ceil(cdata.AllowableAmpacity);
    diam=cdata.DiamCompleteCable*0.0254;
    beta=(cdata.ResistanceACHighdegcMeter-cdata.ResistanceACLowdegcMeter)/(cdata.HighTemp-cdata.LowTemp);
    alpha=cdata.ResistanceACHighdegcMeter-beta*cdata.HighTemp;  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for counter=1:weatherPermutationCount
        [roott,~,~,~] = GetTempIEEEBisection(currents(counter)*maxcurrent,ambtemps(counter),H,diam,phi,winds(counter),alpha,beta,epsilons,psols(counter),f,ff,ffinv);
        root(counter)=roott;           
    end
    x=[(((currents.*maxcurrent).^2)).*alpha,(((currents.*maxcurrent).^2)).*beta,psols.*diam,ambtemps,1./(winds+1),(1./(winds+1)).^2];
    err=realmax;
    i=5;
    %for i=5:8
        mdl=MultiPolyRegress(x,root,i);
        if(mdl.CVMAE<err)
            err=mdl.CVMAE;
            bestpoly(c)=i;
            conductorData(c,:).polymodels=func2str(mdl.PolynomialExpression);
        end
    %end
 end
    %disp(strcat(foldersource,num2str(c1),'matlab.mat'))
    save(strcat('C:\Users\ctc\Documents\GitHub\Chapter3\June2020\',num2str(startCond),'matlab.mat'),'conductorData')
    %writetable(conductorData,strcat(foldersource,'conductorData.csv'));


% writetable(conductorData,strcat(foldersource,'conductorData.csv'));