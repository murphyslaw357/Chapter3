if(ispc==1)
    splinesource='C:\Users\ctc\Documents\GitHub\Chapter3\';
    foldersource='C:\Users\ctc\Documents\GitHub\Chapter3\Step1_PolyTrain\';
    startCond = 215
    endCond = 215
elseif(ismac==1)
    foldersource='/Users/Shaun/Documents/GitHub/Chapter3/';
elseif(isunix==1)
    foldersource='/mnt/HA/groups/nieburGrp/Shaun/Chapter3/Step1_PolyTrain/';
    startCond = str2num(getenv('SGE_TASK_ID'))
    if(isempty(startCond))
        startCond=1
        endCond=415
    else
        endCond = startCond
    end
end
restoredefaultpath
addpath(genpath(foldersource));
    
conductorInfo=importfile1(strcat(foldersource,'ConductorInfo.csv'));
[conductorCount,~]=size(conductorInfo);

conductorInfo.ResistanceACLowdegc=conductorInfo.ResistanceDCLowdegc;
conductorInfo.ResistanceACLowdegcMeter=conductorInfo.ResistanceACLowdegc./conductorInfo.MetersperResistanceInterval;
conductorInfo.ResistanceACHighdegcMeter=conductorInfo.ResistanceACHighdegc./conductorInfo.MetersperResistanceInterval;
conductorInfo.simulated=zeros(conductorCount,1);
conductorInfo.polymodels = strings(conductorCount,1); 
conductorInfo.polyorder = zeros(conductorCount,1); 

epsilons=0.9;
H=0;
phi=pi/2;
maxpsol=1050;
alphas=0.9;
spacer=15;

psols=0:maxpsol/spacer:maxpsol;
winds=0:10/spacer:10;
ambtemps=-33:98/spacer:65;
currents=1.502:-0.01:0.002;
inputCombo = allcomb(currents,psols,winds,ambtemps);
currents=inputCombo(:,1);
psols=inputCombo(:,2);
winds=inputCombo(:,3);
ambtemps=inputCombo(:,4);
weatherPermutationCount = size(inputCombo,1);
fitData=zeros(endCond,39);

if endCond>conductorCount
    endCond=conductorCount;
end
for c=startCond:endCond
    disp(c)
    cdata=conductorInfo(c,:);
    maxcurrent=ceil(cdata.AllowableAmpacity);
    diam=cdata.DiamCompleteCable*0.0254;
    beta=(cdata.ResistanceACHighdegcMeter-cdata.ResistanceACLowdegcMeter)/(cdata.HighTemp-cdata.LowTemp);
    alpha=cdata.ResistanceACHighdegcMeter-beta*cdata.HighTemp;  
    rootts=zeros(weatherPermutationCount,1);
    Pcons=zeros(weatherPermutationCount,1);
    Prads=zeros(weatherPermutationCount,1);
    I2Rs=zeros(weatherPermutationCount,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for counter=1:weatherPermutationCount
        [roott,I2R,Prad,Pcon,~,~,~,~,~,~,~,~,~] = GetTempMorgan(currents(counter).*maxcurrent,ambtemps(counter),H,diam,phi,winds(counter),alpha,beta,epsilons,alphas,psols(counter));
        rootts(counter)=roott;     
        Pcons(counter)=Pcon;
        Prads(counter)=Prad;
        I2Rs(counter)=I2R;
    %    if root(counter)==ambtemps(counter)
    %        error('error with bisection method: Tc = Ta')
    %    end
    end
    x=[(((currents.*maxcurrent).^2)).*alpha,(((currents.*maxcurrent).^2)).*beta,psols.*diam.*alphas,ambtemps,1./(winds+1),(1./(winds+1)).^2];
    err=realmax;
    %i=5;
    for i=1:6
        mdl=MultiPolyRegress(x,rootts-ambtemps,i);
        if(mdl.MAE<err)
            err=mdl.MAE;
            conductorInfo(c,:).polyorder=i;
            conductorInfo(c,:).polymodels=func2str(mdl.PolynomialExpression);
        end
        fitData(c,(i-1)*5+1)=min(mdl.Residuals);
        fitData(c,(i-1)*5+2)=max(mdl.Residuals);
        fitData(c,(i-1)*5+3)=mean(mdl.Residuals);
        fitData(c,(i-1)*5+4)=std(mdl.Residuals);
        fitData(c,(i-1)*5+5)=mdl.CVMAE;
    end
    fitData(c,36)=min(rootts);
    fitData(c,37)=max(rootts);
    fitData(c,38)=mean(rootts);
    fitData(c,39)=std(rootts);
end
conductorInfo = conductorInfo(startCond:endCond,:);
conductorInfo.Index = (startCond:endCond)';
disp(strcat(foldersource,num2str(startCond),'matlab.mat'))
save(strcat(foldersource,num2str(startCond),'matlab.mat'))