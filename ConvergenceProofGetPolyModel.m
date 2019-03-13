clear
clc
close all

if(ispc==1)
    foldersource='C:\Users\ctc\Documents\GitHub\NewtonRaphsonHeatBalance\';
elseif(ismac==1)
    foldersource='/Users/Shaun/Documents/GitHub/NewtonRaphsonHeatBalance/';
elseif(isunix==1)
    foldersource='/mnt/HA/groups/nieburGrp/Shaun/NewtonRaphsonHeatBalance/';
end

load(strcat(foldersource,'GrPrSpline.mat'))
load(strcat(foldersource,'ReNuSpline.mat'))
load(strcat(foldersource,'NuReSpline.mat'))

conductorData=importfileAB(strcat(foldersource,'conductorData.csv'));
[conductorCount,~]=size(conductorData);

conductorData.ResistanceACLowdegc=conductorData.ResistanceDCLowdegc;
conductorData.ResistanceACLowdegcMeter=conductorData.ResistanceACLowdegc./conductorData.MetersperResistanceInterval;
conductorData.ResistanceACHighdegcMeter=conductorData.ResistanceACHighdegc./conductorData.MetersperResistanceInterval;
conductorData.simulated=zeros(conductorCount,1);

Tref=25;
epsilons=0.9;
H=0;
phi=90*pi/180;
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

% polymodels=importfileAB(strcat(foldersource,'polymodels.csv'));
% [polymodelrow,~]=size(polymodels);
bestpoly=zeros(conductorCount,1);
counter=0;

for imagnitude=0.05:(1.5)/spacer:1.55
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

for c1=1:12:conductorCount
    increment=11;
    if(c1+11>conductorCount)
        increment=conductorCount-c1;
    end
    for c=c1:c1+increment
        disp(c)
%         if(c<=polymodelrow) 
             if(~(conductorData(c,:).polymodels==""))
%                 conductorData.polymodels(c)=polymodels.polymodels(c);
                 continue;
             end
%         end
        cdata=conductorData(c,:);
        
        polymodel=[];
        
        convergeCurrent=0;
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
            [roott,Pj,~,Prad,~,Pcon,~] = GetTempNewton(currents(counter)*maxcurrent,ambtemps(counter),H,diam,phi,winds(counter),alpha,beta,epsilons,psols(counter),f,ff,ffinv,polymodel);
            root(counter,1)=roott;           
        end
        IR2s=(((currents.*maxcurrent).^2)).*(alpha+25*beta);
        x=[(((currents.*maxcurrent).^2)).*alpha,(((currents.*maxcurrent).^2)).*beta,psols.*diam,ambtemps,1./(winds+1),(1./(winds+1)).^2];
        err=realmax;
        for i=5:7
            mdl=MultiPolyRegress(x,root,i);
            if(mdl.CVMAE<err)
                err=mdl.CVMAE;
                bestpoly(c)=i;
                conductorData(c,:).polymodels=func2str(mdl.PolynomialExpression);
            end
        end
        rootinfo(:,c)=root;
        deltainfo(:,c)=delta;
        delta1info(:,c)=delta1;
        cinfo(:,c)=cs;
    end
    csvwrite(strcat(foldersource,'rootinfo.csv'),rootinfo);
    csvwrite(strcat(foldersource,'deltainfo.csv'),deltainfo);
    csvwrite(strcat(foldersource,'delta1info.csv'),delta1info);
    csvwrite(strcat(foldersource,'cinfo.csv'),cinfo);
    csvwrite(strcat(foldersource,'convergeCurrents.csv'),convergeCurrents);
    writetable(conductorData,strcat(foldersource,'conductorData.csv'));
end

csvwrite(strcat(foldersource,'psols.csv'),psols);
csvwrite(strcat(foldersource,'winds.csv'),winds);
csvwrite(strcat(foldersource,'ambtemps.csv'),ambtemps);
csvwrite(strcat(foldersource,'currents.csv'),currents);
csvwrite(strcat(foldersource,'polymodels.csv'),conductorData.polymodels);