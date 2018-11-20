clear
clc
 
conductorData=importfile('ConductorInfo.csv');
[conductorCount,~]=size(conductorData);
conductorData.Validated=zeros(conductorCount,1);
conductorData.minfirstderivative=zeros(conductorCount,1);
conductorData.ResistanceACLowdegc=conductorData.ResistanceDCLowdegc;
conductorData.ResistanceACLowdegcMeter=conductorData.ResistanceACLowdegc./conductorData.MetersperResistanceInterval;
conductorData.ResistanceACHighdegcMeter=conductorData.ResistanceACHighdegc./conductorData.MetersperResistanceInterval;

Tref=25;
epsilons=0.9;
H=0;
phi=12;
Vw=2;
sigmab=5.6697e-8;
alphas=0.9;
cond=zeros(conductorCount,1);
minFirstDer=zeros(conductorCount,1);
minSecondDer=zeros(conductorCount,1);
maxCurrent=zeros(conductorCount,1);
counters=zeros(conductorCount,1);
spacer=2;
murphyTemp=zeros(conductorCount,1);
cecchiTemp=zeros(conductorCount,1);
blackTemp=zeros(conductorCount,1);
std738Temp=zeros(conductorCount,1);
for c1=1:83:conductorCount
    for c=c1:c1+82
        cdata=conductorData(c,:);
        A=cellstr(cdata.Type);
        cond(c)=1;
        maxcurrent=ceil(1.5*cdata.AllowableAmpacity);
        diam=cdata.DiamCompleteCable*0.0254;
        disp(strcat(num2str(100*c/conductorCount),cellstr(cdata.CodeWord)))
        minFirstDer(c)=-1*realmax;       
        minSecondDer(c)=-1*realmax;
        maxpsol=1050*diam*alphas;

        beta=(cdata.ResistanceACHighdegcMeter-cdata.ResistanceACLowdegcMeter)/(cdata.HighTemp-cdata.LowTemp);
        alpha=cdata.ResistanceACHighdegcMeter-beta*cdata.HighTemp;
        for psol=0:maxpsol/spacer:maxpsol
            for imagnitude=10:(maxcurrent-10)/spacer:maxcurrent
                IIstar=abs(imagnitude)^2; 
                for ambtemp=-33:65
                    for Vw=0:0.1:10
                         [Tc,~,~,~,~,~,~] =GetTempNewtonFullDiagnostic(imagnitude,ambtemp,H,diam,phi,Vw,cdata.ResistanceACHighdegcMeter,cdata.ResistanceACLowdegcMeter, cdata.HighTemp, cdata.LowTemp,epsilons,psol);
                         murphyTemp(c)=Tc;
                         cecchiTemp(c)=ambtemp;
                         [Tc] = GetTempBlack(imagnitude,ambtemp,diam,phi,Vw,cdata.ResistanceACHighdegcMeter,cdata.ResistanceACLowdegcMeter, cdata.HighTemp, cdata.LowTemp,epsilons,psol);
                         blackTemp(c)=Tc;
                         [Tc] = GetTempStd738(imagnitude,ambtemp,diam,phi,Vw,cdata.ResistanceACHighdegcMeter,cdata.ResistanceACLowdegcMeter, cdata.HighTemp, cdata.LowTemp,epsilons,psol,H);
                         std738Temp(c)= Tc;
                    end
                end
            end
        end
    end
    conductorData.Counters=counters;
    conductorData.Validated=cond;
    conductorData.minfirstderivative=minFirstDer;
    conductorData.minsecondderivative=minSecondDer;
    conductorData.Maxcurrent=maxCurrent;
    writetable(conductorData,'ConductorValidationResults.csv'); 
end
