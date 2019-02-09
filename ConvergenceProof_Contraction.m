clear
clc

%foldersource='C:\Users\ctc\Documents\GitHub\NewtonRaphsonHeatBalance\';
foldersource='/Users/Shaun/Documents/GitHub/NewtonRaphsonHeatBalance/';

conductorData=importfileAA(strcat(foldersource,'ConductorInfo.csv'));
[conductorCount,~]=size(conductorData);
deltainfo=importfile42(strcat(foldersource,'deltainfo.csv'));
delta1info=importfile42(strcat(foldersource,'delta1info.csv'));
rootinfo=importfile42(strcat(foldersource,'rootinfo.csv'));

conductorData.ResistanceACLowdegc=conductorData.ResistanceDCLowdegc;
conductorData.ResistanceACLowdegcMeter=conductorData.ResistanceACLowdegc./conductorData.MetersperResistanceInterval;
conductorData.ResistanceACHighdegcMeter=conductorData.ResistanceACHighdegc./conductorData.MetersperResistanceInterval;

Tref=25;
epsilons=0.9;
H=0;
phi=90*pi/180;
sigmab=5.6697e-8;
alphas=0.9;
a=-33:65;
[~,a3]=size(a);
spacer=10;
weatherPermutationCount=a3*(spacer+1)^3;
cinfo=(realmax*-1).*ones(weatherPermutationCount,conductorCount);
for c=1:conductorCount
    
    cdata=conductorData(c,:);
    maxcurrent=ceil(1.5*cdata.AllowableAmpacity);
    diam=cdata.DiamCompleteCable*0.0254;
    maxpsol=1050*diam*alphas;

    beta=(cdata.ResistanceACHighdegcMeter-cdata.ResistanceACLowdegcMeter)/(cdata.HighTemp-cdata.LowTemp);
    alpha=cdata.ResistanceACHighdegcMeter-beta*cdata.HighTemp;    
    weathercounter=0;
    for psol=0:maxpsol/spacer:maxpsol
        disp(psol)
        for imagnitude=10:(maxcurrent-10)/spacer:maxcurrent
            IIstar=abs(imagnitude)^2; 
            for ambtemp=-33:65
                GuessTc=((psol+IIstar*(alpha+25*beta))/(pi*diam*sigmab*epsilons)+((ambtemp+273)^4))^(1/4)-273; 
                for Vw=0:10/spacer:10
                    weathercounter=weathercounter+1;

                        for Tcc=rootinfo(weathercounter,c)-deltainfo(weathercounter,c):0.1:GuessTc+delta1info(weathercounter,c)
                            [Tc,I2R,I2Rprime,Prad,Pradprime,Pcon,Pconprime] =GetTempNewtonFullDiagnosticFirstIteration(imagnitude,ambtemp,H,diam,phi,Vw,alpha,beta,epsilons,psol,Tcc);
                            h=I2R+psol-Prad-Pcon;
                            hprime=I2Rprime-Pradprime-Pconprime;
                            updatecc=Tcc-h/hprime;
                            for Tcc1=Tcc+0.1:0.1:GuessTc+delta1info(weathercounter,c)
                                [Tc1,I2R1,I2Rprime1,Prad1,Pradprime1,Pcon1,Pconprime1] =GetTempNewtonFullDiagnosticFirstIteration(imagnitude,ambtemp,H,diam,phi,Vw,alpha,beta,epsilons,psol,Tcc1);
                                h1=I2R1+psol-Prad1-Pcon1;
                                hprime1=I2Rprime1-Pradprime1-Pconprime1;
                                updatecc1=Tcc1-h1/hprime1;
                                C=abs(updatecc-updatecc1)/abs(Tcc-Tcc1);
                                if(isnan(C))
                                end
                                if(C>cinfo(weathercounter,c))
                                    cinfo(weathercounter,c)=C;
                                end
                            end
                        end              
                end
            end
        end
    end
    %end
    disp(strcat(num2str(delta(c)),',',num2str(delta1(c)),',',num2str(100*c/conductorCount),',',cellstr(cdata.CodeWord)));
    rootinfo(:,c)=root;
end
%conductorData.delta=delta;
%conductorData.delta1=delta1;
%writetable(conductorData,strcat(foldersource,'ConductorValidationResults.csv')); 
%writetable(rootinfo,strcat(foldersource,'rootinfo.csv')); 

