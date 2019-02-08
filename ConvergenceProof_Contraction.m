clear
clc

conductorData=importfileAA('C:\Users\ctc\Documents\GitHub\NewtonRaphsonHeatBalance\ConductorInfo.csv');
%conductorData=importfileAA('/Users/Shaun/Documents/GitHub/NewtonRaphsonHeatBalance/ConductorInfo.csv');
[conductorCount,~]=size(conductorData);
convergenceData=importfile('C:\Users\ctc\Documents\GitHub\NewtonRaphsonHeatBalance\ConductorValidationResults.csv');

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
delta=0.5.*ones(weatherPermutationCount,1);
delta1=2.*ones(weatherPermutationCount,1);
rootinfo=realmax.*ones(weatherPermutationCount,conductorCount+4);

for c=1:conductorCount
    root=realmax.*ones(a3*(spacer+1)^3,1);
    cdata=conductorData(c,:);
    maxcurrent=ceil(1.5*cdata.AllowableAmpacity);
    diam=cdata.DiamCompleteCable*0.0254;
    maxpsol=1050*diam*alphas;

    beta=(cdata.ResistanceACHighdegcMeter-cdata.ResistanceACLowdegcMeter)/(cdata.HighTemp-cdata.LowTemp);
    alpha=cdata.ResistanceACHighdegcMeter-beta*cdata.HighTemp;    
    counter=0;
    for psol=0:maxpsol/spacer:maxpsol
        disp(psol)
        for imagnitude=10:(maxcurrent-10)/spacer:maxcurrent
            IIstar=abs(imagnitude)^2; 
            for ambtemp=-33:65
                GuessTc=((psol+IIstar*(alpha+25*beta))/(pi*diam*sigmab*epsilons)+((ambtemp+273)^4))^(1/4)-273; 
                for Vw=0:10/spacer:10
                    counter=counter+1;
                    
                        
                      
                        %dist=ceil(10*((GuessTc+delta1(counter))-(root(counter,1)-delta(counter))));
                        %ftracker=zeros(dist,1);
                        %fprimetracker=zeros(dist,1);
                        %pcontracker=zeros(dist,1);
                        %pconprimetracker=zeros(dist,1);
                        %temp=zeros(dist,1); 
                        %fcounter=0;

                        for Tcc=root(counter,1)-conductorData.delta(counter):0.1:GuessTc+conductorData.delta1(counter)
                            [Tc,I2R,I2Rprime,Prad,Pradprime,~,Pcon,Pconprime,~] =GetTempNewtonFullDiagnosticFirstIteration(imagnitude,ambtemp,H,diam,phi,Vw,alpha,beta,epsilons,psol,Tcc);
                            h=I2R+psol-Prad-Pcon;
                            hprime=I2Rprime-Pradprime-Pconprime;
                            updatecc=Tcc-h/hprime;
                            for Tcc1=Tcc:0.1:GuessTc+conductorData.delta1(counter)
                                [Tc1,I2R1,I2Rprime1,Prad1,Pradprime1,~,Pcon1,Pconprime1,~] =GetTempNewtonFullDiagnosticFirstIteration(imagnitude,ambtemp,H,diam,phi,Vw,alpha,beta,epsilons,psol,Tcc1);
                                h1=I2R1+psol-Prad1-Pcon1;
                                hprime1=I2Rprime1-Pradprime1-Pconprime1;
                                updatecc1=Tcc1-h1/hprime1;
                                C=abs(updatecc-updatecc1)/abs(Tcc-Tcc1);
                            end
                        end
                    
                 end
%                      if(rerun)
%                          break
%                      end
            end
%                  if(rerun)
%                      break
%                  end
        end
%             if(rerun)
%                 break
%             end
    end
    %end
    disp(strcat(num2str(delta(c)),',',num2str(delta1(c)),',',num2str(100*c/conductorCount),',',cellstr(cdata.CodeWord)));
    rootinfo(:,c)=root;
end
conductorData.delta=delta;
conductorData.delta1=delta1;
writetable(conductorData,'ConductorValidationResults.csv'); 
writetable(rootinfo,'rootinfo.csv'); 

