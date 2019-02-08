clear
clc

%conductorData=importfileAA('C:\Users\ctc\Documents\GitHub\NewtonRaphsonHeatBalance\ConductorInfo.csv');
conductorData=importfileAA('/Users/Shaun/Documents/GitHub/NewtonRaphsonHeatBalance/ConductorInfo.csv');
[conductorCount,~]=size(conductorData);
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
delta=0.5.*ones(conductorCount,1);
delta1=2.*ones(conductorCount,1);
rootinfo=realmax.*ones(a3*(spacer+1)^3,conductorCount+4);

%Cs=zeros(conductorCount,1);
for c=1:conductorCount
    root=realmax.*ones(a3*(spacer+1)^3,1);

    cdata=conductorData(c,:);
    maxcurrent=ceil(1.5*cdata.AllowableAmpacity);
    diam=cdata.DiamCompleteCable*0.0254;
    maxpsol=1050*diam*alphas;

    beta=(cdata.ResistanceACHighdegcMeter-cdata.ResistanceACLowdegcMeter)/(cdata.HighTemp-cdata.LowTemp);
    alpha=cdata.ResistanceACHighdegcMeter-beta*cdata.HighTemp;
%     a=0:maxpsol/spacer:maxpsol;
%     [~,a1]=size(a);
%     a=10:(maxcurrent-10)/spacer:maxcurrent;
%     [~,a2]=size(a);
%     
    rerun=1;
    
    while(rerun)
        counter=0;
        rerun=0;
        for psol=0:maxpsol/spacer:maxpsol
            disp(psol)
            for imagnitude=10:(maxcurrent-10)/spacer:maxcurrent
                IIstar=abs(imagnitude)^2; 
                for ambtemp=-33:65
                    GuessTc=((psol+IIstar*(alpha+25*beta))/(pi*diam*sigmab*epsilons)+((ambtemp+273)^4))^(1/4)-273; 
                    for Vw=0:10/spacer:10
                        counter=counter+1;
                        if(root(counter,1)==realmax)
                            [roott,~,~,~,~,~,~] =GetTempNewtonFullDiagnostic(imagnitude,ambtemp,H,diam,phi,Vw,alpha,beta,epsilons,psol);
                            root(counter,1)=roott;
                        end
                        for Tcc=root(counter,1)-delta(c):0.1:GuessTc+delta1(c)
                            [Tc,I2R,I2Rprime,Prad,Pradprime,~,Pcon,Pconprime,~] =GetTempNewtonFullDiagnosticFirstIteration(imagnitude,ambtemp,H,diam,phi,Vw,alpha,beta,epsilons,psol,Tcc);
                            h=I2R+psol-Prad-Pcon;
                            hprime=I2Rprime-Pradprime-Pconprime;
                            g=Tcc-h/hprime;
                            if(g<root(counter,1)-delta(c))
                                delta(c)=root(counter,1)-g;
                                rerun=1;
                            elseif(g>GuessTc+delta1(c))
                                delta1(c)=g-GuessTc;
                                rerun=1;
                            end
                            if(rerun) 
                                break 
                            end
                         end
                     end
                     if(rerun)
                         break
                     end
                end
                 if(rerun)
                     break
                 end
            end
            if(rerun)
                break
            end
        end
    end
    disp(strcat(num2str(delta(c)),',',num2str(delta1(c)),',',num2str(100*c/conductorCount),',',cellstr(cdata.CodeWord)));
    rootinfo(:,c)=root;
end
conductorData.delta=delta;
conductorData.delta1=delta1;
writetable(conductorData,'ConductorValidationResults.csv'); 
writetable(rootinfo,'rootinfo.csv'); 