clear
clc
%parpool('local',16) 
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
spacer=2;
deltas=zeros(conductorCount,1);
delta1s=zeros(conductorCount,1);
Cs=zeros(conductorCount,1);
for c1=1:5:conductorCount
    parfor c=c1:c1+4
        cdata=conductorData(c,:);
        %A=cellstr(cdata.Type);
        maxcurrent=ceil(1.5*cdata.AllowableAmpacity);
        diam=cdata.DiamCompleteCable*0.0254;
        disp(strcat(num2str(100*c/conductorCount),cellstr(cdata.CodeWord)))
        maxpsol=1050*diam*alphas;

        beta=(cdata.ResistanceACHighdegcMeter-cdata.ResistanceACLowdegcMeter)/(cdata.HighTemp-cdata.LowTemp);
        alpha=cdata.ResistanceACHighdegcMeter-beta*cdata.HighTemp;
        a=0:maxpsol/spacer:maxpsol;
        [~,a1]=size(a);
        a=10:(maxcurrent-10)/spacer:maxcurrent;
        [~,a2]=size(a);
        a=-33:65;
        [~,a3]=size(a);
        limits=zeros(a1*a2*a2,3);
        counter=0;
        crash=1;
        delta=0.5;
        delta1=2;
        Cs(c)=-1*realmax;
        while(crash)
        crash=0;
        for psol=0:maxpsol/spacer:maxpsol
            for imagnitude=10:(maxcurrent-10)/spacer:maxcurrent
                IIstar=abs(imagnitude)^2; 
                for ambtemp=-33:65
                    GuessTc=((psol+IIstar*(alpha+25*beta))/(pi*diam*sigmab*epsilons)+((ambtemp+273)^4))^(1/4)-273; 
                    %for Vw=0:0.1:10
                         Vw=10;
                         [root,~,~,~,~,~,~] =GetTempNewtonFullDiagnostic(imagnitude,ambtemp,H,diam,phi,Vw,cdata.ResistanceACHighdegcMeter,cdata.ResistanceACLowdegcMeter, cdata.HighTemp, cdata.LowTemp,epsilons,psol);
                         counter=counter+1;
                         %limits(counter,1)=root-2;
                         %limits(counter,2)=GuessTc+2;
                         %limits(counter,3)=-1*realmax;
                         for Tcc=root-delta:GuessTc+delta1
                            [Tc,I2R,I2Rprime,Prad,Pradprime,Pradprimeprime,Pcon,Pconprime,Pconprimeprime] =GetTempNewtonFullDiagnosticFirstIteration(imagnitude,ambtemp,H,diam,phi,Vw,cdata.ResistanceACHighdegcMeter,cdata.ResistanceACLowdegcMeter, cdata.HighTemp, cdata.LowTemp,epsilons,psol,Tcc);
                            h=I2R+psol-Prad-Pcon;
                            hprime=I2Rprime-Pradprime-Pconprime;
                            hprimeprime=-1*Pradprimeprime-Pconprimeprime;
                            gprime=abs(1-((hprime*hprime-h*hprimeprime)/(hprime^2)));
                            g=Tcc-h/hprime;
                            if(g<root-delta)
                                delta=root-g;
                                crash=1;
                            elseif(g>GuessTc+delta1)
                                delta1=g;
                                crash=1;
                            end
                            if(crash) 
                                break 
                            end
                            if(gprime>Cs(c))
                                Cs(c)=gprime;
                            end
                            if(gprime>1)
                               disp('doh') 
                            end
                         end
                     %end
                     if(crash)
                         break
                     end
                end
                 if(crash)
                     break
                 end
            end
            if(crash)
                break
            end
        end
        end
        deltas(c,1)=delta;
        delta1s(c,2)=delta1;
        disp(delta)
        disp(delta1)
    end

    conductorData.deltas=deltas;
    conductorData.delta1s=delta1s;
    writetable(conductorData,'ConductorValidationResults.csv'); 
end