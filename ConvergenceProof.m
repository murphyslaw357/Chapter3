clear
clc
%parpool('local',16) 
conductorData=importfile9('ConductorInfo3.csv');
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
intervals=zeros(conductorCount,2);
for c1=1:83:conductorCount
    for c=c1:c1+82
        cdata=conductorData(c,:);
        A=cellstr(cdata.Type);
        cond(c)=1;
        maxcurrent=ceil(1.5*cdata.AllowableAmpacity);
        diam=cdata.DiamCompleteCable*0.0254;
        disp(strcat(num2str(100*c/conductorCount),cellstr(cdata.CodeWord)))
%         minFirstDer(c)=-1*realmax;       
%         minSecondDer(c)=-1*realmax;
        maxpsol=1050*diam*alphas;
        %tracker=zeros(maxcurrent+1,176);

        beta=(cdata.ResistanceACHighdegcMeter-cdata.ResistanceACLowdegcMeter)/(cdata.HighTemp-cdata.LowTemp);
        alpha=cdata.ResistanceACHighdegcMeter-beta*cdata.HighTemp;
        a=0:maxpsol/spacer:maxpsol;
        [~,a1]=size(a);
        a=10:(maxcurrent-10)/spacer:maxcurrent;
        [~,a2]=size(a);
        a=-33:65;
        [~,a3]=size(a);
        a=0:0.1:10;
        %[~,a4]=size(a);
        a4=1;
        limits=zeros(a1*a2*a3*a4,3);
        counter=0;
        crash=1;
        delta=0.5;
        delta1=2;
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
                         %dist=abs(GuessTc-root)+1;
                         counter=counter+1;
                         limits(counter,1)=root-2;
                         limits(counter,2)=GuessTc+2;
                         limits(counter,3)=-1*realmax;
                         for Tcc=root-delta:GuessTc+delta1
                            %counters(c)=counters(c)+1;
                            [Tc,I2R,I2Rprime,Prad,Pradprime,Pradprimeprime,Pcon,Pconprime,Pconprimeprime] =GetTempNewtonFullDiagnosticFirstIteration(imagnitude,ambtemp,H,diam,phi,Vw,cdata.ResistanceACHighdegcMeter,cdata.ResistanceACLowdegcMeter, cdata.HighTemp, cdata.LowTemp,epsilons,psol,Tcc);
                            h=I2R+psol-Prad-Pcon;
                            hprime=I2Rprime-Pradprime-Pconprime;
                            hprimeprime=-1*Pradprimeprime-Pconprimeprime;
                            gprime=abs(1-((hprime*hprime-h*hprimeprime)/(hprime^2)));
                            g=Tcc-h/hprime;
                            if(g<root-delta || g>GuessTc+delta1)
                                delta=root-g;
                                crash=1;
                            elseif(g>GuessTc+2)
                                delta1=g;
                                crash=1;
                            end
                            if(crash) 
                                break 
                            end
                            if(gprime>limits(counter,3))
                                limits(counter,3)=gprime;
                            end
                            if(gprime>1)
                               disp('doh') 
                            end
                            %tracker(imagnitude+1,ambtemp+101)=hprime;
%                             hprimeprime=-Pradprimeprime-Pconprimeprime;
%                             if(hprime>minFirstDer(c))
%                                 minFirstDer(c)=hprime;
%                             end
%                             if(abs(hprimeprime)>minSecondDer(c))
%                                 minSecondDer(c)=abs(hprimeprime);
%                             end
                            if(hprime>0 || isnan(Tc))
                                cond(c)=0;
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
        intervals(c,1)=delta;
        intervals(c,2)=delta1;
        disp(delta)
        disp(delta1)
    end

    conductorData.Counters=counters;
    conductorData.Validated=cond;
    conductorData.minfirstderivative=minFirstDer;
    conductorData.minsecondderivative=minSecondDer;
    conductorData.Maxcurrent=maxCurrent;
    writetable(conductorData,'ConductorValidationResults.csv'); 
end