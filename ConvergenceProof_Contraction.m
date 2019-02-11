clear
clc

%foldersource='C:\Users\ctc\Documents\GitHub\NewtonRaphsonHeatBalance\';
foldersource='/Users/Shaun/Documents/GitHub/NewtonRaphsonHeatBalance/';

conductorData=importfileAA(strcat(foldersource,'ConductorInfo.csv'));
%[conductorCount,~]=size(conductorData);
deltainfo=importfile42(strcat(foldersource,'deltainfo.csv'));
delta1info=importfile42(strcat(foldersource,'delta1info.csv'));
rootinfo=importfile42(strcat(foldersource,'rootinfo.csv'));
psolinfo=importfile42(strcat(foldersource,'psolinfo.csv'));
windinfo=importfile42(strcat(foldersource,'windinfo.csv'));
ambtempinfo=importfile42(strcat(foldersource,'ambtempinfo.csv'));
currentinfo=importfile42(strcat(foldersource,'currentinfo.csv'));
[weatherPermutationCount,conductorCount]=size(deltainfo);

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
cinfo=(-1*realmax).*ones(weatherPermutationCount,conductorCount);

for c=1:conductorCount    
    cdata=conductorData(c,:);
    diam=cdata.DiamCompleteCable*0.0254;
    beta=(cdata.ResistanceACHighdegcMeter-cdata.ResistanceACLowdegcMeter)/(cdata.HighTemp-cdata.LowTemp);
    alpha=cdata.ResistanceACHighdegcMeter-beta*cdata.HighTemp;    
    counter3=0;
    for w=1090:weatherPermutationCount  
        counter3=counter3+1;
        GuessTc=((psolinfo(w,c)+(currentinfo(w,c)^2)*(alpha+25*beta))/(pi*diam*sigmab*epsilons)+((ambtempinfo(w,c)+273)^4))^(1/4)-273; 
        counter2=0;
        for Tcc=rootinfo(w,c)-deltainfo(w,c):0.01:GuessTc+delta1info(w,c)-0.01
            counter2=counter2+1;
            [Tc,~,~,~,~,~,~] =GetTempNewtonFirstIteration(currentinfo(w,c),ambtempinfo(w,c),H,diam,phi,windinfo(w,c),alpha,beta,epsilons,psolinfo(w,c),Tcc);
            counter1=0;
            for Tcc1=Tcc+0.01:0.01:GuessTc+delta1info(w,c)
                counter1=counter1+1;
                [Tc1,~,~,~,~,~,~] =GetTempNewtonFirstIteration(currentinfo(w,c),ambtempinfo(w,c),H,diam,phi,windinfo(w,c),alpha,beta,epsilons,psolinfo(w,c),Tcc1);
                C=abs(Tc-Tc1)/abs(Tcc-Tcc1);
                if(C>1 && windinfo(w,c)~=0)
                end
                if(C>cinfo(w,c))
                    cinfo(w,c)=C;
                end
            end
        end
    end
    disp(strcat(num2str(delta(c)),',',num2str(delta1(c)),',',num2str(100*c/conductorCount),',',cellstr(cdata.CodeWord)));
end
%conductorData.delta=delta;
%conductorData.delta1=delta1;
%writetable(conductorData,strcat(foldersource,'ConductorValidationResults.csv')); 
%writetable(rootinfo,strcat(foldersource,'rootinfo.csv')); 

