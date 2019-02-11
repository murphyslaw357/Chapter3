clear
clc

%foldersource='C:\Users\ctc\Documents\GitHub\NewtonRaphsonHeatBalance\';
foldersource='/Users/Shaun/Documents/GitHub/NewtonRaphsonHeatBalance/';
%foldersource='/mnt/HA/groups/nieburGrp/Shaun/NewtonRaphsonHeatBalance/';

conductorData=importfileAA(strcat(foldersource,'ConductorInfo.csv'));
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
alphas=0.9;
a=-33:65;
[~,a3]=size(a);
spacer=10;
weatherPermutationCount=a3*(spacer+1)^3;
deltainfo=zeros(weatherPermutationCount,conductorCount);
delta1info=zeros(weatherPermutationCount,conductorCount);
rootinfo=zeros(weatherPermutationCount,conductorCount+4);
weatherPermutations=zeros(weatherPermutationCount,4);

        maxpsol=1050*alphas;

        for psol=0:maxpsol/spacer:maxpsol
            disp(psol)
            for imagnitude=0:maxcurrent/spacer:maxcurrent
                for ambtemp=-33:65
                    for Vw=0:10/spacer:10
                        
                        counter=counter+1;
                        
                            weatherPermutations(counter,1)=psol;
                            weatherPermutations(counter,2)=imagnitude;
                            weatherPermutations(counter,3)=ambtemp;
                            weatherPermutations(counter,4)=Vw;
                       
                       
                     end
                end
            end
        end
        