clear
clc
close all

if(ispc==1)
    startCond = 1
    endCond = 132
    foldersource='C:\Users\ctc\Documents\GitHub\Chapter3\';
elseif(ismac==1)
    startCond = 1
    endCond = 132
    foldersource='/Volumes/THESIS/GitHub/Chapter3/';
elseif(isunix==1)
    startCond = str2num(getenv('SGE_TASK_ID'))
    endCond = startCond
    foldersource='/mnt/HA/groups/nieburGrp/Shaun/Chapter3/';
end

%if(~isfile(strcat(foldersource,'Step3_ConvergenceProof/',num2str(startCond),'matlab.mat')))
    %% Load conductor info
    load(strcat(foldersource,'conductorInfoStep2.mat'))

    [conductorCount,~]=size(conductorInfo);
    conductorInfo.ResistanceACLowdegc=conductorInfo.ResistanceDCLowdegc;
    conductorInfo.simulated=zeros(conductorCount,1);
    conductorInfo.Cmax=zeros(conductorCount,1);
    conductorInfo.Cmin=zeros(conductorCount,1);
    conductorInfo.guessMAPE = zeros(conductorCount,1);
    conductorInfo.guessMIN = zeros(conductorCount,1);
    conductorInfo.guessMAX = zeros(conductorCount,1);
    conductorInfo.guessMEAN = zeros(conductorCount,1);
    conductorInfo.guessSTD = zeros(conductorCount,1);
    conductorInfo.fPrimeCheck = zeros(conductorCount,1);
    %% Setup weather data
    epsilons=0.9;
    H=0;
    phi=90*pi/180;
    maxpsol=1050;
    alphas=0.9;
    spacer=30;

    psols=0:maxpsol/spacer:maxpsol;
    winds=0:10/spacer:10;
    ambtemps=-33:98/spacer:65;
    currents=[1.5:-0.005:0.02, 0.019:-0.001:0.002];
    inputCombo = allcomb(currents,psols,winds,ambtemps);
    currents=inputCombo(:,1);
    psols=inputCombo(:,2);
    winds=inputCombo(:,3);
    ambtemps=inputCombo(:,4);
    weatherPermutationCount = size(inputCombo,1);

    %% Run conductor simulation 
    if endCond>conductorCount
        endCond=conductorCount;
    end
    for c=startCond:endCond
        disp(c)
        cdata=conductorInfo(c,:);  
        rootts = ones(weatherPermutationCount,1);
        maxcurrent=ceil(cdata.AllowableAmpacity);
        diam=cdata.DiamCompleteCable*0.0254;
        beta=(cdata.ResistanceACHighdegcMeter-...
            cdata.ResistanceACLowdegcMeter)/(cdata.HighTemp-cdata.LowTemp);
        alpha=cdata.ResistanceACHighdegcMeter-beta*cdata.HighTemp;  
        polymodel=str2func(cdata.polymodels);
        GuessTcs=GetGuessTemp(currents.*maxcurrent,ambtemps,diam,phi,winds,...
            alpha,beta,epsilons,alphas,psols,polymodel); 
        I2Rs = zeros(weatherPermutationCount,1);
        Prads = zeros(weatherPermutationCount,1);
        Pcons = zeros(weatherPermutationCount,1);
        cmax=zeros(weatherPermutationCount,1);
        fullRun=zeros(weatherPermutationCount,1);
        fPrimeCheck=zeros(weatherPermutationCount,1);
        fPrimeAvg=zeros(weatherPermutationCount,1);
        fPrimePrimeAvg=zeros(weatherPermutationCount,1);
        sim=0;
        tic
        for counter=1:weatherPermutationCount
            if(mod(counter,1000)==0)
                toc
                disp(strcat(num2str(counter/weatherPermutationCount),'_',num2str(counter)))
            end
            GuessTc=GuessTcs(counter);        
            GuessTcRise=GuessTc-ambtemps(counter);
            if(GuessTcRise<=conductorInfo(c,:).minGuessRise)
                continue;
            end

            [roott,I2R,~,Prad,~,Pcon,~] = GetTempNewton(currents(counter)*...
                maxcurrent,ambtemps(counter),H,diam,phi,winds(counter),...
                alpha,beta,epsilons,alphas,psols(counter),GuessTc);

            I2Rs(counter)=I2R;
            Prads(counter)=Prad;
            Pcons(counter)=Pcon;
            rootts(counter)=roott;

            lilTopEnd=max(roott,GuessTc)+0.05;
            lilBottomEnd=min(roott,GuessTc);
            bigTopEnd=lilTopEnd+10;
            bigBottomEnd=lilBottomEnd-10;
            rerun=1;
            reruncounter=0;

            while(rerun)
                rerun=0;
                reruncounter=reruncounter+1;

                if(reruncounter>5000)
                    msg='error condition: rerun counter exceeded limit';
                    error(msg)
                end

                searchIncrement = (lilTopEnd-lilBottomEnd)/50;
                temps=[(bigBottomEnd:searchIncrement:bigTopEnd)'; lilTopEnd; lilBottomEnd];
                temps(temps<=ambtemps(counter))=[];
                searchCount=size(temps,1);

                Tc=zeros(searchCount,1);
                I2R=zeros(searchCount,1);
                I2Rprime=zeros(searchCount,1);
                Prad=zeros(searchCount,1);
                PradPrime=zeros(searchCount,1);
                PradPrimePrime=zeros(searchCount,1);
                Pcon=zeros(searchCount,1);
                PconPrime=zeros(searchCount,1);
                PconPrimePrime=zeros(searchCount,1);

                for i =1:searchCount
                    [Tc(i),I2R(i),I2Rprime(i),Prad(i),PradPrime(i),PradPrimePrime(i),Pcon(i),PconPrime(i),...
                        PconPrimePrime(i),~,~,~,~,~,~,~,~,~] =GetTempNewtonFirstIteration(...
                        currents(counter)*maxcurrent,ambtemps(counter),H,diam,phi,...
                        winds(counter),alpha,beta,epsilons,alphas,psols(counter),...
                        temps(i),[]);
                end
                h=I2R+psols(counter)*diam*alphas-Pcon-Prad;
                hprime=I2Rprime-PconPrime-PradPrime;
                hprimeprime=-1*PconPrimePrime-PradPrimePrime;

                bigSearch=horzcat(temps,Tc,abs((h.*hprimeprime)./(hprime.^2)),hprime);

                searchRes=bigSearch(bigSearch(:,1)>=lilBottomEnd & ...
                bigSearch(:,1)<= lilTopEnd,:);
                if(max(searchRes(:,2))>bigTopEnd)
                    bigTopEnd = max(searchRes(:,2))+10;
                    rerun=1;
                end
                if(min(searchRes(:,2))<bigBottomEnd && min(searchRes(:,2))>ambtemps(counter))
                    bigBottomEnd = max([ambtemps(counter),min(searchRes(:,2))-10]);
                    rerun=1;
                end
                if(min(searchRes(:,2))<=ambtemps(counter))
                    error(strcat('Convergence current too low - update to less than Ta: ',num2str(currents(counter))))
                end
                if(max(searchRes(:,2))>lilTopEnd)
                    lilTopEnd=max(searchRes(:,2))+0.05;
                    rerun=1;
                end
                if(min(searchRes(:,2))<lilBottomEnd)
                    lilBottomEnd=min(searchRes(:,2));
                    rerun=1;
                end
                if(max(searchRes(:,3))>1)
                    if(GuessTcRise>conductorInfo(c,:).minGuessRise)
                        conductorInfo(c,:).minGuessRise=GuessTcRise;
                    end
                    disp(strcat('Guess Trise too low - |c|>1: ',num2str(GuessTcRise)))
                    sim=sim+1;
                    break
                end
            end
    %         fPrimeAvg(counter)=mean(hprime);
    %         fPrimePrimeAvg(counter)=mean();
            if(all(searchRes(:,4)<0) || all(searchRes(:,4)>0))
                fPrimeCheck(counter)=1;
            end
            fullRun(counter)=1;
            cmax(counter)=max(searchRes(:,3));
         end
         if(all(fPrimeCheck(fullRun==1)))
             conductorInfo(c,:).fPrimeCheck=1;
         end
         guessErr=GuessTcs(fullRun==1)-rootts(fullRun==1);
         conductorInfo(c,:).guessMIN=min(guessErr);
         conductorInfo(c,:).guessMAX=max(guessErr);
         conductorInfo(c,:).guessMEAN=mean(guessErr);
         conductorInfo(c,:).guessSTD=std(guessErr);
         conductorInfo(c,:).guessMAPE=mean(abs((guessErr)./(rootts(fullRun==1)+273)));

         conductorInfo(c,:).Cmax=max(cmax(fullRun==1));
         conductorInfo(c,:).Cmin=min(cmax(fullRun==1));
         conductorInfo(c,:).simulated=sim;
    end

    cdata=conductorInfo(c,:);  
    save(strcat(foldersource,'Step3_ConvergenceProof/',num2str(startCond),'matlab.mat'),'cdata')
%end