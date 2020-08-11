clear
clc
close all

if(ispc==1)
    c = 15
    foldersource='C:\Users\ctc\Documents\GitHub\Chapter3\';
elseif(ismac==1)
    c = 132
    foldersource='/Volumes/THESIS/GitHub/Chapter3/';
elseif(isunix==1)
    c = str2num(getenv('SGE_TASK_ID'))
    foldersource='/lustre/scratch/srm329/Chapter3/';
end

restoredefaultpath
addpath(genpath(strcat(foldersource,'Step2_ConvergenceRegion')))
%% Load conductor info
% for c=1:415
    load(strcat(foldersource,'conductorInfoPolyIndv/conductorInfoPoly_',num2str(c),'.mat'),'cdata')
%     cdata=conductorInfo(c,:);  
    maxcurrent=ceil(cdata.AllowableAmpacity);
    %% Setup weather data
    epsilons = 0.9;
    H = 0;
    phi = 90*pi/180;
    maxpsol = 1050;
    alphas = 0.9;
    spacer = 15;

    psols = 0:maxpsol/spacer:maxpsol;
    winds = 0:10/spacer:10;
    ambtemps = -33:98/spacer:65;
    currents=[1.5:-0.005:0.005, 0.00499:-0.00001:0.00001];
    inputCombo = allcomb(currents,psols,winds,ambtemps);
    currents = inputCombo(:,1);
    psols = inputCombo(:,2);
    winds = inputCombo(:,3);
    ambtemps = inputCombo(:,4);
    weatherPermutationCount = size(inputCombo,1);

    %% Run conductor simulation 
    Tc2s = realmax.*ones(weatherPermutationCount,1);
    Tc3s = realmax.*ones(weatherPermutationCount,1);
    diam = cdata.DiamCompleteCable*0.0254;
    beta = (cdata.ResistanceACHighdegcMeter-...
        cdata.ResistanceACLowdegcMeter)/(cdata.HighTemp-cdata.LowTemp);
    alpha = cdata.ResistanceACHighdegcMeter-beta*cdata.HighTemp;  
    polymodel=str2func(cdata.polymodels);
    GuessTcs = GetGuessTemp(maxcurrent.*currents,ambtemps,diam,phi,...
        winds,alpha,beta,alphas,psols,polymodel); 

    t1 = tic;
    for i=3255001:weatherPermutationCount
        if(mod(i,10000)==0)
            toc(t1)
            disp(strcat(num2str(i/weatherPermutationCount),'_',num2str(i)))
        end
        if(GuessTcs(i)>ambtemps(i))
            [Tc2,~,~,~,~,~,~,~,~,~,~,~,A,m,Cinv,ninv,C,n] = ...
                GetTempNewtonFirstIteration(maxcurrent*currents(i),...
                ambtemps(i),H,diam,phi,winds(i),alpha,beta,epsilons,...
                alphas,psols(i),GuessTcs(i),[]);
            Tc2s(i)=Tc2;
            if(Tc2>ambtemps(i))
                [Tc3,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = ...
                    GetTempNewtonFirstIteration(maxcurrent*currents(i),...
                    ambtemps(i),H,diam,phi,winds(i),alpha,beta,epsilons,...
                    alphas,psols(i),Tc2,[A,m,Cinv,ninv,C,n]);
                Tc3s(i)=Tc3;
            end
        end
    end

    criteria=(GuessTcs<=ambtemps)|(Tc2s<=ambtemps)|(Tc3s<=ambtemps);
    cdata.minGuessRise=1.25*max(GuessTcs(criteria)-ambtemps(criteria));
    disp(cdata.minGuessRise)
    save(strcat(foldersource,'Step2_ConvergenceRegion/',num2str(c),'matlab.mat'),'cdata')
% end