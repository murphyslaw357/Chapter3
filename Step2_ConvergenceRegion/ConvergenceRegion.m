clear
clc
close all

if(ispc==1)
    c = 215
    foldersource='C:\Users\ctc\Documents\GitHub\Chapter3\';
elseif(ismac==1)
    c = 132
    foldersource='/Volumes/THESIS/GitHub/Chapter3/';
elseif(isunix==1)
    c = str2num(getenv('SGE_TASK_ID'))
    foldersource='/mnt/HA/groups/nieburGrp/Shaun/Chapter3/';
end

restoredefaultpath
addpath(genpath(strcat(foldersource,'Step2_ConvergenceRegion')))
%% Load conductor info
load(strcat(foldersource,'conductorInfoPoly.mat'))
[conductorCount,~]=size(conductorInfo);
if c>conductorCount
    error('Searching for conductor with index outside library')
else
    cdata=conductorInfo(c,:);  
end
maxcurrent=ceil(cdata.AllowableAmpacity);
%% Setup weather data
epsilons=0.9;
H=0;
phi=90*pi/180;
maxpsol=1050;
alphas=0.9;
spacer=20;

psols=0:maxpsol/spacer:maxpsol;
winds=0:10/spacer:10;
ambtemps=-33:98/spacer:65;
currents=[1.5:-0.01:0.02, 0.019:-0.001:0.001];
inputCombo = allcomb(currents,psols,winds,ambtemps);
currents=inputCombo(:,1);
psols=inputCombo(:,2);
winds=inputCombo(:,3);
ambtemps=inputCombo(:,4);

weatherPermutationCount = size(inputCombo,1);

%% Run conductor simulation 
% output = ones(weatherPermutationCount,3);
output2 = ones(weatherPermutationCount,1);
output3 = ones(weatherPermutationCount,1);
diam=cdata.DiamCompleteCable*0.0254;
beta=(cdata.ResistanceACHighdegcMeter-...
    cdata.ResistanceACLowdegcMeter)/(cdata.HighTemp-cdata.LowTemp);
alpha=cdata.ResistanceACHighdegcMeter-beta*cdata.HighTemp;  
% if(c==215)
%     polymodel=str2func(conductorInfo(c-1,:).polymodels);
% else
%     polymodel=str2func(cdata.polymodels);
% end
GuessTcs=GetGuessTemp(maxcurrent.*currents,ambtemps,diam,phi,winds,...
    alpha,beta,epsilons,alphas,psols,polymodel); 
% if(any(c==[108]))
%      polymodel1=str2func(conductorInfo(c-1,:).polymodels);
%      GuessTcs1=GetGuessTemp(currents,ambtemps,diam,phi,winds,...
%          alpha,beta,epsilons,alphas,psols,polymodel1);
%      polymodel2=str2func(conductorInfo(c+1,:).polymodels);
%      GuessTcs2=GetGuessTemp(currents,ambtemps,diam,phi,winds,...
%          alpha,beta,epsilons,alphas,psols,polymodel2);
%      GuessTcs=(GuessTcs1+GuessTcs2)./2;
%      disp('avg_fn')
%  end
% timeout=zeros(weatherPermutationCount,2);
output(:,1)=GuessTcs;
t1 = tic;
for i=1:weatherPermutationCount
    if(mod(i,10000)==0)
        toc(t1)
        disp(strcat(num2str(i/weatherPermutationCount),'_',num2str(i)))
    end
    if(GuessTcs(i)>ambtemps(i))
        [Tc2,~,~,~,~,~,~,~,~,~,~,~,A,m,Cinv,ninv,C,n] = ...
            GetTempNewtonFirstIteration(maxcurrent*currents(i),ambtemps(i),...
            H,diam,phi,winds(i),alpha,beta,epsilons,alphas,...
            psols(i),GuessTcs(i),[]);
        output2(i)=Tc2;
        if(Tc2>ambtemps(i))
            [Tc3,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = ...
                GetTempNewtonFirstIteration(maxcurrent*currents(i),ambtemps(i),...
                H,diam,phi,winds(i),alpha,beta,epsilons,alphas,...
                psols(i),Tc2,[A,m,Cinv,ninv,C,n]);
            output3(i)=Tc3;
        end
    end
end
criteria=(GuessTcs<=ambtemps)|(output2<=ambtemps)|(output3<=ambtemps);
cdata.minGuessRise=max(GuessTcs(criteria)-ambtemps(criteria))+0.1;
% nonConvergeRegion=horzcat(currents(criteria),ambtemps(criteria),winds(criteria),psols(criteria),GuessTcs(criteria),output2(criteria),output3(criteria));
% r1=max(currents(criteria)-min(currents(criteria)));
% r2=max(ambtemps(criteria)-min(ambtemps(criteria)));
% s=1;
% elp=@(x,y,s) ((x.^2)./((s*r1)^2) + ((y).^2)./((s*r2)^2));
% iterationCount=0;
% while(~all(elp(currents(criteria)-min(currents(criteria)),ambtemps(criteria)-min(ambtemps(criteria)),s)<1))
%     iterationCount=iterationCount+1;
%     s=s+0.01;
%     if(iterationCount>500)
%         error('Runaway looping looking for ellipse parameters!')
%     end    
% end
% forbiddenInputs=inputCombo(elp(currents-min(currents(criteria)),ambtemps-min(ambtemps(criteria)),s)<1,:);
% %scatter(nonConvergeRegion(:,1),nonConvergeRegion(:,2))
save(strcat(foldersource,'Step2_ConvergenceRegion/',num2str(c),'matlab.mat'))