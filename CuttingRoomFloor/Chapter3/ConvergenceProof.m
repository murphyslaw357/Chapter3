clear
clc

if(ispc==1)
    delim='\';
    folderstart='E:\Chapter3\';
elseif(ismac==1)
    delim='/';
    folderstart='/Volumes/Thesis/Chapter3/';
elseif(isunix==1)
    delim='/';
    folderstart='/mnt/HA/groups/nieburGrp/Shaun/WeatherData/';
end
    sigmab=5.6697e-8;
    
conductorData=importfile(strcat(folderstart,'ConductorInfo.csv'));
[conductorCount,~]=size(conductorData);
conductorData.Validated=zeros(conductorCount,1);
conductorData.minfirstderivative=zeros(conductorCount,1);
conductorData.Tc0hot=zeros(conductorCount,1);
%conductorData.I2Rhot=zeros(conductorCount,1);
%conductorData.Pradhot=zeros(conductorCount,1);
%conductorData.Pconhot=zeros(conductorCount,1);
conductorData.Tc0cold=zeros(conductorCount,1);
%conductorData.I2Rcold=zeros(conductorCount,1);
%conductorData.Pradcold=zeros(conductorCount,1);
%conductorData.Pconcold=zeros(conductorCount,1);
Tref=25;
epsilons=0.9;
H=0;
phi=12;
Vw=2;
for c=1:conductorCount
    cdata=conductorData(c,:);
    A=cellstr(cdata.Type);
    cont=0;
    if(strcmp(A{1},'ACSR'))
        maxcurrent=ceil(1183.8*(cdata.MetalOD^1.3129));
        cont=1;
    elseif(strcmp(A{1},'ACSS'))
        maxcurrent=ceil(2182.6*(cdata.MetalOD^1.3434));
        cont=1;
    elseif(strcmp(A{1},'ACAR'))
        maxcurrent=ceil(0.3709*(cdata.MetalOD^3)-14.305*(cdata.MetalOD^2)+196.37*cdata.MetalOD+552.5);
        cont=1;
    elseif(strcmp(A{1},'AAAC'))
        maxcurrent=ceil(2.4107*(cdata.MetalOD^2)+53.925*cdata.MetalOD+106.71);
        cont=1;
    elseif(strcmp(A{1},'AAC'))
        maxcurrent=ceil(0.0352*(cdata.MetalOD^3)-2.3859*(cdata.MetalOD^2)+76.563*cdata.MetalOD+100.58);
        cont=1;
    elseif(strcmp(A{1},'ACSSTW'))
        maxcurrent=ceil(0.0453*(cdata.MetalOD^3)-4.3276*(cdata.MetalOD^2)+160.23*cdata.MetalOD+1297);
        cont=1;
    elseif(strcmp(A{1},'AACTW'))
        maxcurrent=ceil(0.0599*(cdata.MetalOD^3)-2.4935*(cdata.MetalOD^2)+71.093*cdata.MetalOD+634.12);
        cont=1;
    elseif(strcmp(A{1},'ACSRTW'))
        maxcurrent=ceil(0.009*(cdata.MetalOD^3)-0.9364*(cdata.MetalOD^2)+44.741*cdata.MetalOD+826.83);
        cont=1;
    end
    
    if(cont)
        diam=cdata.MetalOD*0.0254;
        disp(strcat(num2str(100*c/conductorCount),cellstr(cdata.Name)))
        %Simulate worst-case hot conditions
        %[GuessTc] = GetTc0(maxcurrent,75,diam,cdata.HiResistancemeter,cdata.LowResistancemeter, cdata.HiTempC, cdata.LowTempC,epsilons,1050);
        conductorData.minfirstderivative(c)=-1*realmax;                                                                                                                                                  %I,Ta,H,D,Vw,R_T_high,R_T_low, T_high, T_low,epsilons
        [GuessTc,I2R,I2Rprime,Prad,Pradprime,Pcon,Pconprime] =GetTempNewtonFirstIteration(maxcurrent,75,H,diam,phi,0.1,cdata.HiResistancemeter,cdata.LowResistancemeter, cdata.HiTempC, cdata.LowTempC,epsilons,1050);
        conductorData.Tc0hot(c,1)=GuessTc;
        %conductorData.I2Rhot(c,1)=I2R;
        %conductorData.Pconhot(c,1)=Pcon;
        %conductorData.Pradhot(c,1)=Prad;
        
        %Simulate worst-case cold conditions
        %[GuessTc] = GetTc0(0,-100,diam,cdata.HiResistancemeter,cdata.LowResistancemeter, cdata.HiTempC, cdata.LowTempC,epsilons,0);
        [GuessTc,I2R,I2Rprime,Prad,Pradprime,Pcon,Pconprime] =GetTempNewtonFirstIteration(0,-100,H,diam,phi,10,cdata.HiResistancemeter,cdata.LowResistancemeter, cdata.HiTempC, cdata.LowTempC,epsilons,0);
        conductorData.Tc0cold(c,1)=GuessTc;
        
        %
        psol=0;
         tracker=zeros(maxcurrent+1,176);
         for imagnitude=0:maxcurrent
             %disp(imagnitude)
             imagnitude2=imagnitude^2;
             for ambtemp=-100:75
%                 %for psol=0:1050
%                 psol=0;
%                 %val=ones(100,1);
%                 for wind=1:100
                     [GuessTc,I2R,I2Rprime,Prad,Pradprime,Pcon,Pconprime] =GetTempNewtonFirstIteration(imagnitude,ambtemp,H,diam,phi,1/10,cdata.HiResistancemeter,cdata.LowResistancemeter, cdata.HiTempC, cdata.LowTempC,epsilons,psol);
                     tracker(imagnitude+1,ambtemp+101)=imagnitude2*conductorData.alpha(c,1)+imagnitude2*conductorData.beta(c,1)*25+((pi*diam*sigmab*epsilons)*((ambtemp+273)^3))-...
                         27.9904*(imagnitude^(8/3))*(conductorData.beta(c,1)^(4/3))/((diam*epsilons)^(1/3));
                     
                     firstder=I2Rprime-Pradprime-Pconprime;
                     if(firstder>conductorData.minfirstderivative(c))
                         conductorData.minfirstderivative(c)=firstder;
                     end
%                     
%                     %if(I2Rprime>(Pradprime+Pconprime))
%                     %    val(wind,1)=0;
%                     %    %condcuctorData.Valdiated(c)=0;
%                     %    disp('doh') 
%                     %end
%                 end
%                 %condcuctorData.Valdiated(c)=(conductorData.Valdiated(c)==1&&all(val == 1));
%                 %end
             end
         end
         if(min(min(tracker))>0)
             conductorData.Validated(c)=1;
         end
    end
end

writetable(conductorData,strcat(folderstart,'ConductorValidationResults.csv'));