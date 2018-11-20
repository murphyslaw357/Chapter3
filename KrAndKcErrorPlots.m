clear
clc
close all 
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
cond=zeros(conductorCount,1);
minFirstDer=zeros(conductorCount,1);
maxCurrent=zeros(conductorCount,1);
counters=zeros(conductorCount,1);
c=65;

cdata=conductorData(c,:);
A=cellstr(cdata.Type);
cond(c)=1;
maxcurrent=ceil(1.5*cdata.AllowableAmpacity);
diam=cdata.DiamCompleteCable*0.0254;
disp(strcat(num2str(100*c/conductorCount),cellstr(cdata.CodeWord)))
minFirstDer(c)=-1*realmax;        
psol=0;
kctracker=zeros(ceil(maxcurrent/10),99,101);
krtracker=zeros(ceil(maxcurrent/10),99,101);
kcsantos=zeros(ceil(maxcurrent/10),99,101);
krblack=zeros(ceil(maxcurrent/10),99,101);
kcblack=zeros(ceil(maxcurrent/10),99,101);
currentcounter=0;
for imagnitude=10:10:maxcurrent
    currentcounter=currentcounter+1;
    tempcounter=0; 
    for ambtemp=-33:65
        tempcounter=tempcounter+1; 
        windcounter=0;
         for Vw=0:0.1:10
            windcounter=windcounter+1;
            counters(c)=counters(c)+1;
            [Tc,I2R,I2Rprime,Prad,Pradprime,Pcon,Pconprime] =GetTempNewtonFullDiagnostic(imagnitude,ambtemp,H,diam,phi,Vw,cdata.ResistanceACHighdegcMeter,cdata.ResistanceACLowdegcMeter, cdata.HighTemp, cdata.LowTemp,epsilons,psol);
            hprime=I2Rprime-Pradprime-Pconprime;
            kctracker(currentcounter,tempcounter,windcounter)=Pcon/(Tc-ambtemp);
            kcsantos(currentcounter,tempcounter,windcounter)=2.8636*pi*diam*sqrt(Vw/diam);

            krtracker(currentcounter,tempcounter,windcounter)=Prad/(Tc-ambtemp);
            E=1.38e8+(1.39e6)*ambtemp;
            krblack(currentcounter,tempcounter,windcounter)=pi*diam*sigmab*epsilons*E;
            h=
            kcblack(currentcounter,tempcounter,windcounter)=pi*diam*h;
            if(hprime>minFirstDer(c))
                minFirstDer(c)=hprime;
                %conductorData.minfirstderivative(c)=hprime;
            end
            if(hprime>0 || isnan(Tc))
                cond(c)=0;
            end
         end
     end
end

Vw=0:0.1:10;
ambtemp=-33:65;
krblackmax=squeeze(krblack(203,:,:));
krtrackermax=squeeze(krtracker(203,:,:));
kcsantosmax=squeeze(kcsantos(203,:,:));
kctrackermax=squeeze(kctracker(203,:,:));
surf(Vw,ambtemp,100.*(krblackmax-krtrackermax)./krtrackermax)
zlabel('Kr Error (%)')
ylabel('Ambient Temperature (°C)')
xlabel('Wind Speed (m/s)')
figure
surf(Vw,ambtemp,100.*(kcsantosmax-kctrackermax)./kctrackermax)
zlabel('Kc Error (%)')
ylabel('Ambient Temperature (°C)')
xlabel('Wind Speed (m/s)')