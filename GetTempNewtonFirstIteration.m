function [GuessTc,I2R,I2Rprime,Prad,PradPrime,PradPrimePrime,Pcon,PconPrime] =GetTempNewtonFirstIteration(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,Psol,GuessTc,fGrPr,fReNu,fNuRe)
    %I - RMS steady-state load current - amps
    %Ta - ambient temperature - degc
    %H - conductor elevation - meters
    %D - conductor diameter - meters
    %phi - angle between the wind direction and conductor axis - radians
    %Vw - Wind velocity - m/s
	%epsilons - conductor emissivity
    %Psol - solar heating - w/m  
    sigmab=5.6697e-8;
    g=9.805;
    IIstar=abs(I)^2;
%     GuessTc=((Psol+IIstar*(alpha+25*beta))/(pi*D*sigmab*epsilons)+((Ta+273)^4))^(1/4)-273;  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RADIATIVE COOLING/HEATING%%%%%%%%%%%
    %Tsky=(0.0552*(Ta+273)^1.5)-273;
    %Tg=Ta+2;
    Prad=pi*D*sigmab*epsilons*(((GuessTc+273)^4)-((Ta+273)^4));
    %Prad=pi*D*sigmab*epsilons*(((GuessTc+273)^4)-0.5*((Tsky+273)^4)-0.5*((Tg+273)^4));
    PradPrime=4*pi*D*sigmab*epsilons*(GuessTc+273)^3;
    PradPrimePrime=12*pi*D*sigmab*epsilons*(GuessTc+273)^2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%JOULE HEATING%%%%%%%%%%%%%%%%%%%%%%%
    I2Rprime=beta*IIstar;
    Resistance=beta*GuessTc+alpha;
    if(Resistance<0)
        Resistance=0;
    end
    I2R=IIstar*Resistance;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CONVECTIVE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%COOLING/HEATING%%%%%%%%%
%     TfilmkPrime=1/2;
    PrPrime=-1.25e-4;
    GuessTfilm=(GuessTc+Ta)/2;
    vf=((1.32e-5)+(9.5e-8)*GuessTfilm)*((1-((6.5e-3)*H)/288.16)^-5.2561);
    vfPrime=(4.75e-8)*((1-((6.5e-3)*H)/288.16)^-5.2561);
    
    Tfilmk=GuessTfilm+273;

    Gr=(g*(D^3)*abs(GuessTc-Ta))/(Tfilmk*(vf^2));    
    Lambdaf=(2.42e-2)+(7.2e-5)*GuessTfilm;
    LambdafPrime=3.6e-5;
    TfilmkPrime=1/2;
    GrPrime=g*(D^3)*(Tfilmk*(vf)-abs(GuessTc-Ta)*(TfilmkPrime*vf+Tfilmk*2*vfPrime))/...
            ((Tfilmk^2)*(vf^3));
    Pr=0.715-(2.5e-4)*GuessTfilm;
    GrPr=Gr*Pr;
    A=0;
    m=0;
    C=0;
    n=0;
    if(GuessTc~=Ta && Vw==0)
        %pure natural convection 
        Nudf=fGrPr(GrPr);
        [NudfPrime,NudfPrimePrime]=differentiate(fGrPr,GrPr); 
        NudfPrime=NudfPrime*(Gr*PrPrime+GrPrime*Pr);
        NudfPrimePrime=
        Pcon=pi*Nudf*Lambdaf*(GuessTc-Ta);
        PconPrime=pi*(Nudf*Lambdaf+(GuessTc-Ta)*(LambdafPrime*Nudf+NudfPrime*Lambdaf));
        PconPrimePrime=pi*(NudfPrime*Lambdaf+Nudf*LambdafPrime+(LambdafPrime*Nudf+NudfPrime*Lambdaf)+(GuessTc-Ta)*(LambdafPrime*NudfPrime+NudfPrimePrime*Lambdaf+Nudf*LambdafPrime));
    elseif(GuessTc~=Ta && Vw ~=0)
    %mixed convection
        Nudf=fGrPr(GrPr);
        NudfPrime=differentiate(fGrPr,GrPr);
        Req=fNuRe(Nudf);
        ReqPrime=differentiate(fNuRe,Nudf)*NudfPrime*(Gr*PrPrime+GrPrime*Pr);
        Re=sin(phi)*Vw*D/vf;
        RePrime=-(sin(phi)*Vw*D*vfPrime)/(vf^2);
        Reeff=sqrt((Re^2)+(Req^2));
        Top=Re*RePrime+Req*ReqPrime;% (GrPrime*Pr+Gr*PrPrime);
        
        ReeffPrime=Top/Reeff;
        [NueffPrime,NueffPrimePrime]=differentiate(fReNu,Reeff);
        Nueff=fReNu(Reeff);
        Pcon=pi*Nueff*Lambdaf*(GuessTc-Ta);
        PconPrime=pi*(Nueff*Lambdaf+(GuessTc-Ta)*(LambdafPrime*Nueff+NueffPrime*ReeffPrime*Lambdaf));
        %PconPrimePrime=pi*NueffPrimePrime*Lambdaf*(GuessTc-Ta)+2*pi*NueffPrime*(LambdafPrime*(GuessTc-Ta)+Lambdaf)+2*pi*Nueff*LambdafPrime;
    elseif(GuessTc==Ta && Vw~=0)
    %pure forced    
        Re=sin(phi)*Vw*D/vf;  
        RePrime=-(sin(phi)*Vw*D*vfPrime)/(vf^2);
        [NuPrime,NuPrimePrime]=differentiate(fReNu,Re);
        Nu=fReNu(Re);
        Pcon=pi*Nu*Lambdaf*(GuessTc-Ta);
        PconPrime=pi*(Nu*Lambdaf+(GuessTc-Ta)*(LambdafPrime*Nu+NuPrime*RePrime*Lambdaf));
        %PconPrimePrime=pi*NuPrimePrime*Lambdaf*(GuessTc-Ta)+2*pi*NuPrime*(LambdafPrime*(GuessTc-Ta)+Lambdaf)+2*pi*Nu*LambdafPrime;
    else
    %no convection at all    
        Pcon=0;
        PconPrime=0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MISMATCH AND UPDATE%%%%%%%%%%%%%%
    Hprime=I2Rprime-PradPrime-PconPrime;
    if(Hprime>0 || isnan(PconPrime) || PconPrime<0 || PradPrime<0)
        msg='error condition';
        disp(msg);
    end
    
    Mismatch=I2R+Psol-Prad-Pcon;
    update=Mismatch/Hprime;
    GuessTc=GuessTc-update;      
end