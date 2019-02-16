function [GuessTc,I2R,I2Rprime,Prad,PradPrime,PradPrimePrime,Pcon,PconPrime,PconPrimePrime] =GetTempNewtonFirstIteration(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,Psol,GuessTc,fGrPr,fReNu,fNuRe)
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
    vfPrimePrime=0;
    dtop=TfilmkPrime*vf+Tfilmk*vfPrime-abs(GuessTc-Ta)*(TfilmkPrime*vfPrime+2*TfilmkPrime*vfPrime+2*Tfilmk*vfPrimePrime)-(TfilmkPrime*vf+2*Tfilmk*vfPrime);
    GrPrimePrime=g*(D^3)*(((Tfilmk^2)*(vf^3))*dtop-...
        ((Tfilmk*vf-abs(GuessTc-Ta)*(TfilmkPrime*vf+2*Tfilmk*vfPrime)))*(Tfilmk*(vf^3)+3*(vf^2)*vfPrime*(Tfilmk^2)))/...
        ((Tfilmk^4)*(vf^6));
    if(GuessTc~=Ta && Vw==0)
        %pure natural convection 
        Nudf=fGrPr(GrPr);
        [NudfPrime,NudfPrimePrime]=differentiate(fGrPr,GrPr); 
        NudfPrime=NudfPrime*(Gr*PrPrime+GrPrime*Pr);
        NudfPrimePrime=NudfPrimePrime*(Gr*PrPrime+GrPrime*Pr)+NudfPrime*(GrPrime*PrPrime+GrPrimePrime*Pr+GrPrime*PrPrime);
        Pcon=pi*Nudf*Lambdaf*(GuessTc-Ta);
        PconPrime=pi*(Nudf*Lambdaf+(GuessTc-Ta)*(LambdafPrime*Nudf+NudfPrime*Lambdaf));
        PconPrimePrime=pi*(Nudf*LambdafPrime+NudfPrime*Lambdaf+(LambdafPrime*Nudf+NudfPrime*Lambdaf)+(GuessTc-Ta)*(LambdafPrime*NudfPrime+NudfPrimePrime*Lambdaf+NudfPrime*LambdafPrime));
    elseif(GuessTc~=Ta && Vw ~=0)
    %mixed convection
        Nudf=fGrPr(GrPr);
        NudfPrime=differentiate(fGrPr,GrPr);
        NudfPrime=NudfPrime*(Gr*PrPrime+GrPrime*Pr);
        Req=fNuRe(Nudf);
        ReqPrime=differentiate(fNuRe,Nudf);
        ReqPrime=ReqPrime*NudfPrime;
        Re=sin(phi)*Vw*D/vf;
        RePrime=-(sin(phi)*Vw*D*vfPrime)/(vf^2);
        Reeff=sqrt((Re^2)+(Req^2));
        Top=Re*RePrime+Req*ReqPrime;% (GrPrime*Pr+Gr*PrPrime);
        ReeffPrime=Top/Reeff;
        Nueff=fReNu(Reeff);
        [NueffPrime,NueffPrimePrime]=differentiate(fReNu,Reeff);
        NueffPrime=NueffPrime*ReeffPrime;
        Pcon=pi*Nueff*Lambdaf*(GuessTc-Ta);
        PconPrime=pi*(Nueff*Lambdaf+(GuessTc-Ta)*(LambdafPrime*Nueff+NueffPrime*Lambdaf));
        PconPrimePrime=pi*(Nueff*LambdafPrime+NueffPrime*Lambdaf+LambdafPrime*Nueff+NueffPrime*Lambdaf+(GuessTc-Ta)*(LambdafPrime*NueffPrime+NueffPrimePrime*Lambdaf+NueffPrime*LambdafPrime));
    elseif(GuessTc==Ta && Vw~=0)
    %pure forced    
        Re=sin(phi)*Vw*D/vf;  
        RePrime=-(sin(phi)*Vw*D*vfPrime)/(vf^2);
        [NuPrime,NuPrimePrime]=differentiate(fReNu,Re);
        NuPrime=NuPrime*RePrime;
        Nu=fReNu(Re);
        Pcon=pi*Nu*Lambdaf*(GuessTc-Ta);
        PconPrime=pi*(Nu*Lambdaf+(GuessTc-Ta)*(LambdafPrime*Nu+NuPrime*Lambdaf));
        PconPrimePrime=pi*(Nu*LambdafPrime+NuPrime*Lambdaf+LambdafPrime*Nu+NuPrime*Lambdaf+(GuessTc-Ta)*(LambdafPrime*NuPrime+NuPrimePrime*Lambdaf+NuPrime*LambdafPrime));
    else
    %no convection at all    
        Pcon=0;
        PconPrime=0;
        PconPrimePrime=0;
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