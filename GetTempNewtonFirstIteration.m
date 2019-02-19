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
        [NudfPrimedgrpr,NudfPrimePrimedgrpr2]=differentiate(fGrPr,GrPr); 
        NudfPrimedtc=NudfPrimedgrpr*(Gr*PrPrime+GrPrime*Pr);
        NudfPrimePrimedtc2=NudfPrimePrimedgrpr2*(Gr*PrPrime+GrPrime*Pr)+NudfPrimedgrpr*(GrPrime*PrPrime+GrPrimePrime*Pr+GrPrime*PrPrime);
        Pcon=pi*Nudf*Lambdaf*(GuessTc-Ta);
        PconPrime=pi*(Nudf*Lambdaf+(GuessTc-Ta)*(LambdafPrime*Nudf+NudfPrimedtc*Lambdaf));
        PconPrimePrime=pi*(Nudf*LambdafPrime+NudfPrimedtc*Lambdaf+(LambdafPrime*Nudf+NudfPrimedtc*Lambdaf)+(GuessTc-Ta)*(LambdafPrime*NudfPrimedtc+NudfPrimePrimedtc2*Lambdaf+NudfPrimedtc*LambdafPrime));
    elseif(GuessTc~=Ta && Vw ~=0)
    %mixed convection
        Nudf=fGrPr(GrPr);
        NudfPrimedgrpr=differentiate(fGrPr,GrPr);
        NudfPrimedtc=NudfPrimedgrpr*(Gr*PrPrime+GrPrime*Pr);
        Req=fNuRe(Nudf);
        ReqPrimednudf=differentiate(fNuRe,Nudf);
        ReqPrimedtc=ReqPrimednudf*NudfPrimedtc;
        Re=sin(phi)*Vw*D/vf;
        RePrimedtc=-(sin(phi)*Vw*D*vfPrime)/(vf^2);
        Reeff=sqrt((Re^2)+(Req^2));
        Top=Re*RePrimedtc+Req*ReqPrimedtc;% (GrPrime*Pr+Gr*PrPrime);
        ReeffPrimedtc=Top/Reeff;
        ReeffPrimePrimedtc2=(Reeff*dtop-Top*ReeffPrimedtc)/(Reeff^2);
        Nueff=fReNu(Reeff);
        [NueffPrimedreeff,NueffPrimePrimedreeff2]=differentiate(fReNu,Reeff);
        NueffPrimedtc=NueffPrimedreeff*ReeffPrimedtc;
        NueffPrimePrimedtc2=NueffPrimePrimedreeff2*ReeffPrimedtc+NueffPrimedreeff*ReeffPrimePrimedtc2;
        Pcon=pi*Nueff*Lambdaf*(GuessTc-Ta);
        PconPrime=pi*(Nueff*Lambdaf+(GuessTc-Ta)*(LambdafPrime*Nueff+NueffPrimedtc*Lambdaf));
        PconPrimePrime=pi*(Nueff*LambdafPrime+NueffPrimedtc*Lambdaf+LambdafPrime*Nueff+NueffPrimedtc*Lambdaf+(GuessTc-Ta)*(LambdafPrime*NueffPrimedtc+NueffPrimePrimedtc2*Lambdaf+NueffPrimedtc*LambdafPrime));
    elseif(GuessTc==Ta && Vw~=0)
    %pure forced    
        Re=sin(phi)*Vw*D/vf;  
        RePrimedtc=-(sin(phi)*Vw*D*vfPrime)/(vf^2);
        RePrimePrimedtc2=(-1*(vf)*sin(phi)*Vw*D*vfPrimePrime+sin(phi)*Vw*D*vfPrime*2*vfPrime)/(vf^3);
        Nu=fReNu(Re);        
        [NuPrimedre,NuPrimePrimedre2]=differentiate(fReNu,Re);
        NuPrimedtc=NuPrimedre*RePrimedtc;
        NuPrimePrimedtc2=NuPrimePrimedre2*RePrimedtc+NuPrimedre*RePrimePrimedtc2;
        Pcon=pi*Nu*Lambdaf*(GuessTc-Ta);
        PconPrime=pi*(Nu*Lambdaf+(GuessTc-Ta)*(LambdafPrime*Nu+NuPrimedtc*Lambdaf));
        PconPrimePrime=pi*(Nu*LambdafPrime+NuPrimedtc*Lambdaf+LambdafPrime*Nu+NuPrimedtc*Lambdaf+(GuessTc-Ta)*(LambdafPrime*NuPrimedtc+NuPrimePrimedtc2*Lambdaf+NuPrimedtc*LambdafPrime));
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