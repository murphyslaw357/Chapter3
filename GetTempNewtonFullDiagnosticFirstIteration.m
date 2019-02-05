function [GuessTc,I2R,I2Rprime,Prad,PradPrime,PradPrimePrime,Pcon,PconPrime,PconPrimePrime] =GetTempNewtonFullDiagnosticFirstIteration(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,Psol,GuessTc)
    phi=90*pi/180;
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
    TfilmkPrime=1/2;
    PrPrime=-1.25e-4;
    vfPrime=(4.75e-8)*((1-((6.5e-3)*H)/288.16)^-5.2561);
    vfPrimePrime=0;
    Relim1=0.437*(4e-3)^0.0895;
    Relim2=0.565*(9e-2)^0.136;
    Relim3=0.800*(1)^0.280;
    Relim4=0.795*(35)^0.384;
    Relim5=0.583*(5e3)^0.471;
    Relim6=0.148*(5e4)^0.633;

    I2Rprime=beta*IIstar;
    
    GuessTfilm=(GuessTc+Ta)/2;
    vf=((1.32e-5)+(9.5e-8)*GuessTfilm)*((1-((6.5e-3)*H)/288.16)^-5.2561);

    Tfilmk=GuessTfilm+273;

    Gr=(g*(D^3)*abs(GuessTc-Ta))/(Tfilmk*(vf^2));        
    GrPrime=g*(D^3)*(Tfilmk*(vf)-abs(GuessTc-Ta)*(TfilmkPrime*vf+Tfilmk*2*vfPrime))/...
        ((Tfilmk^2)*(vf^3));

    Pr=0.715-(2.5e-4)*GuessTfilm;
    GrPr=Gr*Pr;
    %lookup A and m
    A=0;
    m=0;
    if(GrPr>=0 && GrPr<=1e-4) %1e-10
        A=0.675;
        m=0.058;
    elseif(Gr*Pr>1e-4 && Gr*Pr<=1e-1)
        A=0.889;
        m=0.088;
    elseif(Gr*Pr>1e-1 && Gr*Pr<=1e2)
        A=1.02;
        m=0.148;
    elseif(Gr*Pr>1e2 && Gr*Pr<=1e4)
        A=0.85;
        m=0.188;
    elseif(Gr*Pr>1e4 && Gr*Pr<=1e7)
        A=0.48;
        m=0.25;
    elseif(Gr*Pr>1e7 && Gr*Pr<=1e12)
        A=0.125;
        m=0.333;
    end
    if(A==0)
        msg='error';
        error(msg)
    end
    Nudf=A*(GrPr^m);
    %lookup C and n
    C=0;
    n=0;
    if(Nudf<=Relim1)
        C=0.437;
        n=0.0895;
    elseif(Nudf>Relim1 && Nudf<=Relim2)
        C=0.565;
        n=0.136;
    elseif(Nudf>Relim2 && Nudf<=Relim3)
        C=0.8;
        n=0.28;
    elseif(Nudf>Relim3 && Nudf<=Relim4)
        C=0.795;
        n=0.384;
    elseif(Nudf>Relim4 && Nudf<=Relim5)
        C=0.583;
        n=0.471;
    elseif(Nudf>Relim5 && Nudf<=Relim6)
        C=0.148;
        n=0.633;
    elseif(Nudf>Relim6)
        C=0.0208;
        n=0.814;
    end
    if(C==0)
        msg='error';
        error(msg)
    end
    Req=(Nudf/C)^(1/n);

    %k=((A/C)*Pr^m)^(1/n);
    %Req=k*(Gr^(m/n));
    %k1=((A/C)*Pr^m)^(1/n)
    %Req1=((A/C)*GrPr^m)^(1/n);%*(Gr^(m/n))

    Re=sin(phi)*Vw*D/vf;

    RePrime=-(sin(phi)*Vw*D*vfPrime)/(vf^2);
    RePrimePrime=sin(phi)*Vw*D*2*(vfPrime^2)/(vf^3);
    Reeff=sqrt((Re^2)+(Req^2));

    %%
    %re-lookup C and n for effective Reynolds number
    C=0;
    n=0;
    if(Reeff<=4e-3)
        C=0.437;
        n=0.0895;
    elseif(Reeff>4e-3 && Reeff<=9e-2)
        C=0.565;
        n=0.136;
    elseif(Reeff>9e-2 && Reeff<=1)
        C=0.8;
        n=0.28;
    elseif(Reeff>1 && Reeff<=35)
        C=0.795;
        n=0.384;
    elseif(Reeff>35 && Reeff<=5e3)
        C=0.583;
        n=0.471;
    elseif(Reeff>5e3 && Reeff<=5e4)
        C=0.148;
        n=0.633;
    elseif(Reeff>5e4)% && Re<2e5)
        C=0.0208;
        n=0.814;
    end
    if(C==0)
        msg='error';
        error(msg)
    end
    Resistance=beta*GuessTc+alpha;
    if(Resistance<0)
        Resistance=0;
    end
    I2R=IIstar*Resistance;

    %Tsky=(0.0552*(Ta+273)^1.5)-273;
    %Tg=Ta+2;
    Prad=pi*D*sigmab*epsilons*(((GuessTc+273)^4)-((Ta+273)^4));
    %Prad=pi*D*sigmab*epsilons*(((GuessTc+273)^4)-0.5*((Tsky+273)^4)-0.5*((Tg+273)^4));
    PradPrime=4*pi*D*sigmab*epsilons*(GuessTc+273)^3;
    PradPrimePrime=12*pi*D*sigmab*epsilons*(GuessTc+273)^2;

    dtop=TfilmkPrime*vf+Tfilmk*vfPrime-abs(GuessTc-Ta)*(TfilmkPrime*vfPrime+2*TfilmkPrime*vfPrime+2*Tfilmk*vfPrimePrime)-(TfilmkPrime*vf+2*Tfilmk*vfPrime);
    GrPrimePrime=g*(D^3)*(((Tfilmk^2)*(vf^3))*dtop-...
        ((Tfilmk*vf-abs(GuessTc-Ta)*(TfilmkPrime*vf+2*Tfilmk*vfPrime)))*(Tfilmk*(vf^3)+3*(vf^2)*vfPrime*(Tfilmk^2)))/...
        ((Tfilmk^4)*(vf^6));


    Top=Re*RePrime+((A/C)^(2/n))*(m/n)*(GrPr^((2*m-n)/n))*(GrPrime*Pr+Gr*PrPrime);
    TopPrime=Re*RePrimePrime+RePrime*RePrime+((A/C)^(2/n))*(m/n)*(((2*m-n)/n)*(GrPr^((2*m-2*n)/n))*(GrPrime*Pr+Gr*PrPrime)+...
        (GrPr^((2*m-n)/n))*(GrPrimePrime-2*PrPrime*GrPrime));
    ReeffPrime=Top/Reeff;
    ReeffPrimePrime=(Reeff*TopPrime-Top*ReeffPrime)/(Reeff^2);
    Nueff=C*Reeff^n;
    NueffPrime=C*n*(Reeff^(n-1))*ReeffPrime;        

    NueffPrimePrime=C*n*((n-1)*(Reeff^(n-2))*(ReeffPrime^2)+(Reeff^(n-1))*ReeffPrimePrime);

    Lambdaf=(2.42e-2)+(7.2e-5)*GuessTfilm;
    LambdafPrime=3.6e-5;
    Pcon=pi*Nueff*Lambdaf*(GuessTc-Ta);
    PconPrime=pi*(Nueff*Lambdaf+LambdafPrime*Nueff*(GuessTc-Ta)+NueffPrime*Lambdaf*(GuessTc-Ta));
    PconPrimePrime=pi*(NueffPrimePrime*Lambdaf*(GuessTc-Ta)+2*NueffPrime*(LambdafPrime*(GuessTc-Ta)+Lambdaf)+2*Nueff*LambdafPrime);
    Hprime=I2Rprime-PradPrime-PconPrime;

    if(Hprime>0)
        disp('doh')
    end
    Mismatch=I2R+Psol-Prad-Pcon;
    update=Mismatch/Hprime;
    GuessTc=GuessTc-update;      
    
end