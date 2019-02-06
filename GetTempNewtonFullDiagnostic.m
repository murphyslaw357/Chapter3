function [GuessTc,I2R,I2Rprime,Prad,PradPrime,Pcon,PconPrime] =GetTempNewtonFullDiagnostic(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,Psol)
    %phi=90*pi/180;
    %I - RMS steady-state load current - amps
    %Ta - ambient temperature - degc
    %H - conductor elevation - meters
    %D - conductor diameter - meters
    %phi - angle between the wind direction and conductor axis - radians
    %Vw - Wind velocity - m/s
	%epsilons - conductor emissivity
    %Psol - solar heating - w/m  
    %beta=(R_T_high-R_T_low)/(T_high-T_low);
    %alpha=R_T_high-beta*T_high;
    counter=0;
    sigmab=5.6697e-8;
    g=9.805;
    IIstar=abs(I)^2;
    GuessTc=((Psol+IIstar*(alpha+25*beta))/(pi*D*sigmab*epsilons)+((Ta+273)^4))^(1/4)-273;  
    Tolerance=0.001; %tolerance criteria for Newton's method
    epsilon=1e-9;
    update=realmax;
    TfilmkPrime=1/2;
    PrPrime=-1.25e-4;
    vfPrime=(4.75e-8)*((1-((6.5e-3)*H)/288.16)^-5.2561);
    Relim1=0.437*(4e-3)^0.0895;
    Relim2=0.565*(9e-2)^0.136;
    Relim3=0.800*(1)^0.280;
    Relim4=0.795*(35)^0.384;
    Relim5=0.583*(5e3)^0.471;
    Relim6=0.148*(5e4)^0.633;
    if(GuessTc<Ta)
        GuessTc=Ta+Tolerance;
    end
    I2R=0;
    Prad=0;
    Pcon=0;
    I2Rprime=beta*IIstar;
    if(isnan(GuessTc))
    
    end
    while(abs(update)>Tolerance)
        counter=counter+1;
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
        elseif(GrPr>1e-4 && GrPr<=1e-1)
            A=0.889;
            m=0.088;
        elseif(GrPr>1e-1 && GrPr<=1e2)
            A=1.02;
            m=0.148;
        elseif(GrPr>1e2 && GrPr<=1e4)
            A=0.85;
            m=0.188;
        elseif(GrPr>1e4 && GrPr<=1e7)
            A=0.48;
            m=0.25;
        elseif(GrPr>1e7 && GrPr<=1e12)
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
        Re=sin(phi)*Vw*D/vf;
        
        RePrime=-(sin(phi)*Vw*D*vfPrime)/(vf^2);
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
        
        Prad=pi*D*sigmab*epsilons*(((GuessTc+273)^4)-((Ta+273)^4));
        PradPrime=4*pi*D*sigmab*epsilons*(GuessTc+273)^3;
       
        
        Top=Re*RePrime+((A/C)^(2/n))*(m/n)*(GrPr^((2*m-n)/n))*(GrPrime*Pr+Gr*PrPrime);
        %TopPrime=Re*RePrimePrime+RePrime*RePrime+((A/C)^(2/n))*(m/n)*(((2*m-n)/n)*(GrPr^((2*m-2*n)/n))*(GrPrime*Pr+Gr*PrPrime)+...
        %    (GrPr^((2*m-n)/n))*(GrPrimePrime-2*PrPrime*GrPrime));
        ReeffPrime=Top/Reeff;
        %ReeffPrimePrime=(Reeff*TopPrime-Top*ReeffPrime)/(Reeff^2);
        Nueff=C*Reeff^n;
        NueffPrime=C*n*(Reeff^(n-1))*ReeffPrime;        
        
        %NueffPrimePrime=C*n*((n-1)*(Reeff^(n-2))*(ReeffPrime^2)+(Reeff^(n-1))*ReeffPrimePrime);
 
        Lambdaf=(2.42e-2)+(7.2e-5)*GuessTfilm;
        LambdafPrime=3.6e-5;
        Pcon=pi*Nueff*Lambdaf*(GuessTc-Ta);
        PconPrime=pi*(Nueff*Lambdaf+LambdafPrime*Nueff*(GuessTc-Ta)+NueffPrime*Lambdaf*(GuessTc-Ta));
        %PconPrimePrime=pi*(NueffPrimePrime*Lambdaf*(GuessTc-Ta)+2*NueffPrime*(LambdafPrime*(GuessTc-Ta)+Lambdaf)+2*Nueff*LambdafPrime);
        Hprime=I2Rprime-PradPrime-PconPrime;
        if(abs(Hprime)<epsilon) 
            GuessTc=nan; %denominator is too small- exit loop
            break;
        end
        if(Hprime>0)
            disp('doh')
        end
        Mismatch=I2R+Psol-Prad-Pcon;
        update=Mismatch/Hprime;
        if(isnan(update))
           disp('error') 
        end
        GuessTc=GuessTc-update;

        if(GuessTc<Ta)
            %disp('case')
        elseif(GuessTc>1500)
            GuessTc=1500;
        end
        
        if(counter>=5000)
            GuessTc=nan;
            break;
        end
    end
    if(isnan(GuessTc))
        msg='converge failed!';
        disp(msg);
    elseif(round(I2R,1)<0 || Psol <0 || round(GuessTc-Ta,1)<0)
       msg='error condition!';
       error(msg)
    end
    %disp(counter)
    %GuessTc=GuessTc*1.5;
    %GuessTc=GuessTc+10;
end