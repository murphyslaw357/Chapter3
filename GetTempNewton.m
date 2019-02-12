function [GuessTc,I2R,I2Rprime,Prad,PradPrime,Pcon,PconPrime] =GetTempNewton(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,Psol)
    %I - RMS steady-state load current - amps
    %Ta - ambient temperature - degc
    %H - conductor elevation - meters
    %D - conductor diameter - meters
    %phi - angle between the wind direction and conductor axis - radians
    %Vw - Wind velocity - m/s
	%epsilons - conductor emissivity
    %Psol - solar heating - w/m  
    counter=0;
    sigmab=5.6697e-8;
    g=9.805;
    IIstar=abs(I)^2;
    GuessTc=(Psol+IIstar*(alpha+25*beta))/(pi*D*sigmab*epsilons*((1.38e8)+Ta*(1.39e6))+pi*(2.42e-2)*0.645)+Ta;
    %GuessTc=((Psol+IIstar*(alpha+25*beta))/(pi*D*sigmab*epsilons)+((Ta+273)^4))^(1/4)-273;  
    Tolerance=0.001; %tolerance criteria for Newton's method
    epsilon=1e-9;
    update=realmax;
    TfilmkPrime=1/2;
    PrPrime=-1.25e-4;
    vfPrime=(4.75e-8)*((1-((6.5e-3)*H)/288.16)^-5.2561);

%     Relim1=0.437*(4e-3)^0.0895;
%     Relim2=0.565*(9e-2)^0.136;
%     Relim3=0.800*(1)^0.280;
%     Relim4=0.795*(35)^0.384;
%     Relim5=0.583*(5e3)^0.471;
%     Relim6=0.148*(5e4)^0.633;
    Relim1=(0.437/0.565)^(1/(0.136-0.0895));
    Relim2=(0.565/0.8)^(1/(0.280-0.136));
    Relim3=(0.8/0.795)^(1/(0.384-0.280));
    Relim4=(0.795/0.583)^(1/(0.471-0.384));
    Relim5=(0.583/0.148)^(1/(0.633-0.471));
    Relim6=(0.148/0.0208)^(1/(0.814-0.633));
    Nulim1=0.437*(Relim1^0.0895);
    Nulim2=0.565*(Relim2^0.136);
    Nulim3=0.800*(Relim3^0.280);
    Nulim4=0.795*(Relim4^0.384);
    Nulim5=0.583*(Relim5^0.471);
    Nulim6=0.148*(Relim6^0.633);
    GrPrlim1=(0.675/0.889)^(1/(0.088-0.058));
    GrPrlim2=(0.889/1.02)^(1/(0.148-0.088));
    GrPrlim3=(1.02/0.850)^(1/(0.188-0.148));
    GrPrlim4=(0.850/0.480)^(1/(0.250-0.188));
    GrPrlim5=(0.480/0.125)^(1/(0.333-0.250));
    
%     if(GuessTc<Ta)
%         GuessTc=Ta+Tolerance;
%     end
%     I2R=0;
%     Prad=0;
%     Pcon=0;
    IIstar=abs(I)^2;
    I2Rprime=beta*IIstar;
    
    while(abs(update)>Tolerance)
        counter=counter+1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RADIATIVE COOLING/HEATING%%%%%%%%%%%
        %Tsky=(0.0552*(Ta+273)^1.5)-273;
        %Tg=Ta+2;
        Prad=pi*D*sigmab*epsilons*(((GuessTc+273)^4)-((Ta+273)^4));
        %Prad=pi*D*sigmab*epsilons*(((GuessTc+273)^4)-0.5*((Tsky+273)^4)-0.5*((Tg+273)^4));
        PradPrime=4*pi*D*sigmab*epsilons*(GuessTc+273)^3;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%JOULE HEATING%%%%%%%%%%%%%%%%%%%%%%%
        Resistance=beta*GuessTc+alpha;
        if(Resistance<0)
            Resistance=0;
        end
        I2R=IIstar*Resistance;
        %%%%%%%%%%%%%%%%%%%%%%%%CONVECTIVE COOLING/HEATING%%%%%%%%%%%%%%%%%
        GuessTfilm=(GuessTc+Ta)/2;
        vf=((1.32e-5)+(9.5e-8)*GuessTfilm)*((1-((6.5e-3)*H)/288.16)^-5.2561);

        Tfilmk=GuessTfilm+273;

        Gr=(g*(D^3)*abs(GuessTc-Ta))/(Tfilmk*(vf^2));    
        Lambdaf=(2.42e-2)+(7.2e-5)*GuessTfilm;
        LambdafPrime=3.6e-5;
        GrPrime=g*(D^3)*(Tfilmk*(vf)-abs(GuessTc-Ta)*(TfilmkPrime*vf+Tfilmk*2*vfPrime))/...
                ((Tfilmk^2)*(vf^3));
        Pr=0.715-(2.5e-4)*GuessTfilm;
        GrPr=Gr*Pr;
        if(GuessTc~=Ta && Vw==0)
            %pure natural convection
            %lookup A and m
            A=0;
            m=0;
            if(GrPr>=0 && GrPr<=GrPrlim1) %1e-10
                A=0.675;
                m=0.058;
            elseif(GrPr>GrPrlim1 && GrPr<=GrPrlim2)
                A=0.889;
                m=0.088;
            elseif(GrPr>GrPrlim2 && GrPr<=GrPrlim3)
                A=1.02;
                m=0.148;
            elseif(GrPr>GrPrlim3 && GrPr<=GrPrlim4)
                A=0.85;
                m=0.188;
            elseif(GrPr>GrPrlim4 && GrPr<=GrPrlim5)
                A=0.48;
                m=0.25;
            elseif(GrPr>GrPrlim5 && GrPr<=1e12)
                A=0.125;
                m=0.333;
            end
            if(A==0)
                msg='error';
                error(msg)
            end
            Nudf=A*(GrPr^m);
            NudfPrime=m*A*(GrPr^(m-1))*(Gr*PrPrime+Pr*GrPrime);    
            Pcon=pi*Nudf*Lambdaf*(GuessTc-Ta);
            PconPrime=pi*(Nudf*Lambdaf+LambdafPrime*Nudf*(GuessTc-Ta)+NudfPrime*Lambdaf*(GuessTc-Ta));
            
%             f2=Gr*Pr;
%             f2prime=Gr*PrPrime+GrPrime*Pr;
%             f3=(1+((0.559/Pr)^(9/16)))^(16/9);
%             f3prime=((1+(0.559/Pr)^(9/16))^(7/9))*((0.559/Pr)^(-7/16))*(-1)*(0.559/(Pr^2));
%             f1=f2/f3;
%             f1prime=(f3*f2prime-f2*f3prime)/(f3^2);
%             Nu2=(0.6+0.387*f1^(1/6))^2;
%             Nu2Prime=2*(0.6+0.387*f1^(1/6))*(1/6)*(0.387*f1^(-5/6))*f1prime;
%             Pcon=pi*Nu2*Lambdaf*(GuessTc-Ta);
%             PconPrime=pi*(Nu2*Lambdaf+LambdafPrime*Nu2*(GuessTc-Ta)+Nu2Prime*Lambdaf*(GuessTc-Ta));
        elseif(GuessTc~=Ta && Vw ~=0)
        %mixed convection
            GrPrime=g*(D^3)*(Tfilmk*(vf)-abs(GuessTc-Ta)*(TfilmkPrime*vf+Tfilmk*2*vfPrime))/...
                ((Tfilmk^2)*(vf^3));
            Pr=0.715-(2.5e-4)*GuessTfilm;
            GrPr=Gr*Pr;
            %lookup A and m
            A=0;
            m=0;
            if(GrPr>=0 && GrPr<=GrPrlim1) %1e-10
                A=0.675;
                m=0.058;
            elseif(GrPr>GrPrlim1 && GrPr<=GrPrlim2)
                A=0.889;
                m=0.088;
            elseif(GrPr>GrPrlim2 && GrPr<=GrPrlim3)
                A=1.02;
                m=0.148;
            elseif(GrPr>GrPrlim3 && GrPr<=GrPrlim4)
                A=0.85;
                m=0.188;
            elseif(GrPr>GrPrlim4 && GrPr<=GrPrlim5)
                A=0.48;
                m=0.25;
            elseif(GrPr>GrPrlim5 && GrPr<=1e12)
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
            if(Nudf<=Nulim1)
                C=0.437;
                n=0.0895;
            elseif(Nudf>Nulim1 && Nudf<=Nulim2)
                C=0.565;
                n=0.136;
            elseif(Nudf>Nulim2 && Nudf<=Nulim3)
                C=0.8;
                n=0.28;
            elseif(Nudf>Nulim3 && Nudf<=Nulim4)
                C=0.795;
                n=0.384;
            elseif(Nudf>Nulim4 && Nudf<=Nulim5)
                C=0.583;
                n=0.471;
            elseif(Nudf>Nulim5 && Nudf<=Nulim6)
                C=0.148;
                n=0.633;
            elseif(Nudf>Nulim6)
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
            if(Reeff<=Relim1)
                C=0.437;
                n=0.0895;
            elseif(Reeff>Relim1 && Reeff<=Relim2)
                C=0.565;
                n=0.136;
            elseif(Reeff>Relim2 && Reeff<=Relim3)
                C=0.8;
                n=0.28;
            elseif(Reeff>Relim3 && Reeff<=Relim4)
                C=0.795;
                n=0.384;
            elseif(Reeff>Relim4 && Reeff<=Relim5)
                C=0.583;
                n=0.471;
            elseif(Reeff>Relim5 && Reeff<=Relim6)
                C=0.148;
                n=0.633;
            elseif(Reeff>Relim6)
                C=0.0208;
                n=0.814;
            end
            if(C==0)
                msg='error';
                error(msg)
            end

            Top=Re*RePrime+((A/C)^(2/n))*(m/n)*(GrPr^((2*m-n)/n))*(GrPrime*Pr+Gr*PrPrime);

            ReeffPrime=Top/Reeff;
            Nueff=C*Reeff^n;
            NueffPrime=C*n*(Reeff^(n-1))*ReeffPrime;        

            Pcon=pi*Nueff*Lambdaf*(GuessTc-Ta);
            PconPrime=pi*(Nueff*Lambdaf+LambdafPrime*Nueff*(GuessTc-Ta)+NueffPrime*Lambdaf*(GuessTc-Ta));
        elseif(GuessTc==Ta && Vw~=0)
        %pure forced    
            Re=sin(phi)*Vw*D/vf;
            RePrime=-(sin(phi)*Vw*D*vfPrime)/(vf^2);
            %%
            %lookup C and n for Reynolds number
            C=0;
            n=0;
            if(Re<=Relim1)
                C=0.437;
                n=0.0895;
            elseif(Re>Relim1&& Re<=Relim2)
                C=0.565;
                n=0.136;
            elseif(Re>Relim2 && Re<=Relim3)
                C=0.8;
                n=0.28;
            elseif(Re>Relim3 && Re<=Relim4)
                C=0.795;
                n=0.384;
            elseif(Re>Relim4 && Re<=Relim5)
                C=0.583;
                n=0.471;
            elseif(Re>Relim5 && Re<=Relim6)
                C=0.148;
                n=0.633;
            elseif(Re>Relim6)
                C=0.0208;
                n=0.814;
            end
            if(C==0)
                msg='error';
                error(msg)
            end

            Nu=C*Re^n;        
            Pcon=pi*Nu*Lambdaf*(GuessTc-Ta);
            NuPrime=C*n*(Re^(n-1))*RePrime;   
            PconPrime=pi*(Nu*Lambdaf+LambdafPrime*Nu*(GuessTc-Ta)+NuPrime*Lambdaf*(GuessTc-Ta));
        else
        %no convection at all    
            Pcon=0;
            PconPrime=0;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MISMATCH AND UPDATE%%%%%%%%%%%%%%
        Hprime=I2Rprime-PradPrime-PconPrime;
        if(Hprime>0 || isnan(PconPrime) || PconPrime<0 || PradPrime<0)
            msg='error condition';
            error(msg);
        end

        Mismatch=I2R+Psol-Prad-Pcon;
        update=Mismatch/Hprime;
        GuessTc=GuessTc-update;   

        %if(counter>=4000)
        %end
        if(counter>=5000)
            GuessTc=nan;
            break;
        end
    end
    if(isnan(GuessTc))
        msg='converge failed!';
        error(msg);
    elseif(round(I2R,1)<0 || Psol <0 || round(GuessTc-Ta,1)<0)
       msg='error condition!';
       error(msg)
    end
end