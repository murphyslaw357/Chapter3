function [GuessTc,I2R,I2Rprime,Prad,PradPrime,Pcon,PconPrime] =GetTempNewton(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,Psol,fGrPr,fReNu,fNuRe,NudfPrimedgrpr,NudfPrimePrimedgrpr2,mdl)
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
    [GuessTc]=GetGuessTemp(I,Ta,D,phi,Vw,alpha,beta,epsilons,Psol*D,mdl);
    GuessTc2=GuessTc;
    Tolerance=0.001; %tolerance criteria for Newton's method
    update=realmax;
    TfilmkPrime=1/2;
    PrPrime=-1.25e-4;
    vfPrime=(4.75e-8)*((1-((6.5e-3)*H)/288.16)^-5.2561);

    IIstar=abs(I)^2;
    I2Rprime=beta*IIstar;

    while(abs(update)>Tolerance)
        counter=counter+1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RADIATIVE COOLING/HEATING%%%%%%%%%%%
        Prad=pi*D*sigmab*epsilons*(((GuessTc+273)^4)-((Ta+273)^4));
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
        Gr=(g*(D^3).*abs(GuessTc-Ta))./(Tfilmk.*(vf.^2));  
        GrPrime=g*(D^3).*((Tfilmk.*(vf.^2)).*((GuessTc-Ta)./abs(GuessTc-Ta))-abs(GuessTc-Ta).*(TfilmkPrime*vf.^2+...
        Tfilmk*2.*vf.*vfPrime))/((Tfilmk^2).*(vf^4));

        Lambdaf=(2.42e-2)+(7.2e-5)*GuessTfilm;
        LambdafPrime=3.6e-5;
        
        Pr=0.715-(2.5e-4)*GuessTfilm;
        GrPr=Gr*Pr;
        %Natural convection
        Nudf=fGrPr(GrPr);
        NudfPrimedtc=NudfPrimedgrpr(GrPr).*(Gr.*PrPrime+GrPrime.*Pr);
    
        %Mixed convection
        Re=(sin(phi)*Vw*D)./vf;  
        RePrimedtc=-((sin(phi)*Vw*D).*vfPrime)./(vf.^2);
        Req=fNuRe(Nudf);
        ReqPrimednudf=differentiate(fNuRe,Nudf);
        ReqPrimedtc=ReqPrimednudf.*NudfPrimedtc;
        Reeff=sqrt((Re.^2)+(Req.^2));
        Top=Re.*RePrimedtc+Req.*ReqPrimedtc;
        ReeffPrimedtc=Top./Reeff;
        Nueff=fReNu(Reeff);
        NueffPrimedreeff=differentiate(fReNu,Reeff);
        NueffPrimedtc=NueffPrimedreeff.*ReeffPrimedtc;
        
        %Forced convection
        Nu=fReNu(Re);  
        NuPrimedre=differentiate(fReNu,Re);
        NuPrimedtc=NuPrimedre.*RePrimedtc;
        [row,~]=size(GuessTc);
        Pcon=zeros(row,1);
        PconPrime=zeros(row,1);
        
        if(GuessTc~=Ta && Vw==0)
            %pure natural convection         
            Pcon=pi*Nudf*Lambdaf*(GuessTc-Ta);
            PconPrime=pi*(Nudf*Lambdaf+(GuessTc-Ta)*(LambdafPrime*Nudf+NudfPrimedtc*Lambdaf));
        elseif(GuessTc~=Ta && Vw ~=0)
            %mixed convection
            Pcon=pi*Nueff*Lambdaf*(GuessTc-Ta);
            PconPrime=pi*(Nueff*Lambdaf+(GuessTc-Ta)*(LambdafPrime*Nueff+NueffPrimedtc*Lambdaf));
        elseif(GuessTc==Ta && Vw~=0)
            %pure forced    
            Pcon=pi*Nu*Lambdaf*(GuessTc-Ta);
            PconPrime=pi*(Nu*Lambdaf+(GuessTc-Ta)*(LambdafPrime*Nu+NuPrimedtc*Lambdaf));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MISMATCH AND UPDATE%%%%%%%%%%%%%%
        Hprime=I2Rprime-PradPrime-PconPrime;
        if(Hprime>0 || isnan(PconPrime))
            msg='error condition';
            error(msg);
        end

        Mismatch=I2R+Psol*D-Prad-Pcon;
        update=Mismatch/Hprime;
        GuessTc=GuessTc-update;   

        if(counter>=5000)
            GuessTc=nan;
            break;
        end
    end
    if(abs(GuessTc2-GuessTc)>100)
    end
    
    if(isnan(GuessTc))
        msg='converge failed!';
        error(msg);
    elseif(round(I2R,1)<0 || round(GuessTc-Ta,1)<0)
       msg='error condition!';
       error(msg)
    end
end