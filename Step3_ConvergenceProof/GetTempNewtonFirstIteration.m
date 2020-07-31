function [GuessTcOutput,I2R,I2Rprime,Prad,PradPrime,PradPrimePrime,Pcon,...
    PconPrime,PconPrimePrime,Gr,GrPrime,Nun,A,m,Cinv,ninv,C,n] = GetTempNewtonFirstIteration(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,alphas,...
    Psol,GuessTc,AmCinvninvCn)
    %% API
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
    gplim1=1e-10;
    gplim2=1e-4;
    gplim3=1e-1;
    gplim4=1e2;
    gplim5=1e4;
    gplim6=1e7;
    gplim7=1e12;
    
    IIstar=abs(I)^2;
    %%                          RADIATIVE COOLING/HEATING                %%
    Prad=(pi*D*sigmab*epsilons)*(((GuessTc+273)^4)-((Ta+273)^4));
    PradPrime=(4*pi*D*sigmab*epsilons)*(GuessTc+273)^3;
    PradPrimePrime=(12*pi*D*sigmab*epsilons)*(GuessTc+273)^2;
    %%                          JOULE HEATING                            %%
    I2Rprime=beta*IIstar;
    Resistance=beta*GuessTc+alpha;
    if(Resistance<0)
        Resistance=0;
    end
    I2R=IIstar*Resistance;
    %%                      CONVECTIVE COOLING/HEATING                   %%
    vfPrimePrime=0;
    PrPrime=-1.25e-4;
    GuessTfilm=(GuessTc+Ta)/2;
    vf=((1.32e-5)+(9.5e-8)*GuessTfilm)*((1-((6.5e-3)*H)/288.16)^-5.2561);
    vfPrime=(4.75e-8)*((1-((6.5e-3)*H)/288.16)^-5.2561);
    Tfilmk=GuessTfilm+273;
    Lambdaf=(2.42e-2)+(7.2e-5)*GuessTfilm;
    LambdafPrime=3.6e-5;
    
    gD3 = g*(D^3);
    Gr = gD3*(GuessTc-Ta)/(Tfilmk.*(vf^2));
    GrPrime = gD3*(Ta*vf-((GuessTc^2)-Ta^2)*vfPrime)/((Tfilmk^2)*(vf^3));
    GrPrimePrime=(-2*gD3*GuessTc/((Tfilmk.^2).*(vf.^3)))*vfPrime-...
        (gD3/((Tfilmk^3)*(vf^4)))*(vf+3*Tfilmk*vfPrime)*...
        (Ta*vf-((GuessTc^2)-Ta^2)*vfPrime);
    Pr=0.715-(2.5e-4)*GuessTfilm;
    GrPr=Gr*Pr;
    
    if(isempty(AmCinvninvCn))
        %Natural convection
        if(GrPr>gplim1 && GrPr<=gplim2)
            A=0.675;
            m=0.058;
        elseif(GrPr>gplim2 && GrPr<=gplim3)
            A=0.889;
            m=0.088;
        elseif(GrPr>gplim3 && GrPr<=gplim4)
            A=1.02;
            m=0.148;
        elseif(GrPr>gplim4 && GrPr<=gplim5)
            A=0.85;
            m=0.188;
        elseif(GrPr>gplim5 && GrPr<=gplim6)
            A=0.48;
            m=0.25;
        elseif(GrPr>gplim6 && GrPr<=gplim7)
            A=0.125;
            m=0.333;
        end
    else
        A=AmCinvninvCn(1);
        m=AmCinvninvCn(2);
    end
    Nun=A*GrPr^m;
    NunPrimedgrpr=m*A*GrPr^(m-1);
    NunPrimePrimedgrpr2=(m-1)*m*A*GrPr^(m-2);
    NunPrimedtc=NunPrimedgrpr*(Gr*PrPrime+GrPrime*Pr);
    NunPrimePrimedtc2=NunPrimePrimedgrpr2*(Gr*PrPrime+GrPrime*Pr)+...
        NunPrimedgrpr*(GrPrime*PrPrime+GrPrimePrime*Pr+GrPrime*PrPrime);

    %Mixed convection
    Re=(sin(phi)*Vw*D)/vf;  
    RePrimedtc=-((sin(phi)*Vw*D)*vfPrime)/(vf^2);
    RePrimePrimedtc2=(sin(phi)*Vw*D)*...
        (2.*(vfPrime.^2)-1.*vf.*vfPrimePrime)./(vf.^3);
    if(isempty(AmCinvninvCn))
        if(Nun>0.1916 && Nun<=0.2666)
            Cinv=0.437;
            ninv=0.0895;
        elseif(Nun>0.2666 && Nun<=0.40685)
           Cinv=0.565;
           ninv=0.136;
        elseif(Nun>0.40685 && Nun<=0.8136)
           Cinv=0.8;
           ninv=0.28;
        elseif(Nun>0.8136 && Nun<=3.1253)
           Cinv=0.795;
           ninv=0.384;
        elseif(Nun>3.1253 && Nun<=31.45)
           Cinv=0.583;
           ninv=0.471;
        elseif(Nun>31.45 && Nun<=141.449)
           Cinv=0.148;
           ninv=0.633;
        elseif(Nun>141.449 && Nun<=429.6371)
           Cinv=0.0208;
           ninv=0.814;
        else
            error('out of bounds')
        end  
    else
        Cinv=AmCinvninvCn(3);
        ninv=AmCinvninvCn(4);
    end
    Ren=(Nun/Cinv)^(1/ninv);
    RenPrimednudf=(1/ninv)*(Nun/Cinv)^(1/ninv-1);    
    RenPrimePrimednudf=(1/ninv-1)*(1/ninv)*(Nun/Cinv)^(1/ninv-2); 
    RenPrimedtc=RenPrimednudf.*NunPrimedtc;
    RenPrimePrimedtc2=RenPrimePrimednudf.*NunPrimedtc+RenPrimednudf.*...
        NunPrimePrimedtc2;
    Reeff=sqrt((Re.^2)+(Ren.^2));
    Top=Re.*RePrimedtc+Ren.*RenPrimedtc;
    dtop=RePrimedtc.*RePrimedtc+Re.*RePrimePrimedtc2+RenPrimedtc.*...
        RenPrimedtc+Ren.*RenPrimePrimedtc2;
    
    ReeffPrimedtc=Top./Reeff;
    ReeffPrimePrimedtc2=(Reeff.*dtop-Top.*ReeffPrimedtc)./(Reeff.^2);
    if(isempty(AmCinvninvCn))
        if(Reeff>1e-4 && Reeff<=4e-3)
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
        elseif(Reeff>5e4 && Reeff<2e5)
            C=0.0208;
            n=0.814;
        end
    else
        C=AmCinvninvCn(5);
        n=AmCinvninvCn(6);
    end
    Nueff=C*(Reeff^n);
    NueffPrimedreeff=C*n*(Reeff^(n-1));
    NueffPrimePrimedreeff2=(n-1)*C*n*(Reeff^(n-2));
    NueffPrimedtc=NueffPrimedreeff.*ReeffPrimedtc;
    NueffPrimePrimedtc2=NueffPrimePrimedreeff2.*ReeffPrimedtc+...
        NueffPrimedreeff.*ReeffPrimePrimedtc2;
    
    GuessTcOutput=GuessTc;
    %mixed convection
    Pcon=pi*Nueff*Lambdaf*(GuessTc-Ta);
    PconPrime=pi*(Nueff*Lambdaf+(GuessTc-Ta)*(LambdafPrime*Nueff+NueffPrimedtc*Lambdaf));
    PconPrimePrime=pi*(Nueff*LambdafPrime+NueffPrimedtc*Lambdaf+...
        LambdafPrime*Nueff+NueffPrimedtc*Lambdaf+...
        (GuessTc-Ta)*(LambdafPrime*NueffPrimedtc+...
        NueffPrimePrimedtc2*Lambdaf+NueffPrimedtc*LambdafPrime));

    %%                        MISMATCH AND UPDATE                        %%
    Hprime=I2Rprime-PradPrime-PconPrime;
    if(Hprime>0 || isnan(PconPrime))
        msg='Hprime greater than zero or PconPrime is nan';
        error(msg);
    end
    
    Mismatch=I2R+Psol*D*alphas-Prad-Pcon;
    update=Mismatch/Hprime;
    GuessTcOutput=GuessTcOutput-update;  
    if(~isreal(GuessTcOutput))
        error('Temperature process non-real')
    end
end