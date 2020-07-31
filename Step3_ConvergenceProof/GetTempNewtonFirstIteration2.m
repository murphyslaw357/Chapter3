function [GuessTcOutput,I2R,I2Rprime,Prad,PradPrime,PradPrimePrime,Pcon,...
    PconPrime,PconPrimePrime,Gr,GrPrime,Nun,A,m,Cinv,ninv,C,n] = GetTempNewtonFirstIteration2(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,alphas,...
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
    IIstar=abs(I).^2;
    %%                          RADIATIVE COOLING/HEATING                %%
    Prad=(pi*D*sigmab*epsilons).*(((GuessTc+273).^4)-((Ta+273).^4));
    PradPrime=(4*pi*D*sigmab*epsilons).*(GuessTc+273).^3;
    PradPrimePrime=(12*pi*D*sigmab*epsilons).*(GuessTc+273).^2;
    %%                          JOULE HEATING                            %%
    I2Rprime=beta.*IIstar;
    Resistance=beta.*GuessTc+alpha;
    Resistance(Resistance<0)=0;
    I2R=IIstar.*Resistance;
    %%                      CONVECTIVE COOLING/HEATING                   %%
    vfPrimePrime=0;
    PrPrime=-1.25e-4;
    GuessTfilm=(GuessTc+Ta)./2;
    vf=((1.32e-5)+(9.5e-8).*GuessTfilm).*((1-((6.5e-3)*H)/288.16)^-5.2561);
    vfPrime=(4.75e-8)*((1-((6.5e-3)*H)/288.16)^-5.2561);
    Tfilmk=GuessTfilm+273;
    Lambdaf=(2.42e-2)+(7.2e-5).*GuessTfilm;
    LambdafPrime=3.6e-5;
    
    gD3 = g*(D^3);
    Gr = gD3.*(GuessTc-Ta)./(Tfilmk.*(vf.^2));
    GrPrime = gD3.*(Ta.*vf-((GuessTc.^2)-(Ta.^2)).*vfPrime)./((Tfilmk.^2).*(vf.^3));
    GrPrimePrime=((-2*gD3).*GuessTc./((Tfilmk.^2).*(vf.^3))).*vfPrime-...
        (gD3./((Tfilmk.^3).*(vf.^4))).*(vf+3.*Tfilmk.*vfPrime).*...
        (Ta.*vf-((GuessTc.^2)-Ta.^2).*vfPrime);
    Pr=0.715-(2.5e-4).*GuessTfilm;
    GrPr=Gr.*Pr;
    
    A=AmCinvninvCn(1);
    m=AmCinvninvCn(2);
    Nun=A.*GrPr.^m;
    NunPrimedgrpr=m.*A.*GrPr.^(m-1);
    NunPrimePrimedgrpr2=((m-1)*m*A).*(GrPr.^(m-2));
    NunPrimedtc=NunPrimedgrpr.*(Gr.*PrPrime+GrPrime.*Pr);
    NunPrimePrimedtc2=NunPrimePrimedgrpr2.*(Gr.*PrPrime+GrPrime.*Pr)+...
        NunPrimedgrpr.*(GrPrime.*PrPrime+GrPrimePrime.*Pr+GrPrime.*PrPrime);

    %Mixed convection
    Re=(sin(phi).*Vw.*D)./vf;  
    RePrimedtc=-((sin(phi).*Vw.*D).*vfPrime)./(vf.^2);
    RePrimePrimedtc2=(sin(phi).*Vw.*D).*...
        (2.*(vfPrime.^2)-1.*vf.*vfPrimePrime)./(vf.^3);
    Cinv=AmCinvninvCn(3);
    ninv=AmCinvninvCn(4);
    
    Ren=(Nun./Cinv).^(1/ninv);
    RenPrimednudf=(1/ninv).*(Nun./Cinv).^(1/ninv-1);    
    RenPrimePrimednudf=(1/ninv-1).*(1/ninv).*(Nun./Cinv).^(1/ninv-2); 
    RenPrimedtc=RenPrimednudf.*NunPrimedtc;
    RenPrimePrimedtc2=RenPrimePrimednudf.*NunPrimedtc+RenPrimednudf.*...
        NunPrimePrimedtc2;
    Reeff=sqrt((Re.^2)+(Ren.^2));
    Top=Re.*RePrimedtc+Ren.*RenPrimedtc;
    dtop=RePrimedtc.*RePrimedtc+Re.*RePrimePrimedtc2+RenPrimedtc.*...
        RenPrimedtc+Ren.*RenPrimePrimedtc2;
    
    ReeffPrimedtc=Top./Reeff;
    ReeffPrimePrimedtc2=(Reeff.*dtop-Top.*ReeffPrimedtc)./(Reeff.^2);
    
    C=AmCinvninvCn(5);
    n=AmCinvninvCn(6);
    
    Nueff=C.*(Reeff.^n);
    NueffPrimedreeff=C.*n.*(Reeff.^(n-1));
    NueffPrimePrimedreeff2=(n-1).*C.*n.*(Reeff.^(n-2));
    NueffPrimedtc=NueffPrimedreeff.*ReeffPrimedtc;
    NueffPrimePrimedtc2=NueffPrimePrimedreeff2.*ReeffPrimedtc+...
        NueffPrimedreeff.*ReeffPrimePrimedtc2;
    
    GuessTcOutput=GuessTc;
    %mixed convection
    Pcon=pi.*Nueff.*Lambdaf.*(GuessTc-Ta);
    PconPrime=pi.*(Nueff.*Lambdaf+(GuessTc-Ta).*(LambdafPrime.*Nueff+NueffPrimedtc.*Lambdaf));
    PconPrimePrime=pi.*(Nueff.*LambdafPrime+NueffPrimedtc.*Lambdaf+...
        LambdafPrime.*Nueff+NueffPrimedtc.*Lambdaf+...
        (GuessTc-Ta).*(LambdafPrime.*NueffPrimedtc+...
        NueffPrimePrimedtc2.*Lambdaf+NueffPrimedtc.*LambdafPrime));

    %%                        MISMATCH AND UPDATE                        %%
    Hprime=I2Rprime-PradPrime-PconPrime;
    if(any(Hprime>0) || any(isnan(PconPrime)))
        msg='Hprime greater than zero or PconPrime is nan';
        error(msg);
    end
    
    Mismatch=I2R+Psol.*D.*alphas-Prad-Pcon;
    update=Mismatch./Hprime;
    GuessTcOutput=GuessTcOutput-update;  
    if(any(~isreal(GuessTcOutput)))
        error('Temperature process non-real')
    end
end