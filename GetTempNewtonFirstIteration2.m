function [GuessTcOutput,I2R,I2Rprime,Prad,PradPrime,PradPrimePrime,Pcon,...
    PconPrime,PconPrimePrime,Gr,GrPrime,Nun] = ...
    GetTempNewtonFirstIteration2(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,...
    Psol,GuessTc)%,fGrPr,fReNu,fNuRe
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
    IIstar=abs(I)^2;
    %%                          RADIATIVE COOLING/HEATING                %%
    Prad=(pi*D*sigmab*epsilons).*(((GuessTc+273).^4)-((Ta+273).^4));
    PradPrime=(4*pi*D*sigmab*epsilons).*(GuessTc+273).^3;
    PradPrimePrime=(12*pi*D*sigmab*epsilons).*(GuessTc+273).^2;
    %%                          JOULE HEATING                            %%
    I2Rprime=beta*IIstar;
    Resistance=beta*GuessTc+alpha;
    if(Resistance<0)
        Resistance=0;
    end
    I2R=IIstar*Resistance;
    %%                      CONVECTIVE COOLING/HEATING                   %%
    [row,~]=size(GuessTc);
    GuessTcOutput=zeros(row,1);
    
    GuessTfilm=(GuessTc+Ta)./2;
    vf=((1.32e-5)+(9.5e-8).*GuessTfilm).*((1-((6.5e-3)*H)/288.16)^-5.2561);
    vfPrime=(4.75e-8)*((1-((6.5e-3)*H)/288.16)^-5.2561);
    vfPrimePrime=0;
    Tfilmk=GuessTfilm+273;
    Lambdaf=(2.42e-2)+(7.2e-5)*GuessTfilm;
    LambdafPrime=3.6e-5;
    
    gD3 = g*(D^3);
    Gr = gD3.*(GuessTc-Ta)./(Tfilmk.*(vf.^2));
    GrPrime = gD3.*(Ta.*vf-((GuessTc.^2)-Ta.^2).*vfPrime)./...
        ((vf.^3).*(Tfilmk.^2));
    GrPrimePrime=(-2.*gD3.*GuessTc./((Tfilmk.^2).*(vf.^3))).*vfPrime-...
        (gD3./((Tfilmk.^3).*(vf.^4))).*(vf+3.*Tfilmk.*vfPrime).*...
        (Ta.*vf-((GuessTc.^2)-Ta.^2).*vfPrime);   

    Pr=0.715-(2.5e-4).*GuessTfilm;
    PrPrime=-1.25e-4;
    GrPr = Gr.*Pr;
    GrPrPrimedTc = Gr.*PrPrime+GrPrime.*Pr;
    GrPrPrimePrimedTc2 = GrPrime.*PrPrime+GrPrimePrime.*Pr+GrPrime.*PrPrime;
    
    %% Forced convection
    Ref = (sin(phi)*Vw*D)./vf;  
    RefPrimedTc = -((sin(phi)*Vw*D).*vfPrime)./(vf.^2);
    RePrimePrimedTc2 = ((2*sin(phi)*Vw*D)./(vf.^3)).*(vfPrime.^2);
    %% Natural convection
    Nun=exp(1).^(0.5524.*nthroot(GrPr,9))-0.5798;
    
    NunPrimedTc=exp(1).^(0.5524.*nthroot(GrPr,9)).*(0.5524.*(1/9).*nthroot(GrPr.^-8,9)).*GrPrPrimedTc;
    NunPrimePrimedTc2=exp(1).^(0.5524.*nthroot(GrPr,9)).*((0.5524.*(1/9).*GrPr.^(-8/9)).*GrPrPrimedTc).^2+...
        exp(1).^(0.5524.*nthroot(GrPr,9)).*((-8/9).*0.5524.*(1/9).*GrPr.^(-17/9)).*GrPrPrimedTc.^2+...
        exp(1).^(0.5524.*nthroot(GrPr,9)).*(0.5524.*(1/9).*GrPr.^(-8/9)).*GrPrPrimePrimedTc2;
% % %     NudfPrimedtc=NudfPrimedgrpr.*(Gr.*PrPrime+GrPrime.*Pr);
% % %     NudfPrimePrimedtc2=NudfPrimePrimedgrpr2.*(Gr.*PrPrime+GrPrime.*Pr)+...
% % %         NudfPrimedgrpr.*...
% % %         (GrPrime.*PrPrime+GrPrimePrime.*Pr+GrPrime.*PrPrime);
    %% Mixed convection
    Reeq=(log(Nun+1.91)./1.055).^7;
    ReeqPrimedNun = 7.*((log(Nun+1.91)./1.055).^6).*(1./((Nun+1.91).*1.055));
    ReeqPrimePrimedNun2= 42.*((log(Nun+1.91)./1.055).^5).*((1./((Nun+1.91).*1.055)).^2)-...
        7.*((log(Nun+1.91)./1.055).^6).*(1./(((Nun+1.91).^2).*1.055));
    ReeqPrimedTc=ReeqPrimedNun.*NunPrimedTc;
    ReeqPrimePrimedTc2=ReeqPrimedNun.*NunPrimePrimedTc2+...
        ReeqPrimePrimedNun2.*(NunPrimedTc.^2);
    Reeff=sqrt((Ref.^2)+(Reeq.^2));
    Num=Ref.*RefPrimedTc+Reeq.*ReeqPrimedTc;
    dnum=RefPrimedTc.*RefPrimedTc+Ref.*RePrimePrimedTc2+ReeqPrimedTc.*...
        ReeqPrimedTc+Reeq.*ReeqPrimePrimedTc2;

    ReeffPrimedTc=Num./Reeff;
    ReeffPrimePrimedTc2=(Reeff.*dnum-Num.*ReeffPrimedTc)./(Reeff.^2);
    Nueff=exp(1).^(1.055.*nthroot(Reeff,7))-1.91;
    NueffPrimedReeff = exp(1).^(1.055.*nthroot(Reeff,7)).*((1.055/7).*Reeff.^(-6/7));
    NueffPrimedTc=NueffPrimedReeff.*ReeffPrimedTc;
    NueffPrimePrimedReeff2 = exp(1).^(1.055.*nthroot(Reeff,7)).*(((1.055/7).*Reeff.^(-6/7)).^2)+...
        exp(1).^(1.055.*nthroot(Reeff,7)).*((-6/7).*(1.055/7).*Reeff.^(-13/7));
    NueffPrimePrimedTc2=NueffPrimedReeff.*ReeffPrimePrimedTc2+...
        NueffPrimePrimedReeff2.*(ReeffPrimedTc.^2);

    Pcon=pi.*Nueff.*Lambdaf.*(GuessTc-Ta);
    PconPrime=pi.*(Nueff.*Lambdaf+(GuessTc-Ta).*(LambdafPrime.*Nueff+NueffPrimedTc.*Lambdaf));
    PconPrimePrime=pi.*(Nueff.*LambdafPrime+NueffPrimedTc.*Lambdaf+...
        LambdafPrime.*Nueff+NueffPrimedTc.*Lambdaf+...
        (GuessTc-Ta).*(LambdafPrime.*NueffPrimedTc+NueffPrimePrimedTc2.*Lambdaf+NueffPrimedTc.*LambdafPrime));
        
    %% Tie it all together
    for i=1:row
        %%                        MISMATCH AND UPDATE                    %%
        Hprime=I2Rprime-PradPrime(i)-PconPrime(i);
        if(Hprime>0 || isnan(PconPrime(i)))
            msg='error condition1';
            disp(msg);
        end

        Mismatch=I2R(i)+Psol*D-Prad(i)-Pcon(i);
        update=Mismatch/Hprime;
        GuessTcOutput(i)=GuessTc(i)-update;  
    end
end