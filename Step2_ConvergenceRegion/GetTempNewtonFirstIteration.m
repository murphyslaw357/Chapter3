function [GuessTcOutput,I2R,dI2R_dTc,Prad,dPrad_dTc,d2Prad_dTc2,Pcon,...
    dPcon_dTc,d2Pcon_dTc2,Gr,dGr_dTc,Nun,A,m,Cstar,nstar,C,n] = ...
    GetTempNewtonFirstIteration(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,...
    alphas,Psol,GuessTc,AmCinvninvCn)
    %%%%%%%%%%%%%%%%%%%%%%%%%Input Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %I - RMS steady-state load current - amps
    %Ta - ambient temperature - degc
    %H - conductor elevation - meters
    %D - conductor diameter - meters
    %phi - angle between the wind direction and conductor axis - radians
    %Vw - Wind velocity - m/s
    %alpha - resistance at 0 degc
    %beta - resisance per degc
    %epsilons - radiative absorptivity
    %alphas - solar absorptivity
    %Psol - solar heating - w/m  
    %GuessTc - starting temperature guess
    %AmCinvninvCn - convective cooling parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sigmab=5.6697e-8;
    g=9.805;
    grprlim1=1e-10;
    grprlim2=1e-4;
    grprlim3=1e-1;
    grprlim4=1e2;
    grprlim5=1e4;
    grprlim6=1e7;
    grprlim7=1e12;
    
    IIstar=abs(I).^2;
    dI2R_dTc=beta.*IIstar;
    
    Tak = Ta + 273;
    GuessTck = GuessTc + 273;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%RADIATIVE COOLING/HEATING%%%%%%%%%%%%%%%%%%%
    Prad = (pi*D*sigmab*epsilons).*((GuessTck.^4)-(Tak.^4));
    dPrad_dTc = (4*pi*D*sigmab*epsilons).*(GuessTck.^3);
    d2Prad_dTc2 = (12*pi*D*sigmab*epsilons).*(GuessTck.^2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%JOULE HEATING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Resistance=alpha+beta.*GuessTc;
    if(any(Resistance<0))
        error('Error - negative resistance')
    end
    I2R=IIstar.*Resistance;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%CONVECTIVE COOLING/HEATING%%%%%%%%%%%%%%%%%
    if(GuessTc<Ta)
        error('Error - starting guess less than Ta')
    end
    dPr_dTc=-1.25e-4;
    GuessTfilm=(GuessTc+Ta)./2;
    GuessTfilmk=GuessTfilm+273;
    vf=((1.32e-5)+(9.5e-8).*GuessTfilm)./...
            ((1-((6.5e-3).*H)./288.16).^5.2561);
    dvf_dTc=(4.75e-8)./((1-((6.5e-3).*H)./288.16).^5.2561);
    gD3 = g*(D^3);
    %% Natural convection
    Gr = gD3.*abs(GuessTc-Ta)./(GuessTfilmk.*(vf.^2));
    dGr_dTc = gD3.*(Tak.*vf-(GuessTck.^2-Tak.^2).*dvf_dTc)./...
            ((vf.^3).*(GuessTfilmk.^2)); 
    d2Gr_dTc2 = gD3.*(Tak-2.*GuessTck).*dvf_dTc./((GuessTfilmk.^2).*(vf.^3))-...
        (gD3./((GuessTfilmk.^3).*(vf.^4))).*(vf+3.*GuessTfilmk.*dvf_dTc).*...
        (Tak.*vf-((GuessTck.^2)-(Tak.^2)).*dvf_dTc);
    Lambdaf=(2.42e-2)+(7.2e-5).*GuessTfilm;
    dLambdaf_dTc=3.6e-5;    
    Pr=0.715-(2.5e-4).*GuessTfilm;
    GrPr=Gr*Pr;
    
    if(isempty(AmCinvninvCn))
        if(GrPr>grprlim1 && GrPr<=grprlim2)
            A=0.675;
            m=0.058;
        elseif(GrPr>grprlim2 && GrPr<=grprlim3)
            A=0.889;
            m=0.088;
        elseif(GrPr>grprlim3 && GrPr<=grprlim4)
            A=1.02;
            m=0.148;
        elseif(GrPr>grprlim4 && GrPr<=grprlim5)
            A=0.85;
            m=0.188;
        elseif(GrPr>grprlim5 && GrPr<=grprlim6)
            A=0.48;
            m=0.25;
        elseif(GrPr>grprlim6 && GrPr<=grprlim7)
            A=0.125;
            m=0.333;
        end
    else
        A=AmCinvninvCn(1);
        m=AmCinvninvCn(2);
    end
    Nun=A*GrPr^m;
    dNun_dGrPr=m*A*GrPr^(m-1);
    d2Nun_dGrRr2=(m-1)*m*A*GrPr^(m-2);
    dNun_dTc=dNun_dGrPr*(Gr*dPr_dTc+dGr_dTc*Pr);
    d2Nun_dTc2=d2Nun_dGrRr2*(dGr_dTc*Pr+Gr*dPr_dTc)^2+...
        dNun_dGrPr*(2*dGr_dTc*dPr_dTc+d2Gr_dTc2*Pr);
    if(isempty(AmCinvninvCn))
        if(Nun>0.1916 && Nun<=0.2666)
            Cstar=0.437;
            nstar=0.0895;
        elseif(Nun>0.2666 && Nun<=0.40685)
            Cstar=0.565;
            nstar=0.136;
        elseif(Nun>0.40685 && Nun<=0.8136)
            Cstar=0.8;
            nstar=0.28;
        elseif(Nun>0.8136 && Nun<=3.1253)
            Cstar=0.795;
            nstar=0.384;
        elseif(Nun>3.1253 && Nun<=31.45)
            Cstar=0.583;
            nstar=0.471;
        elseif(Nun>31.45 && Nun<=141.449)
            Cstar=0.148;
            nstar=0.633;
        elseif(Nun>141.449 && Nun<=429.6371)
            Cstar=0.0208;
            nstar=0.814;
        else
            error('out of bounds')
        end  
    else
        Cstar=AmCinvninvCn(3);
        nstar=AmCinvninvCn(4);
    end
    Ren=(Nun/Cstar)^(1/nstar);
    dRen_dNun=(1/(Cstar*nstar))*(Nun/Cstar)^(1/nstar-1);    
    d2Ren_dNun2=((1-nstar)/((Cstar^2)*(nstar^2)))*(Nun/Cstar)^(1/nstar-2); 
    dRen_dTc=dRen_dNun.*dNun_dTc;
    d2Ren_dTc2=dRen_dNun.*d2Nun_dTc2+dNun_dTc.*d2Ren_dNun2;
    
    %% Forced convection
    Ref=(sin(phi)*Vw*D)/vf;  
    dRef_dTc=-((sin(phi)*Vw*D)*dvf_dTc)/(vf^2);
    d2Ref_dTc2=2*(sin(phi)*Vw*D)*(dvf_dTc.^2)./(vf.^3);
    
    %% Mixed convection
    Reeff=sqrt((Ref.^2)+(Ren.^2));
    Num=Ref.*dRef_dTc+Ren.*dRen_dTc;
    dNum=Ren.*d2Ren_dTc2+dRen_dTc.^2+dRef_dTc.^2+Ref.*d2Ref_dTc2;
    
    dReeff_dTc=Num./Reeff;
    d2Reeff_dTc2=(Reeff.*dNum-Num.*dReeff_dTc)./(Reeff.^2);
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
    Nueff = C*(Reeff^n);
    dNueff_dReeff=C*n*(Reeff^(n-1));
    d2Nueff_dReeff2=(n-1)*C*n*(Reeff^(n-2));
    dNueff_dTc=dNueff_dReeff.*dReeff_dTc;
    d2Nueff_dTc2=d2Nueff_dReeff2*dReeff_dTc+dNueff_dReeff*d2Reeff_dTc2;
    
    Pcon=pi*Nueff*Lambdaf*(GuessTc-Ta);
    dPcon_dTc=pi*(Nueff*Lambdaf+(GuessTc-Ta)*...
        (dLambdaf_dTc*Nueff+dNueff_dTc*Lambdaf));
    d2Pcon_dTc2=pi*(d2Nueff_dTc2*Lambdaf*(GuessTc-Ta)+...
        2*dNueff_dTc*(dLambdaf_dTc*(GuessTc-Ta)+Lambdaf)+...
        2*Nueff*dLambdaf_dTc);
    
    %%                        MISMATCH AND UPDATE                        %%
    Hprime=dI2R_dTc-dPrad_dTc-dPcon_dTc;
    if(Hprime>0 || isnan(dPcon_dTc))
        error('Hprime greater than zero or PconPrime is nan');
    end
    
    Mismatch=I2R+Psol*D*alphas-Prad-Pcon;
    update=Mismatch/Hprime;
    GuessTcOutput=GuessTc-update;  
    if(~isreal(GuessTcOutput))
        error('Temperature process non-real')
    end
end