function [GuessTcOutput,I2R,dI2R_dTc,Prad,dPrad_dTc,d2Prad_dTc2,Pcon,...
    dPcon_dTc,d2Pcon_dTc2,Gr,dGr_dTc,Nun,A,m,Cstar,nstar,C,n] = ...
    GetTempNewtonFirstIteration2(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,...
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
    GrPr=Gr.*Pr;
    
    A=AmCinvninvCn(1);
    m=AmCinvninvCn(2);
    Nun=A.*GrPr.^m;
    dNun_dGrPr=(m*A).*GrPr.^(m-1);
    d2Nun_dGrRr2=((m-1)*m*A).*GrPr.^(m-2);
    dNun_dTc=dNun_dGrPr.*(Gr.*dPr_dTc+dGr_dTc.*Pr);
    d2Nun_dTc2=d2Nun_dGrRr2.*(dGr_dTc.*Pr+Gr.*dPr_dTc).^2+...
        dNun_dGrPr.*(2.*dGr_dTc.*dPr_dTc+d2Gr_dTc2.*Pr);
    Cstar=AmCinvninvCn(3);
    nstar=AmCinvninvCn(4);
    Ren=(Nun./Cstar).^(1/nstar);
    dRen_dNun=(1/(Cstar*nstar)).*(Nun./Cstar).^(1/nstar-1);    
    d2Ren_dNun2=((1-nstar)/((Cstar^2)*(nstar^2))).*(Nun./Cstar).^(1/nstar-2); 
    dRen_dTc=dRen_dNun.*dNun_dTc;
    d2Ren_dTc2=dRen_dNun.*d2Nun_dTc2+dNun_dTc.*d2Ren_dNun2;
    
    %% Forced convection
    Ref=(sin(phi)*Vw*D)./vf;  
    dRef_dTc=-((sin(phi)*Vw*D).*dvf_dTc)./(vf.^2);
    d2Ref_dTc2=(2*sin(phi)*Vw*D).*(dvf_dTc.^2)./(vf.^3);
    
    %% Mixed convection
    Reeff=sqrt((Ref.^2)+(Ren.^2));
    Num=Ref.*dRef_dTc+Ren.*dRen_dTc;
    dNum=Ren.*d2Ren_dTc2+dRen_dTc.^2+dRef_dTc.^2+Ref.*d2Ref_dTc2;
    
    dReeff_dTc=Num./Reeff;
    d2Reeff_dTc2=(Reeff.*dNum-Num.*dReeff_dTc)./(Reeff.^2);
    C=AmCinvninvCn(5);
    n=AmCinvninvCn(6);
    Nueff = C.*Reeff.^n;
    dNueff_dReeff=(C*n).*(Reeff.^(n-1));
    d2Nueff_dReeff2=((n-1)*C*n).*(Reeff.^(n-2));
    dNueff_dTc=dNueff_dReeff.*dReeff_dTc;
    d2Nueff_dTc2=d2Nueff_dReeff2.*dReeff_dTc+dNueff_dReeff.*d2Reeff_dTc2;
    
    Pcon=pi.*Nueff.*Lambdaf.*(GuessTc-Ta);
    dPcon_dTc=pi.*(Nueff.*Lambdaf+(GuessTc-Ta).*...
        (dLambdaf_dTc.*Nueff+dNueff_dTc.*Lambdaf));
    d2Pcon_dTc2=pi.*(d2Nueff_dTc2.*Lambdaf.*(GuessTc-Ta)+...
        2.*dNueff_dTc.*(dLambdaf_dTc.*(GuessTc-Ta)+Lambdaf)+...
        2.*Nueff.*dLambdaf_dTc);
    
    %%                        MISMATCH AND UPDATE                        %%
    Hprime=dI2R_dTc-dPrad_dTc-dPcon_dTc;
    if(any(Hprime>0) || any(isnan(dPcon_dTc)))
        error('Hprime greater than zero or dPcon_dTc is nan');
    end
    
    Mismatch=I2R+Psol.*D.*alphas-Prad-Pcon;
    update=Mismatch./Hprime;
    GuessTcOutput=GuessTc-update;  
    if(any(~isreal(GuessTcOutput)))
        error('Temperature process non-real')
    end
end