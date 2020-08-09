function [GuessTc,I2R,dI2R_dTc,Prad,dPrad_dTc,Pcon,dPcon_dTc,GrPr,A,m,...
    Nueff,Cstar,nstar,Reeff,C,n] = ...
    GetTempNewton(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,alphas,Psol,GuessTc)
    %%%%%%%%%%%%%%%%%%%%%%%%%Input Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %I - RMS steady-state load current - amps
    %Ta - ambient temperature - degc
    %H - conductor elevation - meters
    %D - conductor diameter - meters
    %phi - angle between the wind direction and conductor axis - radians
    %Vw - Wind velocity - m/s
    %alpha - resistance at 0 degc
    %beta - resisance per degc
    %alphas - solar absorptivity
    %Psol - solar heating - w/m  
    %GuessTc - starting temperature guess
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    counter=0;
    sigmab=5.6697e-8;
    g=9.805;
    grprlim1=1e-10;
    grprlim2=1e-4;
    grprlim3=1e-1;
    grprlim4=1e2;
    grprlim5=1e4;
    grprlim6=1e7;
    grprlim7=1e12;
    
    Tolerance=1e-6; %tolerance criteria for Newton's method
    update=realmax;
    PrPrime=-1.25e-4;
    vfPrime=(4.75e-8)/((1-((6.5e-3)*H)/288.16)^5.2561);

    IIstar=abs(I)^2;
    dI2R_dTc=beta*IIstar;

    while(abs(update)>Tolerance)
        counter=counter+1;
        Tak = Ta + 273;
        GuessTck = GuessTc + 273;
        %%%%%%%%%%%%%%%%%%%%%%%RADIATIVE COOLING/HEATING%%%%%%%%%%%%%%%%%%%
        Prad=pi*D*sigmab*epsilons*((GuessTck^4)-(Tak^4));
        dPrad_dTc = 4*pi*D*sigmab*epsilons*(GuessTck^3);
        %%%%%%%%%%%%%%%%%%%%%%%%%JOULE HEATING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Resistance=alpha+beta*GuessTc;
        if(Resistance<0)
            error('Error - negative resistance')
        end
        I2R=IIstar*Resistance;
        %%%%%%%%%%%%%%%%%%%%%%%%CONVECTIVE COOLING/HEATING%%%%%%%%%%%%%%%%%
        GuessTfilm=(GuessTc+Ta)/2;
        GuessTfilmk = GuessTfilm + 273;
        vf=((1.32e-5)+(9.5e-8)*GuessTfilm)/...
            ((1-((6.5e-3)*H)/288.16)^5.2561);
        Gr=(g*(D^3).*abs(GuessTc-Ta))./(GuessTfilmk.*(vf.^2));  
        dGr_dTc = g*(D^3).*(Tak.*vf-(GuessTck.^2-Ta.^2).*vfPrime)./...
            ((vf.^3).*(GuessTfilmk.^2)); 
        Lambdaf = (2.42e-2)+(7.2e-5)*GuessTfilm;
        dLambdaf_dTc = 3.6e-5;
        Pr=0.715-(2.5e-4)*GuessTfilm;
        GrPr=Gr*Pr;
        %%Natural convection
        if(counter==1)
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
        end
        Nun=A*GrPr^m;
        dNun_dGrPr=m*A*GrPr^(m-1);
        dNun_dTc=dNun_dGrPr.*(Gr.*PrPrime+dGr_dTc.*Pr);
   
        if(counter==1)
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
                error('GrPr out of bounds')
            end
        end
        
        Ren=(Nun/Cstar)^(1/nstar);
        dRen_dNun=(1/(Cstar*nstar))*(Nun/Cstar)^(1/nstar-1);    
        dRen_dTc=dRen_dNun.*dNun_dTc;
        
        %% Forced Convection
        Ref=(sin(phi)*Vw*D)./vf;  
        dRef_dTc=-((sin(phi)*Vw*D).*vfPrime)./(vf.^2);
        
        %% Mixed convection
        Reeff = sqrt((Ref.^2)+(Ren.^2));
        Num = Ref.*dRef_dTc+Ren.*dRen_dTc;
        dReeff_dTc = Num./Reeff;
        if(counter==1)
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
        end
        Nueff = C*(Reeff^n);
        dNueff_dReeff = n*C*(Reeff^(n-1));
        dNueff_dTc = dNueff_dReeff.*dReeff_dTc;
        
        Pcon=pi*Nueff*Lambdaf*(GuessTc-Ta);
        dPcon_dTc=pi*Nueff*Lambdaf+pi*(GuessTc-Ta)*...
            (dLambdaf_dTc*Nueff+dNueff_dTc*Lambdaf);       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%MISMATCH AND UPDATE%%%%%%%%%%%%%%%%%%%%
        dH_dTc = dI2R_dTc-dPrad_dTc-dPcon_dTc;
        if(dH_dTc>0 || isnan(dPcon_dTc))
            error('Error - positive hPrime or dPcon is Nan');
        end
        
        Mismatch = I2R+Psol*D*alphas-Prad-Pcon;
        update = Mismatch/dH_dTc;
        GuessTc = GuessTc-update; 

        if(counter>=5000)
            GuessTc=nan;
            break;
        end
    end
    
    if(isnan(GuessTc))
        error('Converge failed! GuessTc is Nan');
    elseif(round(I2R,1)<0 || round(GuessTc-Ta,1)<0)
       error('Error condition! GuessTc too close to Ta')
    end
end