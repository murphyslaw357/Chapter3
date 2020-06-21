function [GuessTc,I2R,Prad,Pcon] =GetTempIEEEBisection(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,Psol,fGrPr,fReNu,fNuRe)
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
    phi=pi/2;
    GuessTcLeft = Ta;
    GuessTcRight = 1000;
    Tolerance=0.0001; %tolerance criteria for bisection method

    IIstar=abs(I)^2;

    while(abs(GuessTcRight-GuessTcLeft)>Tolerance)
        counter=counter+1;
        GuessTc=(GuessTcLeft+GuessTcRight)/2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RADIATIVE COOLING/HEATING%%%%%%%%%%%
        Prad=pi*D*sigmab*epsilons*(((GuessTc+273)^4)-((Ta+273)^4));
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
        Lambdaf=(2.42e-2)+(7.2e-5)*GuessTfilm;
        Gr=(g*(D^3).*abs(GuessTc-Ta))./(Tfilmk.*(vf.^2));  
        Pr=0.715-(2.5e-4)*GuessTfilm;
        GrPr=Gr*Pr;
        %Natural convection
        Nun=fGrPr(GrPr);
        %Mixed convection
        Ref=(sin(phi)*Vw*D)./vf;  
        Ren=fNuRe(Nun);
        Reeff=sqrt((Ref.^2)+(Ren.^2));
        Nueff=fReNu(Reeff);
        Pcon=pi*Nueff*Lambdaf*(GuessTc-Ta);
        %Check mismatch and update
        if((I2R + Psol*D)>(Prad+Pcon))
            GuessTcLeft = GuessTc;
        else
            GuessTcRight = GuessTc;
        end
        
        if(counter>=5000)
            GuessTc=nan;
            break;
        end
    end
    
    if(isnan(GuessTc))
        msg='converge failed!';
        error(msg);
    elseif(round(I2R,1)<0 || round(GuessTc-Ta,1)<0)
       msg='error condition!';
       error(msg)
    elseif(abs(I2R + Psol*D - Prad - Pcon)>0.01)
       msg='error condition!';
       error(msg)
    end
end