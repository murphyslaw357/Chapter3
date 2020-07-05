function [GuessTcs,Prads,Pcons,I2R] =GetTempIEEEBisection(Is,Tas,H,D,phi,Vws,alpha,beta,epsilons,alphas,Psols,fGrPr,fReNu,fNuRe)
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
    weatherPermutationCount = size(Tas,1);
    GuessTcLeft = Tas;
    GuessTcRight = 1000.*ones(weatherPermutationCount,1);
    Tolerance=0.0001; %tolerance criteria for bisection method
    IIstar=abs(Is).^2;
        
    while(any(abs(GuessTcRight-GuessTcLeft)>Tolerance))
        counter=counter+1;
        GuessTcs = (GuessTcLeft+GuessTcRight)./2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RADIATIVE COOLING/HEATING%%%%%%%%%%%
        Prads=(pi*D*sigmab*epsilons).*(((GuessTcs+273).^4)-((Tas+273).^4));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%JOULE HEATING%%%%%%%%%%%%%%%%%%%%%%%
        Resistance=beta.*GuessTcs+alpha;
        %if(Resistance<0)
        %    Resistance=0;
        %end
        I2R=IIstar.*Resistance;
        %%%%%%%%%%%%%%%%%%%%%%%%CONVECTIVE COOLING/HEATING%%%%%%%%%%%%%%%%%
        GuessTfilm=(GuessTcs+Tas)./2;
        vf=((1.32e-5)+(9.5e-8).*GuessTfilm).*((1-((6.5e-3)*H)/288.16)^-5.2561);
        Tfilmk=GuessTfilm+273;
        Lambdaf=(2.42e-2)+(7.2e-5).*GuessTfilm;
        Gr=((g*(D^3)).*abs(GuessTcs-Tas))./(Tfilmk.*(vf.^2));  
        Pr=0.715-(2.5e-4).*GuessTfilm;
        GrPr=Gr.*Pr;
        %Natural convection
        Nun=fGrPr(GrPr);
        %Mixed convection
        Ref=((sin(phi)*D).*Vws)./vf;  
        Ren=fNuRe(Nun);
        Reeff=sqrt((Ref.^2)+(Ren.^2));
        Nueff=fReNu(Reeff);
        Pcons=pi.*Nueff.*Lambdaf.*(GuessTcs-Tas);
        %Check mismatch and update
        
        %if((I2R + Psol*D)>(Prad+Pcon))
            GuessTcLeft((I2R + Psols.*D.*alphas)>(Prads+Pcons)) = GuessTcs((I2R + Psols.*D.*alphas)>(Prads+Pcons));
        %else
            GuessTcRight((I2R + Psols.*D.*alphas)<=(Prads+Pcons)) = GuessTcs((I2R + Psols.*D.*alphas)<=(Prads+Pcons));
        %end
        %a=[length(I2R),length(Psols),length(Prad),length(Pcon),length(GuessTcLeft),length(GuessTcRight),length(GuessTcs)]
        %close all
        %plot(GuessTcLeft)
        %hold on
        %plot(GuessTcRight)
        %plot(GuessTcRight-GuessTcLeft)
        %hold on
        %plot(GuessTcLeft)
        if(counter>=5000)
            GuessTcs=nan(weatherPermutationCount);
            break;
        end
    end
    
    if(any(isnan(GuessTcs)))
        msg='converge failed!';
        error(msg);
    elseif(any(round(I2R,1)<0) || any(round(GuessTcs-Tas,1)<0))
       msg='error condition!';
       error(msg)
    elseif(any(abs(I2R + Psols.*D.*alphas - Prads - Pcons)>0.01))
       msg='error condition!';
       error(msg)
    end
end