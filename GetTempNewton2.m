function [GuessTc,I2R,I2Rprime,Prad,PradPrime,Pcon,PconPrime] =GetTempNewton2(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,Psol,fGrPr,fReNu,fNuRe)
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
%     IIstar=abs(I)^2;
    [GuessTc,Pcon2,Prad2,Pj2,Nueff2,lambdaf2,Nre]=GetGuessTemp(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,Psol);

    Tolerance=0.001; %tolerance criteria for Newton's method
    update=realmax;

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
%         GrPrime=g*(D^3)*(Tfilmk*(vf)-abs(GuessTc-Ta)*(TfilmkPrime*vf+Tfilmk*2*vfPrime))/...
%                 ((Tfilmk^2)*(vf^3));
        Pr=0.715-(2.5e-4)*GuessTfilm;
        GrPr=Gr*Pr;
        if(GuessTc~=Ta && Vw==0)
            %pure natural convection 
            Nudf=fGrPr(GrPr);
            NudfPrime=differentiate(fGrPr,GrPr); 
            Pcon=pi*Nudf*Lambdaf*(GuessTc-Ta);
            PconPrime=pi*(Nudf*Lambdaf+LambdafPrime*Nudf*(GuessTc-Ta)+NudfPrime*Lambdaf*(GuessTc-Ta));
        elseif(GuessTc~=Ta && Vw ~=0)
        %mixed convection
            Nudf=fGrPr(GrPr);
            Req=fNuRe(Nudf);
%             %lookup C and n
%             C=0;
%             n=0;
%             if(Nudf<=Nulim1)
%                 C=0.437;
%                 n=0.0895;
%             elseif(Nudf>Nulim1 && Nudf<=Nulim2)
%                 C=0.565;
%                 n=0.136;
%             elseif(Nudf>Nulim2 && Nudf<=Nulim3)
%                 C=0.8;
%                 n=0.28;
%             elseif(Nudf>Nulim3 && Nudf<=Nulim4)
%                 C=0.795;
%                 n=0.384;
%             elseif(Nudf>Nulim4 && Nudf<=Nulim5)
%                 C=0.583;
%                 n=0.471;
%             elseif(Nudf>Nulim5 && Nudf<=Nulim6)
%                 C=0.148;
%                 n=0.633;
%             elseif(Nudf>Nulim6)
%                 C=0.0208;
%                 n=0.814;
%             end
%             if(C==0)
%                 msg='error';
%                 error(msg)
%             end
%             Req=(Nudf/C)^(1/n);
            Re=sin(phi)*Vw*D/vf;
            Reeff=sqrt((Re^2)+(Req^2));
            NueffPrime=differentiate(fReNu,Reeff);
            Nueff=fReNu(Reeff);
            Pcon=pi*Nueff*Lambdaf*(GuessTc-Ta);
            PconPrime=pi*(Nueff*Lambdaf+LambdafPrime*Nueff*(GuessTc-Ta)+NueffPrime*Lambdaf*(GuessTc-Ta));
        elseif(GuessTc==Ta && Vw~=0)
        %pure forced    
            Re=sin(phi)*Vw*D/vf;
            NuPrime=differentiate(fReNu,Re);
            Nu=fReNu(Re);
            Pcon=pi*Nu*Lambdaf*(GuessTc-Ta);
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