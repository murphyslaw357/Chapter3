function [GuessTc,I2R,I2Rprime,Prad,PradPrime,Pcon,PconPrime,GrPr,A,m,Nueff,Cinv,ninv,Reeff,C,n] =GetTempNewton(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,alphas,Psol,GuessTc)
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
    gplim1=1e-10;
    gplim2=1e-4;
    gplim3=1e-1;
    gplim4=1e2;
    gplim5=1e4;
    gplim6=1e7;
    gplim7=1e12;
    
    %[GuessTc]=GetGuessTemp(I,Ta,D,phi,Vw,alpha,beta,epsilons,alphas,Psol,mdl);
    Tolerance=1e-6; %tolerance criteria for Newton's method
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
        if(counter==1)
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
        end
        Nun=A*GrPr^m;
        %%
        NudfPrimedgrpr=m*A*GrPr^(m-1);
        NudfPrimedtc=NudfPrimedgrpr.*(Gr.*PrPrime+GrPrime.*Pr);
    
        %Mixed convection
        Re=(sin(phi)*Vw*D)./vf;  
        RePrimedtc=-((sin(phi)*Vw*D).*vfPrime)./(vf.^2);
        if(counter==1)
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
        end
        
        Ren=(Nun/Cinv)^(1/ninv);
        RenPrimednudf=(1/ninv)*(Nun/Cinv)^(1/ninv-1);      %differentiate(fNuRe,Nudf);
        RenPrimedtc=RenPrimednudf.*NudfPrimedtc;
        Reeff=sqrt((Re.^2)+(Ren.^2));
        Top=Re.*RePrimedtc+Ren.*RenPrimedtc;
        ReeffPrimedtc=Top./Reeff;
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
        Nueff=C*(Reeff^n);
            
            
%         Nueff=fReNu(Reeff);
        NueffPrimedreeff=n*C*(Reeff^(n-1));       %differentiate(fReNu,Reeff);
        NueffPrimedtc=NueffPrimedreeff.*ReeffPrimedtc;
        
        %Forced convection
%         Nu=fReNu(Re);  
%         NuPrimedre=differentiate(fReNu,Re);
%         NuPrimedtc=NuPrimedre.*RePrimedtc;
%         [row,~]=size(GuessTc);
%         Pcon=zeros(row,1);
%         PconPrime=zeros(row,1);
        
%         if(Vw==0)
%             %pure natural convection         
%             Pcon=pi*Nudf*Lambdaf*(GuessTc-Ta);
%             PconPrime=pi*(Nudf*Lambdaf+(GuessTc-Ta)*(LambdafPrime*Nudf+NudfPrimedtc*Lambdaf));
%         elseif(round(GuessTc,4)~=round(Ta,4))
            %mixed convection
            Pcon=pi*Nueff*Lambdaf*(GuessTc-Ta);
            PconPrime=pi*(Nueff*Lambdaf+(GuessTc-Ta)*(LambdafPrime*Nueff+NueffPrimedtc*Lambdaf));
%         else
%             %pure forced    
%             Pcon=pi*Nu*Lambdaf*(GuessTc-Ta);
%             PconPrime=pi*(Nu*Lambdaf+(GuessTc-Ta)*(LambdafPrime*Nu+NuPrimedtc*Lambdaf));
%         end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MISMATCH AND UPDATE%%%%%%%%%%%%%%
        Hprime=I2Rprime-PradPrime-PconPrime;
        if(Hprime>0 || isnan(PconPrime))
            msg='error condition';
            error(msg);
        end
        
        Mismatch=I2R+Psol*D*alphas-Prad-Pcon;
        update=Mismatch/Hprime;
        GuessTc=GuessTc-update; 

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
    end
end