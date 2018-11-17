function [GuessTc,GuessTfilm,Pcon,Prad,I2R] =GetTempNewton(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,Psol)
    phi=12*pi/180;
    %Vw=0.1;
    %Psol=Psol*0.5;
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
    IIstar=abs(I)^2;
    GuessTc=((Psol+IIstar*(alpha+25*beta))/(pi*D*sigmab*epsilons)+((Ta+273)^4))^(1/4)-273;   
    Tolerance=0.05; %tolerance criteria for Newton's method
    epsilon=1e-9;
    update=realmax;
    while(abs(update)>Tolerance)
        counter=counter+1;
        GuessTfilm=(GuessTc+Ta)/2;
        vf=((1.32e-5)+(9.5e-8)*GuessTfilm)*((1-((6.5e-3)*H)/288.16)^-5.2561);
        vfprime=(4.75e-8)*((1-((6.5e-3)*H)/288.16)^-5.2561);
        
        Tfilmk=GuessTfilm+273;
        Tfilmkprime=1/2;
        Gr=(g*(D^3)*(GuessTc-Ta))/(Tfilmk*(vf^2));        
        Grprime=(Tfilmk*(vf^2)*g*(D^3)-g*(D^3)*(GuessTc-Ta)*(Tfilmkprime*(vf^2)+Tfilmk*2*vf*vfprime))/((Tfilmk^2)*(vf^4));
        
        Pr=0.715-(2.5e-4)*GuessTfilm;
        Prprime=-1.25e-4;
        GrPr=Gr*Pr;
        %lookup A and m
        A=0;
        m=0;
        if(GrPr>1e-10 && GrPr<=1e-4)
            A=0.675;
            m=0.058;
        elseif(Gr*Pr>1e-4 && Gr*Pr<=1e-1)
            A=0.889;
            m=0.088;
        elseif(Gr*Pr>1e-1 && Gr*Pr<=1e2)
            A=1.02;
            m=0.148;
        elseif(Gr*Pr>1e2 && Gr*Pr<=1e4)
            A=0.85;
            m=0.188;
        elseif(Gr*Pr>1e4 && Gr*Pr<=1e7)
            A=0.48;
            m=0.25;
        elseif(Gr*Pr>1e7 && Gr*Pr<=1e12)
            A=0.125;
            m=0.333;
        end
        Re=sin(phi)*Vw*D/vf;
        Reprime=-(sin(phi)*Vw*D)*(4.75e-8)*((1-((6.5e-3)*H)/288.16)^-5.2561)/(vf^2);
        %lookup C and n
        C=0;
        n=0;
        if(Re>1e-4 && Re<=4e-3)
            C=0.437;
            n=0.0895;
        elseif(Re>4e-3 && Re<=9e-2)
            C=0.565;
            n=0.136;
        elseif(Re>9e-2 && Re<=1)
            C=0.8;
            n=0.28;
        elseif(Re>1 && Re<=35)
            C=0.795;
            n=0.384;
        elseif(Re>35 && Re<=5e3)
            C=0.583;
            n=0.471;
        elseif(Re>5e3 && Re<=5e4)
            C=0.148;
            n=0.633;
        elseif(Re>5e4)% && Re<2e5)
            C=0.0208;
            n=0.814;
        end
        %Tsky=(0.0552*(Ta+273)^1.5)-273;
        %Tg=Ta+2;
               
        Reeq=sqrt((Re^2)+(((A/C)*(GrPr^m))^(2/n)));
       
        %%
        %re-lookup C and n for effective Reynolds number
        C=0;
        n=0;
        if(Reeq>1e-4 && Reeq<=4e-3)
            C=0.437;
            n=0.0895;
        elseif(Reeq>4e-3 && Reeq<=9e-2)
            C=0.565;
            n=0.136;
        elseif(Reeq>9e-2 && Reeq<=1)
            C=0.8;
            n=0.28;
        elseif(Reeq>1 && Reeq<=35)
            C=0.795;
            n=0.384;
        elseif(Reeq>35 && Reeq<=5e3)
            C=0.583;
            n=0.471;
        elseif(Reeq>5e3 && Reeq<=5e4)
            C=0.148;
            n=0.633;
        elseif(Reeq>5e4)% && Re<2e5)
            C=0.0208;
            n=0.814;
        end
        
        I2R=IIstar*(beta*GuessTc+alpha);
        I2Rprime=beta*IIstar;
        
        %Tsky=(0.0552*(Ta+273)^1.5)-273;
        %Tg=Ta+2;
        Prad=pi*D*sigmab*epsilons*(((GuessTc+273)^4)-((Ta+273)^4));
        %Prad=pi*D*sigmab*epsilons*(((GuessTc+273)^4)-0.5*((Tsky+273)^4)-0.5*((Tg+273)^4));
        Pradprime=4*pi*D*sigmab*epsilons*(GuessTc+273)^3;
        
        Nueff=C*Reeq^n;
        Nueffprime=C*(n/2)*(((Re^2)+(((A/C)*(GrPr^m))^(2/n)))^((n/2)-1))*...
            (2*Re*Reprime+((A/C)^(2/n))*(2*m/n)*(GrPr^((2*m/n)-1))*(Grprime*Pr+Gr*Prprime));        
        
        LambdafTcTa=((2.42e-2)+(7.2e-5)*GuessTfilm)*(GuessTc-Ta);
        LambdafTcTaPrime=((2.42e-2)+(7.2e-5)*GuessTfilm)+(3.6e-5)*(GuessTc-Ta);

        Pcon=pi*Nueff*LambdafTcTa;
        Pconprime=pi*(Nueffprime*LambdafTcTa+Nueff*LambdafTcTaPrime);
%         if(GuessTc>Ta)

             Hprime=I2Rprime-Pradprime-Pconprime;
            if(abs(Hprime)<epsilon) 
                GuessTc=nan; %denominator is too small- exit loop
                break;
            end
             Mismatch=I2R+Psol-Prad-Pcon;
             update=Mismatch/Hprime;
             GuessTc=GuessTc-update;
%         else
%             Hprime=-Pradprime;
%             Mismatch=I2R+Psol-Prad;
%             GuessTc=GuessTc-Mismatch/Hprime;
%         end
        %disp([GuessTc Mismatch])
        if(counter>=500)
            GuessTc=nan;
            break;
        end
    end
    if(isnan(GuessTc))
        msg='converge failed!';
        disp(msg);
    elseif(I2R<0 || Psol <0 || Prad<0 || Pcon<0 || GuessTc<Ta)
       msg='error condition!';
       error(msg)
    end
    %disp(counter)
    %GuessTc=GuessTc*1.5;
    %GuessTc=GuessTc+10;
end