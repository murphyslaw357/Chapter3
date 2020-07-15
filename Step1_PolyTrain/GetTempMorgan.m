function [GuessTc,I2R,Prad,Pcon,GrPr,A,m,Nueff,Cinv,ninv,Reeff,C,n] = GetTempMorgan(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,alphas,Psol)
    g=9.805;
    sigmab=5.6697e-8;
    tolerance=1e-6;
    GuessTcRight=1000;
    GuessTcLeft=Ta;
    IIstar=abs(I)^2;
    counter=0;
    while(abs(GuessTcRight-GuessTcLeft)>tolerance)
        counter=counter+1;
        GuessTc=(GuessTcRight+GuessTcLeft)/2;
        I2R=IIstar*(alpha+beta*GuessTc);
        Prad=pi*D*sigmab*epsilons*(((GuessTc+273)^4)-((Ta+273)^4));
        GuessTfilm=(GuessTc+Ta)/2;
        vf=((1.32e-5)+(9.5e-8)*GuessTfilm)*((1-((6.5e-3)*H)/288.16)^-5.2561);
        
        Gr=(g*(D^3).*abs(GuessTc-Ta))./((GuessTfilm+273).*(vf.^2));  
        Lambdaf=(2.42e-2)+(7.2e-5)*GuessTfilm;
        Pr=0.715-(2.5e-4)*GuessTfilm;
        GrPr=Gr*Pr;
        %then do the table 4.4 to find m and A
        
        gplim1=1e-10;
        gplim2=1e-4;
        gplim3=1e-1;
        gplim4=1e2;
        gplim5=1e4;
        gplim6=1e7;
        gplim7=1e12;
        
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
        
%         scatter(log([(gplim1+gplim2)/2 (gplim2+gplim3)/2 (gplim3+gplim4)/2 (gplim4+gplim5)/2 (gplim5+gplim6)/2 (gplim6+gplim7)/2]),[0.675 0.889 1.02 0.85 0.48 0.125])
%         figure
%         scatter([gplim1 gplim2 gplim3 gplim4 gplim5 gplim6 gplim7],[0.058 0.088 0.148 0.188 0.25 0.333 0.333])
        
        Nueff=A*(GrPr)^m;
        if(Nueff>0.1916 && Nueff<=0.2666)
            Cinv=0.437;
            ninv=0.0895;
        elseif(Nueff>0.2666 && Nueff<=0.40685)
           Cinv=0.565;
           ninv=0.136;
        elseif(Nueff>0.40685 && Nueff<=0.8136)
           Cinv=0.8;
           ninv=0.28;
        elseif(Nueff>0.8136 && Nueff<=3.1253)
           Cinv=0.795;
           ninv=0.384;
        elseif(Nueff>3.1253 && Nueff<=31.45)
           Cinv=0.583;
           ninv=0.471;
        elseif(Nueff>31.45 && Nueff<=141.449)
           Cinv=0.148;
           ninv=0.633;
        elseif(Nueff>141.449 && Nueff<=429.6371)
           Cinv=0.0208;
           ninv=0.814;
        else
            error('out of bounds')
        end
        Ren=(Nueff/Cinv)^(1/ninv);
        
        Re=sin(phi)*Vw*D/vf;
        Reeff=sqrt(Ren^2+Re^2);
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
        Nueff=C*(Reeff^n);
        Pcon=pi*Nueff*Lambdaf*(GuessTc-Ta);
        
        if((I2R + Psol*D*alphas)>(Prad+Pcon))
            GuessTcLeft = GuessTc;
        else
            GuessTcRight = GuessTc;
        end
        if(counter>=5000)
            error('error')
        end
    end
end