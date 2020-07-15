function [GuessTcOutput,I2R,I2Rprime,Prad,PradPrime,PradPrimePrime,Pcon,...
    PconPrime,PconPrimePrime,Gr,GrPrime,Nudf] = ...
    GetTempNewtonFirstIteration2(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,alphas,...
    Psol,GuessTc)
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
    gplim1=1e-10;
    gplim2=1e-4;
    gplim3=1e-1;
    gplim4=1e2;
    gplim5=1e4;
    gplim6=1e7;
    gplim7=1e12;
    
    IIstar=abs(I).^2;
    %%                          RADIATIVE COOLING/HEATING                %%
    Prad=(pi*D*sigmab*epsilons).*(((GuessTc+273).^4)-((Ta+273).^4));
    PradPrime=(4*pi*D*sigmab*epsilons).*(GuessTc+273).^3;
    PradPrimePrime=(12*pi*D*sigmab*epsilons).*(GuessTc+273).^2;
    %%                          JOULE HEATING                            %%
    I2Rprime=beta.*IIstar;
    Resistance=beta.*GuessTc+alpha;
    if(Resistance<0)
        Resistance=0;
    end
    I2R=IIstar*Resistance;
    %%                      CONVECTIVE COOLING/HEATING                   %%
    vfPrimePrime=0;
    PrPrime=-1.25e-4;
    GuessTfilm=(GuessTc+Ta)./2;
    vf=((1.32e-5)+(9.5e-8).*GuessTfilm).*((1-((6.5e-3)*H)/288.16)^-5.2561);
    vfPrime=(4.75e-8)*((1-((6.5e-3)*H)/288.16)^-5.2561);
    Tfilmk=GuessTfilm+273;
    TfilmkPrime=1/2;
    Lambdaf=(2.42e-2)+(7.2e-5)*GuessTfilm;
    LambdafPrime=3.6e-5;
    
    gD3 = g*(D^3);
    Gr = gD3.*(GuessTc-Ta)./(Tfilmk.*(vf.^2));
    GrPrime = gD3.*(Ta.*vf-((GuessTc.^2)-Ta.^2).*vfPrime)./...
        ((Tfilmk.^2).*(vf.^3));
    GrPrimePrime=(-2.*gD3.*GuessTc./((Tfilmk.^2).*(vf.^3))).*vfPrime-...
        (gD3./((Tfilmk.^3).*(vf.^4))).*(vf+3.*Tfilmk.*vfPrime).*...
        (Ta.*vf-((GuessTc.^2)-Ta.^2).*vfPrime);
    Pr=0.715-(2.5e-4).*GuessTfilm;
    GrPr=Gr*Pr;
    
    %Natural convection
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
    Nudf=A*GrPr^m;
    NudfPrimedgrpr=m*A*GrPr^(m-1);
    NudfPrimePrimedgrpr2=(m-1)*m*A*GrPr^(m-2);
%     Nudf=fGrPr(GrPr);
%     NudfPrimedgrpr=fGrPr.a.*fGrPr.b.*(GrPr.^(fGrPr.b-1));
%     NudfPrimePrimedgrpr2=(fGrPr.b-1).*fGrPr.a.*fGrPr.b.*(GrPr.^(fGrPr.b-2));
    
    NudfPrimedtc=NudfPrimedgrpr.*(Gr.*PrPrime+GrPrime.*Pr);
    NudfPrimePrimedtc2=NudfPrimePrimedgrpr2.*(Gr.*PrPrime+GrPrime.*Pr)+...
        NudfPrimedgrpr.*...
        (GrPrime.*PrPrime+GrPrimePrime.*Pr+GrPrime.*PrPrime);

    %Mixed convection
    Re=(sin(phi)*Vw*D)./vf;  
    RePrimedtc=-((sin(phi)*Vw*D).*vfPrime)./(vf.^2);
    RePrimePrimedtc2=(sin(phi)*Vw*D).*...
        (2.*(vfPrime.^2)-1.*vf.*vfPrimePrime)./(vf.^3);
    if(Nudf>0.1916 && Nudf<=0.2666)
        Cinv=0.437;
        ninv=0.0895;
    elseif(Nudf>0.2666 && Nudf<=0.40685)
       Cinv=0.565;
       ninv=0.136;
    elseif(Nudf>0.40685 && Nudf<=0.8136)
       Cinv=0.8;
       ninv=0.28;
    elseif(Nudf>0.8136 && Nudf<=3.1253)
       Cinv=0.795;
       ninv=0.384;
    elseif(Nudf>3.1253 && Nudf<=31.45)
       Cinv=0.583;
       ninv=0.471;
    elseif(Nudf>31.45 && Nudf<=141.449)
       Cinv=0.148;
       ninv=0.633;
    elseif(Nudf>141.449 && Nudf<=429.6371)
       Cinv=0.0208;
       ninv=0.814;
    else
        error('out of bounds')
    end  
    Req=(Nudf/Cinv)^(1/ninv);
    ReqPrimednudf=(1/ninv)*(Nudf/Cinv)^(1/ninv-1);    
    ReqPrimePrimednudf=(1/ninv-1)*(1/ninv)*(Nudf/Cinv)^(1/ninv-2); 
    
    %Req=fNuRe(Nudf);
    %[ReqPrimednudf,ReqPrimePrimednudf]=differentiate(fNuRe,Nudf);
    ReqPrimedtc=ReqPrimednudf.*NudfPrimedtc;
    ReqPrimePrimedtc2=ReqPrimePrimednudf.*NudfPrimedtc+ReqPrimednudf.*...
        NudfPrimePrimedtc2;
    Reeff=sqrt((Re.^2)+(Req.^2));
    Top=Re.*RePrimedtc+Req.*ReqPrimedtc;
    dtop=RePrimedtc.*RePrimedtc+Re.*RePrimePrimedtc2+ReqPrimedtc.*...
        ReqPrimedtc+Req.*ReqPrimePrimedtc2;
    
    ReeffPrimedtc=Top./Reeff;
    ReeffPrimePrimedtc2=(Reeff.*dtop-Top.*ReeffPrimedtc)./(Reeff.^2);
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
    NueffPrimedreeff=C*n*(Reeff^(n-1));
    NueffPrimePrimedreeff2=(n-1)*C*n*(Reeff^(n-2));
    %Nueff=fReNu(Reeff);
    %[NueffPrimedreeff,NueffPrimePrimedreeff2]=differentiate(fReNu,Reeff);
    NueffPrimedtc=NueffPrimedreeff.*ReeffPrimedtc;
    NueffPrimePrimedtc2=NueffPrimePrimedreeff2.*ReeffPrimedtc+...
        NueffPrimedreeff.*ReeffPrimePrimedtc2;
%     %Forced convection
%     Nu=fReNu(Re);  
%     [NuPrimedre,NuPrimePrimedre2]=differentiate(fReNu,Re);
%     NuPrimedtc=NuPrimedre.*RePrimedtc;
%     NuPrimePrimedtc2=NuPrimePrimedre2.*RePrimedtc+...
%         NuPrimedre.*RePrimePrimedtc2;
    
    [row,~]=size(GuessTc);
    GuessTcOutput=GuessTc;
    Pcon=zeros(row,1);
    PconPrime=zeros(row,1);
    PconPrimePrime=zeros(row,1);
    for i=1:row
%         if(Vw==0)
%             %pure natural convection         
%             Pcon(i)=pi*Nudf(i)*Lambdaf(i)*(GuessTc(i)-Ta);
%             PconPrime(i)=pi*(Nudf(i)*Lambdaf(i)+(GuessTc(i)-Ta)*(LambdafPrime*Nudf(i)+NudfPrimedtc(i)*Lambdaf(i)));
%             PconPrimePrime(i)=pi*(Nudf(i)*LambdafPrime+NudfPrimedtc(i)*Lambdaf(i)+(LambdafPrime*Nudf(i)+NudfPrimedtc(i)*Lambdaf(i))+(GuessTc(i)-Ta)*(LambdafPrime*NudfPrimedtc(i)+NudfPrimePrimedtc2(i)*Lambdaf(i)+NudfPrimedtc(i)*LambdafPrime));
%         elseif(round(GuessTc(i),4)~=round(Ta,4))
        %mixed convection
            Pcon(i)=pi*Nueff(i)*Lambdaf(i)*(GuessTc(i)-Ta);
            PconPrime(i)=pi*(Nueff(i)*Lambdaf(i)+(GuessTc(i)-Ta)*(LambdafPrime*Nueff(i)+NueffPrimedtc(i)*Lambdaf(i)));
            PconPrimePrime(i)=pi*(Nueff(i)*LambdafPrime+NueffPrimedtc(i)*Lambdaf(i)+LambdafPrime*Nueff(i)+NueffPrimedtc(i)*Lambdaf(i)+(GuessTc(i)-Ta)*(LambdafPrime*NueffPrimedtc(i)+NueffPrimePrimedtc2(i)*Lambdaf(i)+NueffPrimedtc(i)*LambdafPrime));
%         else
%         %pure forced    
%             Pcon(i)=pi*Nu(i)*Lambdaf(i)*(GuessTc(i)-Ta);
%             PconPrime(i)=pi*(Nu(i)*Lambdaf(i)+(GuessTc(i)-Ta)*(LambdafPrime*Nu(i)+NuPrimedtc(i)*Lambdaf(i)));
%             PconPrimePrime(i)=pi*(Nu(i)*LambdafPrime+NuPrimedtc(i)*Lambdaf(i)+LambdafPrime*Nu(i)+NuPrimedtc(i)*Lambdaf(i)+(GuessTc(i)-Ta)*(LambdafPrime*NuPrimedtc(i)+NuPrimePrimedtc2(i)*Lambdaf(i)+NuPrimedtc(i)*LambdafPrime));
%         end
        
    %%                        MISMATCH AND UPDATE                        %%
    Hprime=I2Rprime-PradPrime(i)-PconPrime(i);
    if(Hprime>0 || isnan(PconPrime(i)))
        msg='error condition1';
        disp(msg);
    end
    
    Mismatch=I2R(i)+Psol*D*alphas-Prad(i)-Pcon(i);
    update=Mismatch/Hprime;
    GuessTcOutput(i)=GuessTcOutput(i)-update;  
    end
end