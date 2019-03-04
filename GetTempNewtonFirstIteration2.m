function [GuessTc,I2R,I2Rprime,Prad,PradPrime,PradPrimePrime,Pcon,PconPrime,PconPrimePrime,Gr,GrPrime,Nudf] =GetTempNewtonFirstIteration2(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,Psol,GuessTc,fGrPr,fReNu,fNuRe)
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
    IIstar=abs(I)^2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RADIATIVE COOLING/HEATING%%%%%%%%%%%
    Prad=(pi*D*sigmab*epsilons).*(((GuessTc+273).^4)-((Ta+273).^4));
    PradPrime=(4*pi*D*sigmab*epsilons).*(GuessTc+273).^3;
    PradPrimePrime=(12*pi*D*sigmab*epsilons).*(GuessTc+273).^2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%JOULE HEATING%%%%%%%%%%%%%%%%%%%%%%%
    I2Rprime=beta*IIstar;
    Resistance=beta*GuessTc+alpha;
    if(Resistance<0)
        Resistance=0;
    end
    I2R=IIstar*Resistance;
    %%%%%%%%%%%%%%%%%%%%%%%%CONVECTIVE COOLING/HEATING%%%%%%%%%%%%%%%%%%%%%
    vfPrimePrime=0;
    PrPrime=-1.25e-4;
    GuessTfilm=(GuessTc+Ta)./2;
    vf=((1.32e-5)+(9.5e-8).*GuessTfilm).*((1-((6.5e-3)*H)/288.16)^-5.2561);
    vfPrime=(4.75e-8)*((1-((6.5e-3)*H)/288.16)^-5.2561);
    Tfilmk=GuessTfilm+273;
    TfilmkPrime=1/2;
    Lambdaf=(2.42e-2)+(7.2e-5)*GuessTfilm;
    LambdafPrime=3.6e-5;
    
    [row,~]=size(GuessTc);
%     Gr=zeros(row,1);
    Gr=g*(D^3).*abs(GuessTc-Ta)./(Tfilmk.*(vf.^2));
    GrPrime=g*(D^3).*((Tfilmk.*(vf.^2)).*((GuessTc-Ta)./abs(GuessTc-Ta))-abs(GuessTc-Ta).*(TfilmkPrime.*vf.^2+...
                2.*Tfilmk.*vf.*vfPrime))./((Tfilmk.^2).*(vf.^4));
    topp=(Tfilmk.*(vf.^2)).*((GuessTc-Ta)./abs(GuessTc-Ta))-abs(GuessTc-Ta).*(TfilmkPrime.*vf.^2+2.*Tfilmk.*vf.*vfPrime);
    t1=Tfilmk.*vf.^2;
    t1prime=TfilmkPrime.*(vf.^2)+2.*Tfilmk.*vf.*vfPrime;
    t2=(GuessTc-Ta)./abs(GuessTc-Ta);
    t2prime=(abs(GuessTc-Ta)-(GuessTc-Ta).*(GuessTc-Ta)./abs(GuessTc-Ta))./(abs(GuessTc-Ta).^2);
    t3=abs(GuessTc-Ta);
    t3prime=(GuessTc-Ta)./abs(GuessTc-Ta);
    t4=TfilmkPrime.*vf.^2+2.*Tfilmk.*vf.*vfPrime;
    t4prime=vf.*vfPrime+2.*TfilmkPrime.*vf.*vfPrime+2.*Tfilmk.*vfPrime.*vfPrime;
    dtopp=t1prime.*t2+t1.*t2prime-t3prime.*t4-t3.*t4prime;            
    GrPrimePrime=g*(D^3).*(((Tfilmk.^2).*(vf.^4)).*dtopp-topp.*(Tfilmk.*(vf.^4)+(Tfilmk.^2).*4.*(vf.^3).*vfPrime))./((Tfilmk.^4).*(vf.^8));
%     GrPrime=zeros(row,1);
%     GrPrimePrime=zeros(row,1);
%     for i=1:row
%         if(GuessTc(i)>Ta)
%              
% %             GrPrime(i)=(g*(D^3)).*(Tfilmk(i)*vf(i)-(GuessTc(i)-Ta)*(TfilmkPrime*vf(i)+...
% %                 Tfilmk(i)*2.*vfPrime))/((Tfilmk(i)^2).*(vf(i)^3));
%             dtop=TfilmkPrime*vf(i)+Tfilmk(i)*vfPrime-(GuessTc(i)-Ta)*...
%                 (TfilmkPrime.*vfPrime+vfPrime+2*Tfilmk(i)*vfPrimePrime)-...
%                 (TfilmkPrime.*vf(i)+2.*Tfilmk(i).*vfPrime);
%             GrPrimePrime(i)=(g*(D^3)).*(((Tfilmk(i)^2).*(vf(i)^3)).*dtop-...
%                 ((Tfilmk(i)*vf(i)-(GuessTc(i)-Ta)*(TfilmkPrime*vf(i)+2*Tfilmk(i)*vfPrime)))*...
%                 (Tfilmk(i)*(vf(i)^3)+3*(vf(i)^2)*vfPrime*(Tfilmk(i)^2)))/...
%                 ((Tfilmk(i)^4).*(vf(i)^6));
%         elseif(GuessTc(i)<Ta)
%             
% %             GrPrime(i)=(-1*g*(D^3))*(Tfilmk(i)*vf(i)-(GuessTc(i)-Ta)*(TfilmkPrime*vf(i)+...
% %                 Tfilmk(i)*2*vfPrime))/((Tfilmk(i)^2)*(vf(i)^3));
%             dtop=TfilmkPrime*vf(i)+Tfilmk(i)*vfPrime-(GuessTc(i)-Ta)*...
%                 (TfilmkPrime*vfPrime+2*TfilmkPrime*vfPrime+2*Tfilmk(i)*vfPrimePrime)-...
%                 (TfilmkPrime*vf(i)+2*Tfilmk(i)*vfPrime);
%             GrPrimePrime(i)=(-1*g*(D^3)).*(((Tfilmk(i)^2).*(vf(i)^3)).*dtop-...
%                 ((Tfilmk(i)*vf(i)-(GuessTc(i)-Ta)*(TfilmkPrime*vf(i)+2*Tfilmk(i)*vfPrime)))*...
%                 (Tfilmk(i)*(vf(i)^3)+3*(vf(i)^2)*vfPrime*(Tfilmk(i)^2)))/...
%                 ((Tfilmk(i)^4).*(vf(i)^6));
%         end
%     end
    Pr=0.715-(2.5e-4).*GuessTfilm;
    GrPr=Gr.*Pr;
    
    %Natural convection
    Nudf=fGrPr(GrPr);
    %%
    [NudfPrimedgrpr,NudfPrimePrimedgrpr2]=differentiate(fGrPr,GrPr); 

    NudfPrimedtc=NudfPrimedgrpr.*(Gr.*PrPrime+GrPrime.*Pr);
    NudfPrimePrimedtc2=NudfPrimePrimedgrpr2.*(Gr.*PrPrime+GrPrime.*Pr)+...
        NudfPrimedgrpr.*(GrPrime.*PrPrime+GrPrimePrime.*Pr+GrPrime.*PrPrime);

    %Mixed convection
    Re=(sin(phi)*Vw*D)./vf;  
    RePrimedtc=-((sin(phi)*Vw*D).*vfPrime)./(vf.^2);
    RePrimePrimedtc2=(sin(phi)*Vw*D).*...
        (2.*(vfPrime.^2)-1.*vf.*vfPrimePrime)./(vf.^3);
    Req=fNuRe(Nudf);
    [ReqPrimednudf,ReqPrimePrimednudf]=differentiate(fNuRe,Nudf);
    ReqPrimedtc=ReqPrimednudf.*NudfPrimedtc;
    ReqPrimePrimedtc2=ReqPrimePrimednudf.*NudfPrimedtc+ReqPrimednudf.*NudfPrimePrimedtc2;
    Reeff=sqrt((Re.^2)+(Req.^2));
    Top=Re.*RePrimedtc+Req.*ReqPrimedtc;
    dtop=RePrimedtc.*RePrimedtc+Re.*RePrimePrimedtc2+ReqPrimedtc.*ReqPrimedtc+Req.*ReqPrimePrimedtc2;
    
    ReeffPrimedtc=Top./Reeff;
    ReeffPrimePrimedtc2=(Reeff.*dtop-Top.*ReeffPrimedtc)./(Reeff.^2);
    Nueff=fReNu(Reeff);
    [NueffPrimedreeff,NueffPrimePrimedreeff2]=differentiate(fReNu,Reeff);
    NueffPrimedtc=NueffPrimedreeff.*ReeffPrimedtc;
    NueffPrimePrimedtc2=NueffPrimePrimedreeff2.*ReeffPrimedtc+...
        NueffPrimedreeff.*ReeffPrimePrimedtc2;
    %Forced convection
    Nu=fReNu(Re);  
    [NuPrimedre,NuPrimePrimedre2]=differentiate(fReNu,Re);
    NuPrimedtc=NuPrimedre.*RePrimedtc;
    NuPrimePrimedtc2=NuPrimePrimedre2.*RePrimedtc+...
        NuPrimedre.*RePrimePrimedtc2;
    
    Pcon=zeros(row,1);
    PconPrime=zeros(row,1);
    PconPrimePrime=zeros(row,1);
    for i=1:row
%         if(GuessTc(i)~=Ta && Vw==0)
%             %pure natural convection         
%             Pcon(i)=pi*Nudf(i)*Lambdaf(i)*(GuessTc(i)-Ta);
%             PconPrime(i)=pi*(Nudf(i)*Lambdaf(i)+(GuessTc(i)-Ta)*(LambdafPrime*Nudf(i)+NudfPrimedtc(i)*Lambdaf(i)));
%             PconPrimePrime(i)=pi*(Nudf(i)*LambdafPrime+NudfPrimedtc(i)*Lambdaf(i)+(LambdafPrime*Nudf(i)+NudfPrimedtc(i)*Lambdaf(i))+(GuessTc(i)-Ta)*(LambdafPrime*NudfPrimedtc(i)+NudfPrimePrimedtc2(i)*Lambdaf(i)+NudfPrimedtc(i)*LambdafPrime));
%         else
            if((GuessTc(i)-Ta)>0.1)% && Vw ~=0)
        %mixed convection
            Pcon(i)=pi*Nueff(i)*Lambdaf(i)*(GuessTc(i)-Ta);
            PconPrime(i)=pi*(Nueff(i)*Lambdaf(i)+(GuessTc(i)-Ta)*(LambdafPrime*Nueff(i)+NueffPrimedtc(i)*Lambdaf(i)));
            PconPrimePrime(i)=pi*(Nueff(i)*LambdafPrime+NueffPrimedtc(i)*Lambdaf(i)+LambdafPrime*Nueff(i)+NueffPrimedtc(i)*Lambdaf(i)+(GuessTc(i)-Ta)*(LambdafPrime*NueffPrimedtc(i)+NueffPrimePrimedtc2(i)*Lambdaf(i)+NueffPrimedtc(i)*LambdafPrime));
        elseif(GuessTc(i)-Ta<=0.1 && Vw~=0)
        %pure forced    
            Pcon(i)=pi*Nu(i)*Lambdaf(i)*(GuessTc(i)-Ta);
            PconPrime(i)=pi*(Nu(i)*Lambdaf(i)+(GuessTc(i)-Ta)*(LambdafPrime*Nu(i)+NuPrimedtc(i)*Lambdaf(i)));
            PconPrimePrime(i)=pi*(Nu(i)*LambdafPrime+NuPrimedtc(i)*Lambdaf(i)+LambdafPrime*Nu(i)+NuPrimedtc(i)*Lambdaf(i)+(GuessTc(i)-Ta)*(LambdafPrime*NuPrimedtc(i)+NuPrimePrimedtc2(i)*Lambdaf(i)+NuPrimedtc(i)*LambdafPrime));
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MISMATCH AND UPDATE%%%%%%%%%%%%%%
    Hprime=I2Rprime-PradPrime(i)-PconPrime(i);
    if(Hprime>0 || isnan(PconPrime(i)))
        msg='error condition1';
        disp(msg);
    end
    
    Mismatch=I2R(i)+Psol*D-Prad(i)-Pcon(i);
    update=Mismatch/Hprime;
    GuessTc(i)=GuessTc(i)-update;  
    end
end