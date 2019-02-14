function [GuessTc,Pcon,Prad,Pj,Nueff,lambdaf,Nre] =GetGuessTemp(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,Psol)
    Relim1=(0.437/0.565)^(1/(0.136-0.0895));
    Relim2=(0.565/0.8)^(1/(0.280-0.136));
    Relim3=(0.8/0.795)^(1/(0.384-0.280));
    Relim4=(0.795/0.583)^(1/(0.471-0.384));
    Relim5=(0.583/0.148)^(1/(0.633-0.471));
    Relim6=(0.148/0.0208)^(1/(0.814-0.633));
    sigmab=5.6697e-8;
    GuessTfilm=Ta;
    lambdaf=2.42e-2+(7.2e-5)*GuessTfilm;
    vf=((1.32e-5)+(9.5e-8)*GuessTfilm)*((1-((6.5e-3)*H)/288.16)^-5.2561);
    Ref=sin(phi)*Vw*D/vf;
    Reeff=sqrt(0.2971+(Ref)^2);
    C=0;
    n=0;
    if(Reeff<=Relim1)
        C=0.437;
        n=0.0895;
    elseif(Reeff>Relim1 && Reeff<=Relim2)
        C=0.565;
        n=0.136;
    elseif(Reeff>Relim2 && Reeff<=Relim3)
        C=0.8;
        n=0.28;
    elseif(Reeff>Relim3 && Reeff<=Relim4)
        C=0.795;
        n=0.384;
    elseif(Reeff>Relim4 && Reeff<=Relim5)
        C=0.583;
        n=0.471;
    elseif(Reeff>Relim5 && Reeff<=Relim6)
        C=0.148;
        n=0.633;
    elseif(Reeff>Relim6)
        C=0.0208;
        n=0.814;
    end
    
    Nueff=C*Reeff^n;
    Prad=pi*D*sigmab*epsilons*((1.38e8)+Ta*(1.39e6));
    kf=2.424e-2+(7.477e-5)*GuessTfilm-(4.407e-9)*(GuessTfilm^2);
    muf=((1.458e-6)*(GuessTfilm+273)^1.5)/(GuessTfilm+383.4);
    rhof=(1.293-(1.525e-4)*H+(6.379e-9)*(H^2))/(1+0.00367*GuessTfilm);
    if(Vw<0.2)
        Vw=0.2;
    end
    Nre=D*rhof*Vw/muf;
    Pcon=(1.01+1.35*Nre^0.52)*kf;
%     Pcon=pi*Nueff*lambdaf;
    GuessTc=(Psol+(I^2)*(alpha+beta*25))/(Prad+pi*Nueff*lambdaf)+Ta;
%     GuessTc=(((I^2)*(alpha))/(Prad+Pcon)+Psol/(Prad+Pcon)+Ta)/(1+((I^2)*beta)/(Prad+Pcon));

    
%      GuessTc=(Psol+(I^2)*alpha+Ta*(Prad+pi*Nueff*lambdaf))/(Prad+pi*Nueff*lambdaf-(I^2)*beta);
%     GuessTc2=(Psol+(I^2)*alpha)/(Prad+pi*Nueff*lambdaf-(I^2)*beta)+Ta*((Prad+pi*Nueff*lambdaf))/(Prad+pi*Nueff*lambdaf-(I^2)*beta);
%     if(round(GuessTc-GuessTc)~=0)
%     end
%     if(GuessTc2<Ta)
%     end
    Pj=((I^2)*(alpha+beta*GuessTc));
end