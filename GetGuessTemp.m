function [GuessTc,Pcon,Prad,Pj] =GetGuessTemp(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,Psol,fGrPr,fReNu,fNuRe)
%     sigmab=5.6697e-8;
%     g=9.805;
%     GuessTfilm=Ta+10;
%     lambdaf=(2.42e-2)+(7.2e-5)*GuessTfilm;    
%     Tfilmk=GuessTfilm+273;
%     vf=((1.32e-5)+(9.5e-8)*GuessTfilm)*((1-((6.5e-3)*H)/288.16)^-5.2561);
%     
%     Gr=(g*(D^3)*abs(10))/(Tfilmk*(vf^2));    
%     Pr=0.715-(2.5e-4).*GuessTfilm;
%     GrPr=Gr.*Pr;
%     Nudf=fGrPr(GrPr);
%     if(Vw==0)
%         %pure natural convection         
%         Pcon=pi*Nudf*lambdaf;
%     else
%     %mixed convection
%         Re=(sin(phi)*Vw*D)./vf;  
%         Req=fNuRe(Nudf);
%         Reeff=sqrt((Re.^2)+(Req.^2));
%         Nueff=fReNu(Reeff);
%         Pcon=pi*Nueff*lambdaf;
%     end
% 
%     %Ref=sin(phi)*Vw*D/vf;
%     %Reeff=sqrt(0.2971+(Ref)^2);
% 
%     %kf=2.424e-2+(7.477e-5)*GuessTfilm-(4.407e-9)*(GuessTfilm^2);
%     %muf=((1.458e-6)*(GuessTfilm+273)^1.5)/(GuessTfilm+383.4);
%     %rhof=(1.293-(1.525e-4)*H+(6.379e-9)*(H^2))/(1+0.00367*GuessTfilm);
% %     if(Vw<0.2)
% %         Vw=0.2;
% %     end
%     %Nre=D*rhof*Vw/muf;
%     %Pcon=(1.01+1.35*Nre^0.52)*kf;
% 
%     Prad=pi*D*sigmab*epsilons*((1.38e8)+Ta*(1.39e6));
%     GuessTc=(Psol+(I^2)*(alpha+beta*25))/(Prad+Pcon)+Ta;
%     GuessTc2=(Psol+(I^2)*alpha+Ta*(Prad+Pcon))/(Prad+Pcon-(I^2)*beta);
%     if(abs(GuessTc2-GuessTc)>25)
%     end
Pcon=0;
Prad=0;
Pj=((I^2)*(alpha+beta*25));
GuessTc=Ta;
if(I>0)
%      a=0.57711;
%      b=1.9637;
%      c=0.56494;
%      d=0.29798;
%      e=1.0489;
a=2.4113;
b=1.8582;
c=0.53887;
d=0.27644;
e=1.0476;
     GuessTc=a+b*Pj/(c+Vw)+d*Pj+e*Ta;
     %x=[((currents.*maxcurrent).^2).*(alpha+25*beta),winds,ambtemps];
end
if(GuessTc<Ta)
    GuessTc=Ta;
end