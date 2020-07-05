function [GuessTc] =GetGuessTemp(I,Ta,D,phi,Vw,alpha,beta,epsilons,alphas,Psol,mdl)
     GuessTc=Ta;
     sigmab=5.6697e-8;
    
%    I2Rapr=(alpha+Ta.*beta).*(I.^2);
%    GuessTc = Ta+(Psol.*D.*alphas+I2Rapr)./(krfit(Ta)+kcfit(Vw));
     %if(isempty(mdl))
     %    Pj=(I.^2).*(alpha+25.*beta);
     %    GuessTc(I>0)=(((Psol(I>0).*D+Pj(I>0))./(pi*D*sigmab*epsilons)+((Ta(I>0)+273).^4)).^(1/4))-273;
     %else
         GuessTc=Ta+mdl(((I.^2)).*alpha,((I.^2)).*beta,Psol.*D.*alphas,Ta,1./(Vw+1),(1./(Vw+1)).^2);
     %end
    
GuessTc(GuessTc<=Ta)=Ta(GuessTc<=Ta)+0.1;    
%if(GuessTc<Ta)
%    GuessTc=Ta+0.1;
%end