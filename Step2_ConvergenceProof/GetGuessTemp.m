function [GuessTc] =GetGuessTemp(I,Ta,D,phi,Vw,alpha,beta,epsilons,Psol,mdl)
    GuessTc=Ta;
    sigmab=5.6697e-8;
    
    
    if(isempty(mdl))
        Pj=(I.^2).*(alpha+25.*beta);
        GuessTc(I>0)=(((Psol(I>0).*D+Pj(I>0))./(pi*D*sigmab*epsilons)+((Ta(I>0)+273).^4)).^(1/4))-273;
    else
        GuessTc(I>0)=Ta(I>0)+mdl(((I(I>0).^2)).*alpha,((I(I>0).^2)).*beta,Psol(I>0).*D,Ta(I>0),1./(Vw(I>0)+1),(1./(Vw(I>0)+1)).^2);
    end
    
GuessTc(GuessTc<=Ta)=Ta(GuessTc<=Ta)+0.1;    
%if(GuessTc<Ta)
%    GuessTc=Ta+0.1;
%end