function [GuessTc] =GetGuessTemp(I,Ta,D,phi,Vw,alpha,beta,epsilons,alphas,Psol,mdl)

    GuessTc=Ta+mdl(((I.^2)).*alpha,((I.^2)).*beta,Psol.*D.*alphas,Ta,1./(Vw+1),(1./(Vw+1)).^2);
    GuessTc(GuessTc<=Ta)=Ta(GuessTc<=Ta)+0.001;