function [GuessTc] =GetGuessTemp(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,Psol,fGrPr,fReNu,fNuRe)
% Pcon=0;
% Prad=0;
Pj=((I^2)*(alpha+beta*25));
GuessTc=Ta;

if(I>0)
    a=0.43012;
    b=21.496;
    c=3.0337;
    d=-0.97255;
    e=1.087;
    GuessTc=a+b*Pj/(c+Vw)+d*Pj+e*Ta;
end
% if(GuessTc<Ta+2)
%     GuessTc=Ta+2;
% end
% IR2s=((currents.*maxcurrent).^2).*(alpha+25*beta);
%x=[IR2s,winds,ambtemps];
%modelfun = @(b,x)b(1) + b(2)*x(:,1)./(b(3)+x(:,2)) + b(4)*x(:,1)+b(5).*x(:,3);
%beta0 = [2.4 1.85 0.5 0.27 1.04];
%mdl = fitnlm(x,root,modelfun,beta0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x=[winds,ambtemps,currents];
% modelfun = @(b,x)(alpha.*(x(:,3).^2) + psols -b(1).*x(:,1)-b(2).*x(:,2).*x(:,1)-b(3))./(b(4).*x(:,1)+b(5)+b(6).*(x(:,3).^2)-(beta.*x(:,3).^2));
% beta0=[0.58 -0.1312 1 0.123 1 37];
% mdl = fitnlm(x,root,modelfun,beta0)