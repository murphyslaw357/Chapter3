function [GuessTc] =GetGuessTemp(I,Ta,D,phi,Vw,alpha,beta,epsilons,Psol,mdl)
    GuessTc=Ta;
    if(I>0)
        if(isempty(mdl))
            Pj=(I.^2).*(alpha+25*beta);
            a=0.43012;
            b=21.496;
            c=3.0337;
            d=-0.97255;
            e=1.087;
            GuessTc=a+b*Pj/(c+Vw)+d*Pj+e*Ta;
        else
            %GuessTc=coef(1)+coef(2).*Pj/(coef(3)+Vw)+coef(4)*Pj+coef(5)*Ta;
            %GuessTc=coef(1)+coef(2).*(Pj+coef(6).*Psol)./(coef(3)+Vw)+coef(4).*Pj+coef(5).*Ta;
            GuessTc=mdl(((I.^2)).*alpha,((I.^2)).*beta,Psol.*D,Ta,1./(Vw+1),(1./(Vw+1)).^2);
            %GuessTc=mdl((I^2)*(alpha+25*beta),Psol*D,Ta,Vw)+Ta;
        end
    end
% if(GuessTc<Ta+2)
%     GuessTc=Ta+2;
% end
% if(GuessTc>200)
%     GuessTc=200;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x=[winds,ambtemps,currents.*maxcurrent,psols.*diam.*alphas];
% pconfun=@(b,x) b(1).*x(:,1);
% pradfun=@(b,x) b(2)+b(3)*(x(:,3).^3)./(x(:,1)+1);
% modelfun = @(b,x)(alpha.*(x(:,3).^2) + x(:,4) +pconfun(b,x).*x(:,2)-pradfun(b,x).*(x(:,2)))./(pconfun(b,x)+pradfun(b,x)-(beta.*x(:,3).^2));
% beta0=[0.1 0.000639 3.02e-08 1 1 1];
% mdl = fitnlm(x,root+273,modelfun,beta0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IR2s=((currents.*maxcurrent).^2).*(alpha+25*beta);
% x=[IR2s,winds,ambtemps];
% modelfun=@(b,x) b(1)+b(2).*x(:,1)./(b(3)+x(:,2))+b(4).*x(:,1)+b(5).*x(:,3);
% beta0=[0.46 21.37 3.02 -0.9667 1.0864];
% mdl = fitnlm(x,root,modelfun,beta0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IR2s=((currents.*maxcurrent).^2).*(alpha+25*beta);
% x=[IR2s,psols*diam,ambtemps,winds];
% mdl=MultiPolyRegress(x,root,3);