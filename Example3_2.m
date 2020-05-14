clear
clc
close all

if(ispc==1)
    foldersource='C:\Users\ctc\Documents\GitHub\Chapter3\';
elseif(ismac==1)
    foldersource='/Users/Shaun/Documents/GitHub/Chapter3/';
elseif(isunix==1)
    foldersource='/mnt/HA/groups/nieburGrp/Shaun/Chapter3/';
end

k=8.85e-12;
D25=5;
Deq=10;
r=0.0296/2;
S=200;
Tref=25;
Lref=8*(D25^2)/(3*S)+S;
alphaas=2.2e-5;

Temp=-30:0.01:100;
Trefloc=find(Temp==Tref);
[~,col]=size(Temp);
L=(Lref.*(1+alphaas.*(Temp-Tref)));
D=(20-sqrt((3*S*(L-S))/8));
WoverH=((20-D).*2./((S/2)^2));
Havg=(2/S).*(WoverH.*((S/2)^3)/3+D.*(S/2));

H1=2.*Havg;
H2=H1;
H3=H1+2*sqrt(Deq^2-(Deq/2)^2);
H12=sqrt(Deq^2+H1.^2);
H23=sqrt((Deq/2)^2+(H1+sqrt(Deq^2-(Deq/2)^2)).^2);
H31=H23;
C1=(2*pi*k)./(log(Deq/r)).*ones(1,col).*S;
C2=(2*pi*k)./(log(Deq/r)).*ones(1,col).*L;
C3=((2*pi*k)./(log(Deq/r)-log(((H12(Trefloc).*H23(Trefloc).*H31(Trefloc)).^(1/3))./((H1(Trefloc).*H2(Trefloc).*H3(Trefloc)).^(1/3))))).*S;
C4=((2*pi*k)./(log(Deq/r)-log(((H12(Trefloc).*H23(Trefloc).*H31(Trefloc)).^(1/3))./((H1(Trefloc).*H2(Trefloc).*H3(Trefloc)).^(1/3))))).*L;
C5=((2*pi*k)./(log(Deq/r)-log(((H12.*H23.*H31).^(1/3))./((H1.*H2.*H3).^(1/3))))).*S;
C6=((2*pi*k)./(log(Deq/r)-log(((H12.*H23.*H31).^(1/3))./((H1.*H2.*H3).^(1/3))))).*L;

figure('Renderer', 'painters', 'Position', [10 10 500 600]);
plot(Temp,100.*abs(C6-C1)./C6,'linewidth',1)
hold on
plot(Temp,100.*abs(C6-C2)./C6,'linewidth',1)
plot(Temp,100.*abs(C6-C3)./C6,'linewidth',1)
plot(Temp,100.*abs(C6-C4)./C6,'linewidth',1)
plot(Temp,100.*abs(C6-C5)./C6,'linewidth',1)
xlabel('Conductor Temperature - ^{\circ}C')
ylabel('Absolute Error - %')
legend('Fixed Length, Neglected Ground','Variable Length, Neglected Ground','Fixed Length, Fixed Ground','Variable Length, Fixed Ground','Fixed Length, Ground Included','Location', 'southoutside')
xlim([-40 110])

set(gca,'FontSize',10)
savefig(strcat(foldersource,'Figure3_5.fig'))
set(gcf, 'Color', 'w');
if(ispc)
    export_fig C:\Users\ctc\Documents\GitHub\Chapter3\Figure3_5.png -m3
elseif(ismac)
    export_fig /Volumes/THESIS/Github/Chapter3/Figure3_5.png -m3
end