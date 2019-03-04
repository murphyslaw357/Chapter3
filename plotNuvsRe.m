clear
clc
close all

grprlim1=0;
grprlim2=(0.889/0.675)^(1/(0.058-0.088));
grprlim3=(1.02/0.889)^(1/(0.088-0.148));
grprlim4=(0.85/1.02)^(1/(0.148-0.188));
grprlim5=(0.48/0.85)^(1/(0.188-0.25));
grprlim6=(0.125/0.48)^(1/(0.25-0.333));
grprlim7=1e12;
grprx1=(grprlim1:(grprlim2-grprlim1)/10:grprlim2);
grprx2=(grprlim2:(grprlim3-grprlim2)/100:grprlim3);
grprx3=(grprlim3:(grprlim4-grprlim3)/1000:grprlim4);
grprx4=(grprlim4:(grprlim5-grprlim4)/1000:grprlim5);
grprx5=(grprlim5:(grprlim6-grprlim5)/1000:grprlim6);
grprx6=(grprlim6:(grprlim7-grprlim6)/100000:grprlim7);

fgrpr1=0.675.*grprx1.^0.058;
fgrpr2=0.889.*grprx2.^0.088;
fgrpr3=1.02.*grprx3.^0.148;
fgrpr4=0.85.*grprx4.^0.188;
fgrpr5=0.48.*grprx5.^0.25;
fgrpr6=0.125.*grprx6.^0.333;

grprx=[grprx1(1:end-1) grprx2(1:end-1) grprx3(1:end-1) grprx4(1:end-1) grprx5(1:end-1) grprx6];
fgrpr=[fgrpr1(1:end-1) fgrpr2(1:end-1) fgrpr3(1:end-1) fgrpr4(1:end-1) fgrpr5(1:end-1) fgrpr6];

%[f,gof,out] = fit(grprx',fgrpr','smoothingspline','SmoothingParam',1);
[f,gof,out] = fit(grprx',fgrpr','a*x^b','startpoint',[1/4,1/3]);
% 
% % df = differentiate(f,grprx');
save('GrPrSpline.mat','f')
% 
% % fgrpr2=0.6+0.387.*((grprx./(1+((0.559)))).^(1/6))
clear
clc
% lim1=1e-4;
% lim2=4e-3;
% lim3=9e-2;
% lim4=1;
% lim5=35;
% lim6=5e3;
% lim7=5e4;
% lim8=2e5;

lim1=0;
lim2=(0.565/0.437)^(1/(0.0895-0.136));
lim3=(0.8/0.565)^(1/(0.136-0.280));
lim4=(0.795/0.8)^(1/(0.280-0.384));
lim5=(0.583/0.795)^(1/(0.384-0.471));
lim6=(0.148/0.583)^(1/(0.471-0.633));
lim7=(0.0208/0.148)^(1/(0.633-0.814));
lim8=2e5;

x1=(lim1:(lim2-lim1)/10:lim2);
x2=(lim2:(lim3-lim2)/100:lim3);
x3=(lim3:(lim4-lim3)/1000:lim4);
x4=(lim4:(lim5-lim4)/1000:lim5);
x5=(lim5:(lim6-lim5)/1000:lim6);
x6=(lim6:(lim7-lim6)/10000:lim7);
x7=(lim7:(lim8-lim7)/100000:lim8);

f1=0.437.*x1.^0.0895;
f2=0.565.*x2.^0.136;
f3=0.800.*x3.^0.280;
f4=0.795.*x4.^0.384;
f5=0.583.*x5.^0.471;
f6=0.148.*x6.^0.633;
f7=0.0208.*x7.^0.814;

x=[x1 x2 x3 x4 x5 x6 x7];
f=[f1 f2 f3 f4 f5 f6 f7];
% [ff,gof,out] = fit(x',f','smoothingspline','SmoothingParam',1)
[ff,gof,out] = fit(x',f','a*x^b','startpoint',[1/4,1/3]);
% [ffinv,gof,out] = fit(f',x','smoothingspline','SmoothingParam',1);
[ffinv,gof,out] = fit(f',x','a*x^b','startpoint',[1/4,1/3]);
save('ReNuSpline.mat','ff')
save('NuReSpline.mat','ffinv')
% loglog(x,f);
% 
% xlabel('Reynolds Number - Re')
% ylabel('Nusselt Number - Nu')
% legend('[1]','[3]')

