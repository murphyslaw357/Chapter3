clear
clc
close all

grprlim1=0;
grprlim2=1e-04;
grprlim3=0.1;
grprlim4=95;
grprlim5=1e4;
grprlim6=1e7;
grprlim7=1e12;
grprx1=(grprlim1:(grprlim2-grprlim1)/100:grprlim2);
grprx2=(grprlim2:(grprlim3-grprlim2)/100:grprlim3);
grprx3=(grprlim3:(grprlim4-grprlim3)/100:grprlim4);
grprx4=(grprlim4:(grprlim5-grprlim4)/100:grprlim5);
grprx5=(grprlim5:(grprlim6-grprlim5)/100:grprlim6);
grprx6=(grprlim6:(grprlim7-grprlim6)/100:grprlim7);

fgrpr1=0.675.*grprx1.^0.058;
fgrpr2=0.889.*grprx2.^0.088;
fgrpr3=1.02.*grprx3.^0.148;
fgrpr4=0.85.*grprx4.^0.188;
fgrpr5=0.48.*grprx5.^0.25;
fgrpr6=0.125.*grprx6.^0.333;

grprx=[grprx1 grprx2 grprx3 grprx4 grprx5 grprx6];
fgrpr=[fgrpr1 fgrpr2 fgrpr3 fgrpr4 fgrpr5 fgrpr6];


f = fit(grprx',fgrpr','smoothingspline');

f2= spline(grprx',fgrpr');
% df = differentiate(f,grprx');
save('GrPrSpline.mat')

fgrpr2=0.6+0.387.*((grprx./(1+((0.559)))).^(1/6))
clear
clc
lim1=1e-4;
lim2=4e-3;
lim3=9e-2;
lim4=1;
lim5=35;
lim6=5e3;
lim7=5e4;
lim8=2e5;
x1=(lim1:(lim2-lim1)/100:lim2);
x2=(lim2:(lim3-lim2)/100:lim3);
x3=(lim3:(lim4-lim3)/100:lim4);
x4=(lim4:(lim5-lim4)/100:lim5);
x5=(lim5:(lim6-lim5)/100:lim6);
x6=(lim6:(lim7-lim6)/100:lim7);
x7=(lim7:(lim8-lim7)/100:lim8);

f1=0.437.*x1.^0.0895;
f2=0.565.*x2.^0.136;
f3=0.800.*x3.^0.280;
f4=0.795.*x4.^0.384;
f5=0.583.*x5.^0.471;
f6=0.148.*x6.^0.633;
f7=0.0208.*x7.^0.814;

x=[x1 x2 x3 x4 x5 x6 x7];
f=[f1 f2 f3 f4 f5 f6 f7];
ff = fit(x',f','smoothingspline');
ffinv = fit(f',x','smoothingspline');

% %IEEE 738
% g1=1.01+1.35.*x.^0.52;
% g2=0.754.*x.^0.6;
% [~,col]=size(g1);
% g=zeros(1,col);
% for i=1:col
%     g(i)=max(g1(i),g2(i));
% end
% loglog(x,g)
% hold on
loglog(x,f);

xlabel('Reynolds Number - Re')
ylabel('Nusselt Number - Nu')
legend('[1]','[3]')

