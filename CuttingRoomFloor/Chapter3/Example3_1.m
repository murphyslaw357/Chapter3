I=500/2;
Ta=21;
H=0;
D=0.029591;
phi=pi/2;
Vw=4;
alpha=5.3686e-5;
beta=2.6843e-7;
epsilons=0.7;
Psol=29.591;

[GuessTc,I2R,I2Rprime,Prad,PradPrime,Pcon,PconPrime] =GetTempNewton(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,Psol)