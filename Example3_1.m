clear
clc
close all

if(ispc==1)
    foldersource='C:\Users\ctc\Documents\GitHub\Chapter3\';
elseif(ismac==1)
    foldersource='/Users/Shaun/Documents/GitHub/Chapter3/';
end

load(strcat(foldersource,'GrPrSpline.mat'))
load(strcat(foldersource,'ReNuSpline.mat'))
load(strcat(foldersource,'NuReSpline.mat'))

I=500/2;
Ta=21;
H=0;
D=0.029591;
phi=pi/2;
Vw=4;
alpha=5.3686e-5;
beta=2.6843e-7;
epsilons=0.7;
Psol=26.6319;
%conductorData=importfileAB(strcat(foldersource,'conductorData.csv'));
%polymodel=str2func(conductorData(strcmp(conductorData.CodeWord,"Rail"),:).polymodels);
polymodel=[];
[GuessTc,I2R,I2Rprime,Prad,PradPrime,Pcon,PconPrime] =GetTempNewton(I,Ta,H,D,phi,Vw,alpha,beta,epsilons,Psol,f,ff,ffinv,polymodel)