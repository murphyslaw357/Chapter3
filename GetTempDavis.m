function [GuessTc] =GetTempDavis(I,Ta,D,phi,Vw,R_T_high,R_T_low, T_high, T_low,epsilons,Psol)
    
    %I - RMS steady-state load current - amps
    %Ta - ambient temperature - degc
    %H - conductor elevation - meters
    %D - conductor diameter - meters
    %phi - angle between the wind direction and conductor axis - radians
    %Vw - Wind velocity - m/s
	%epsilons - conductor emissivity
    %Psol - solar heating - w/m  
    beta=(R_T_high-R_T_low)/(T_high-T_low);
    alpha=R_T_high-beta*T_high;
    
    sigmab=5.6697e-8;
    tolerance=0.15;
    error=realmax;
    GuessTcTop=400;
    GuessTcBottom=Ta;
    IIstar=abs(I)^2;

    J=0.138*D*epsilons/(10^8);
    W=4*273;
    Re=sin(phi)*Vw*D/vf;
    L=10^(a0+a1*log10(NRe)+a2*(log10(NRe)^2));
    p = [-J -J*W 3 -2 -4];
    r = roots(p)
end