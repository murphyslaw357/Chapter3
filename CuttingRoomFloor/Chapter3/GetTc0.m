function [GuessTc] =GetTc0(I,Ta,D,R_T_high,R_T_low, T_high, T_low,epsilons,Psol)
    %I - RMS steady-state load current - amps
    %Ta - ambient temperature - degc
    %D - conductor diameter - meters
    %Psol - solar heating - w/m
    beta=(R_T_high-R_T_low)/(T_high-T_low);
    alpha=R_T_high-beta*T_high;
    sigmab=5.6697e-8;
    IIstar=abs(I)^2;
    GuessTc=((Psol+IIstar*(alpha+25*beta))/(pi*D*sigmab*epsilons)+((Ta+273)^4))^(1/4)-273;
end