function [GuessTc] =GetTempBlack(I,Ta,D,phi,Vw,R_T_high,R_T_low, T_high, T_low,epsilons,Psol)
    phi=90*pi/180;
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
    E=1.38e8+(1.39e6)*Ta;

    if(Vw==0)    
        counter=0;
        while(abs(error)>tolerance)
            counter=counter+1;
            GuessTc=(GuessTcTop+GuessTcBottom)/2;
            h=1.287*((GuessTc-Ta)/D)^0.25;
            error=GuessTc-(Ta+(IIstar*(alpha+beta*Ta)+Psol)/(pi*h*D+epsilons*pi*D*sigmab*E-IIstar*beta));
            if(error>0)
                GuessTcTop=GuessTc;
            else
                GuessTcBottom=GuessTc;
            end
        if(counter>=5000)
            disp('error')
        end
        end
    else
        h=(0.0272/D)*10^(2.217+0.652*log10(Vw*D)+0.0355*(log10(Vw*D)^2));
        
        GuessTc=Ta+(IIstar*(alpha+beta*Ta)+Psol)/(pi*h*D+epsilons*pi*D*sigmab*E-IIstar*beta);
    end

end