function [GuessTc] =GetTempStd738(I,Ta,D,phi,Vw,R_T_high,R_T_low, T_high, T_low,epsilons,Psol,H)
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
    tolerance=1;
    error=realmax;
    GuessTcTop=400;
    GuessTcBottom=Ta;
    IIstar=abs(I)^2;
    
    Kangle=1.194-cos(phi)+0.194*cos(2*phi)+0.368*sin(2*phi);
    counter=0;
    while(abs(error)>tolerance)
        counter=counter+1;
        GuessTc=(GuessTcTop+GuessTcBottom)/2;
        %qj=IIstar*(alpha+beta*GuessTc);
        qr=pi*D*sigmab*epsilons*(((GuessTc+273)^4)-((Ta+273)^4));
        Tfilm=(GuessTc+Ta)/2;
        muf=((1.458e-6)*((Tfilm+273)^1.5))/(Tfilm+383.4);
        rhof=(1.293-(1.525e-4)*H+(6.379e-9)*(H^2))/(1+0.00367*H);
        qcn=3.645*sqrt(rhof)*(D^0.75)*((GuessTc-Ta)^1.25);
        Nre=(D*rhof*Vw)/muf;
        kf=(2.424e-2)+(7.477e-5)*Tfilm-(4.407e-9)*(Tfilm^2);
        qc1=Kangle*(1.01+1.35*(Nre^0.52))*kf*(GuessTc-Ta);
        qc2=Kangle*(0.754*(Nre^0.6))*kf*(GuessTc-Ta);
        qc=max([qcn,qc1,qc2]);
        error=I-sqrt((qc+qr-Psol)/(alpha+beta*GuessTc));
        
        if(error>0)
            GuessTcBottom=GuessTc;
        else
            GuessTcTop=GuessTc;
        end
        if(counter>=5000)
            disp('error')
        end
    end

end