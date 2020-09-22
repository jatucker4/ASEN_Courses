function velocity = pitot_static(pressure_change, T_atm, P_atm)
%pitot_static function calculates the exit velocity given a pressure change
%measured from a pitot static probe
velocity = sqrt(2*pressure_change*(287*T_atm/P_atm));
end

