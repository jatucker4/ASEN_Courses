syms T_s

temp = solve(8==4.7945*.2*.15*(T_s-20) + .8*.2*.15*(5.67e-8)*((T_s+273)^4 - (20+273)^4));
temp = double(temp);

temp_1 = solve(8==7.317*.2*.15*(T_s-20) + .8*.2*.15*(5.67e-8)*((T_s+273)^4 - (20+273)^4));
temp_1 = double(temp_1);

temp_2 = solve(8==3.3486*.2*.15*(T_s-20) + .8*.2*.15*(5.67e-8)*((T_s+273)^4 - (20+273)^4));
temp_2 = double(temp_2);

temp_3 = solve(54==(7.854)*(.020106)*(T_s-25) + .9*(.020106)*(5.67e-8)*((T_s+273)^4 - (25+273)^4));
temp_3 = double(temp_3);