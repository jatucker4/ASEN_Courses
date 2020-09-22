function FOS = factor_of_safety(desired)
f_dsr = icdf('normal', desired, 4.8, 0.4);
n = (4.8-f_dsr)/0.4;
FOS = 1/(1-(n*(0.4/4.8)));

