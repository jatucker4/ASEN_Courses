function [F0,F1,F2,F3D,LVDT] = data_process(data_mat)
%AVG_DATA This function takes in the data matrices and outputs a matrix of
%the averaged data for each 
F0_raw = data_mat(:,2)*4.44822;
F1_raw = data_mat(:,3)*4.44822;
F2_raw = data_mat(:,4)*4.44822; 
F3D_raw = data_mat(:,5)*4.44822; 
LVDT_raw = data_mat(:,6)*0.0254;

F0(1) = mean(F0_raw(1:10,1));
F0(2) = mean(F0_raw(11:20,1));
F0(3) = mean(F0_raw(21:30,1));
F0(4) = mean(F0_raw(31:40,1));
F0(5) = mean(F0_raw(41:50,1));
F0(6) = mean(F0_raw(51:60,1));


F1(1) = mean(F1_raw(1:10,1));
F1(2) = mean(F1_raw(11:20,1));
F1(3) = mean(F1_raw(21:30,1));
F1(4) = mean(F1_raw(31:40,1));
F1(5) = mean(F1_raw(41:50,1));
F1(6) = mean(F1_raw(51:60,1));

F2(1) = mean(F2_raw(1:10,1));
F2(2) = mean(F2_raw(11:20,1));
F2(3) = mean(F2_raw(21:30,1));
F2(4) = mean(F2_raw(31:40,1));
F2(5) = mean(F2_raw(41:50,1));
F2(6) = mean(F2_raw(51:60,1));

F3D(1) = mean(F3D_raw(1:10,1));
F3D(2) = mean(F3D_raw(11:20,1));
F3D(3) = mean(F3D_raw(21:30,1));
F3D(4) = mean(F3D_raw(31:40,1));
F3D(5) = mean(F3D_raw(41:50,1));
F3D(6) = mean(F3D_raw(51:60,1));

LVDT(1) = mean(LVDT_raw(1:10,1));
LVDT(2) = mean(LVDT_raw(11:20,1));
LVDT(3) = mean(LVDT_raw(21:30,1));
LVDT(4) = mean(LVDT_raw(31:40,1));
LVDT(5) = mean(LVDT_raw(41:50,1));
LVDT(6) = mean(LVDT_raw(51:60,1));
end

