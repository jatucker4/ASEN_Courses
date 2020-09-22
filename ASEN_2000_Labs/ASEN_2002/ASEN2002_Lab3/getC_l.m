function [C_l_true,C_l_true_2,C_l_true_3,C_d_true_1,C_d_true_2,...
    C_d_true_3] = getC_l(cpData,x_position_of_ports,y_position_of_ports)
C_l_true = [];
C_d_true_1 = [];
[C_d,C_l] = drag_lift_coefficients(x_position_of_ports,...
    y_position_of_ports,cpData(1,3:19),cpData(1,1));
C_l_true = [C_l_true,C_l];
C_d_true_1 = [C_d_true_1,C_d];
[C_d,C_l] = drag_lift_coefficients(x_position_of_ports,...
    y_position_of_ports,cpData(4,3:19),cpData(4,1));
C_l_true = [C_l_true,C_l];
C_d_true_1 = [C_d_true_1,C_d];
[C_d,C_l] = drag_lift_coefficients(x_position_of_ports,...
    y_position_of_ports,cpData(7,3:19),cpData(7,1));
C_l_true = [C_l_true,C_l];
C_d_true_1 = [C_d_true_1,C_d];
[C_d,C_l] = drag_lift_coefficients(x_position_of_ports,...
    y_position_of_ports,cpData(10,3:19),cpData(10,1));
C_l_true = [C_l_true,C_l];
C_d_true_1 = [C_d_true_1,C_d];
% C_l_true = fliplr(C_l_true)*-1;

% next velocity
C_l_true_2 = [];
C_d_true_2 = [];
[C_d,C_l] = drag_lift_coefficients(x_position_of_ports,...
    y_position_of_ports,cpData(2,3:19),cpData(2,1));
C_l_true_2 = [C_l_true_2,C_l];
C_d_true_2 = [C_d_true_2,C_d];
[C_d,C_l] = drag_lift_coefficients(x_position_of_ports,...
    y_position_of_ports,cpData(5,3:19),cpData(5,1));
C_l_true_2 = [C_l_true_2,C_l];
C_d_true_2 = [C_d_true_2,C_d];
[C_d,C_l] = drag_lift_coefficients(x_position_of_ports,...
    y_position_of_ports,cpData(8,3:19),cpData(8,1));
C_l_true_2 = [C_l_true_2,C_l];
C_d_true_2 = [C_d_true_2,C_d];
[C_d,C_l] = drag_lift_coefficients(x_position_of_ports,...
    y_position_of_ports,cpData(11,3:19),cpData(11,1));
C_l_true_2 = [C_l_true_2,C_l];
C_d_true_2 = [C_d_true_2,C_d];
% C_l_true_2 = fliplr(C_l_true_2)*-1;
% Last velocity for the file
C_l_true_3 = [];
C_d_true_3 = [];
[C_d,C_l] = drag_lift_coefficients(x_position_of_ports,...
    y_position_of_ports,cpData(3,3:19),cpData(3,1));
C_l_true_3 = [C_l_true_3,C_l];
C_d_true_3 = [C_d_true_3,C_d];
[C_d,C_l] = drag_lift_coefficients(x_position_of_ports,...
y_position_of_ports,cpData(6,3:19),cpData(6,1));
C_l_true_3 = [C_l_true_3,C_l];
C_d_true_3 = [C_d_true_3,C_d];
[C_d,C_l] = drag_lift_coefficients(x_position_of_ports,...
    y_position_of_ports,cpData(9,3:19),cpData(9,1));
C_l_true_3 = [C_l_true_3,C_l];
C_d_true_3 = [C_d_true_3,C_d];
[C_d,C_l] = drag_lift_coefficients(x_position_of_ports,...
    y_position_of_ports,cpData(12,3:19),cpData(12,1));
C_l_true_3 = [C_l_true_3,C_l];
C_d_true_3 = [C_d_true_3,C_d];
% C_l_true_3 = fliplr(C_l_true_3)*-1;
end