%% Calculate the external load at each joint
Joint_1 = calc_joint_force([6, 8, 12.8], 1);
Joint_2 = calc_joint_force([8, 8, 10, 6, 10, 12.8], 1);
Joint_3 = calc_joint_force([8, 10, 6, 10, 10], 1);
Joint_4 = calc_joint_force([6, 8, 10, 12.8], 1);
Joint_5 = calc_joint_force([6, 8, 8, 10, 10, 12.8], 1);
Joint_6 = calc_joint_force([6, 8, 10, 12.8], 1);
Joint_7 = calc_joint_force([10, 10, 13.8], 1);
Joint_8 = calc_joint_force([8, 10, 10, 12.8, 12.8], 2);
Joint_9 = calc_joint_force([8, 10, 10, 12.8, 12.8, 13.8], 3);
Joint_total = Joint_1 + Joint_2 + Joint_3 + Joint_4 + Joint_5 + Joint_6 + Joint_7 + Joint_8 + Joint_9;