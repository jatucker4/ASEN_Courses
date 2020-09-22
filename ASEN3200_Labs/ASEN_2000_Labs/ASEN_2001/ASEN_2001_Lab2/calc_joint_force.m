function joint_force = calc_joint_force(lengths, sleeve_count)
for i = 1:length(lengths)
    lengths(i) = lengths(i)*0.0254*0.04767*9.81;
end
magnet_weight = (0.0018*9.81);
ball_weight = 0.0084*9.81;
sleeve_weight = 0.00535*9.81;
joint_force = ((sleeve_count*sleeve_weight) + (sum(lengths)))/2 + magnet_weight + ball_weight; 
end
