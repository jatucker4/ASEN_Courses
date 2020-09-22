function delta = get_half_width(deficit_vec,y_vec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GET_HALF_WIDTH: 
% This function takes in velocity deficit and y location vectors and
% outputs the corresponding half width vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loc_plus_y = find(deficit_vec./max(deficit_vec) > 0.5,1, 'first');
loc_minus_y = find(deficit_vec./max(deficit_vec) > 0.5,1, 'last');

y_plus = y_vec(loc_plus_y);
y_minus = y_vec(loc_minus_y);

delta = 0.5.*(y_minus - y_plus);
end

