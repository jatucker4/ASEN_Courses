function b = create_b(F, M, Fcoord)
%this function will create the b vector for the linear system of
%equilibrium equations, where b is (negative) sum of external force and moments

%F is a *x4 cell array with format F, dx, dy, dz
%M is a *x4 cell array with format M, ux, uy, uz
[rF,~] = size(F);
[rM,~] = size(M);

%b should be a 6x1 matrix with format Fx, Fy, Fz, Mx, My, Mz for external
%forces and moments
b = zeros(6,1);

%initialize necessary temporary variables
Ftemp = zeros(rF,3);
Mtemp = zeros(rM,3);
i = 1;
r = 1;

%fill Ftemp and Mtemp with true directional force/moment components from
%mag*dx,dy,dz
for r = 1:rF
    Ftemp(r,:) = F(r,1)*F(r,2:4); %multiplies magnitude of force by each component
end

for r = 1:rM
    Mtemp(r,:) = M(r,1)*M(r,2:4);%multiplies magnitude of moment by each component
end

b_temp = zeros(size(b));
%fill b_temp vector with the sum of the external forces and partial of sum of external moments
for i = 1:2
    if i == 1
        if rF == 1 %if only one external force, then no sum needed
            b_temp(1:3) = Ftemp;
        elseif rF ~= 1
            b_temp(1:3) = sum(Ftemp);
        end
    elseif i == 2
        if rM == 1 %if only one external moment, then no sum needed
            b_temp(4:6) = Mtemp;
        elseif rM ~= 1
            b_temp(4:6) = sum(Mtemp);
        end
    end
end

%external moments due to external forces - use vector determination of
%moments from the x, y, and z axis 
extM_fromF = zeros(rF,3);

for i = 1:rF  
    r = Fcoord(i,:); %radius from origin
    Force = F(i,1)*F(i,2:4); %true x, y, z components of force
    extM_fromF(i,:) = cross(r,Force); %cross produces vector of x,y,z components of moment
end

%splitting/summing moments x, y, z components
if rF == 1 %if only one force, no sum is needed
    extM_fromF = extM_fromF'; 
else
    extM_fromF = sum(extM_fromF)';
end

%for rows 1-3 (sum of forces), b_temp has correct values
%for row 4-5 (sum of moments), b_temp only has moments from the external
%moments and not moments produced by external forces
b =[b_temp(1:3);b_temp(4:6)+extM_fromF];

%negate b to simulate subtracting external forces/moments from 0 (on RHS)
b = -b;
end