function disp = beamApprox()
    %% constants
    rho = 2700;                 %kg/m^3; density of 6061-T6
    E = 69e9;                   %Pa; Young's Modulus of 6061-T6
    L_b = 0.25;                 %m; length of single bay
    L_beam = L_b*16;            %m; length of entire beam
    P = 222.4;                  %N; load applied
    
    %% Moment of inertia
    % calculate mass of individual strut
    R = ((3/8) / 2) * 0.0254;   %m; outer radius
    r = R - ((1/16)*0.0254);    %m; inner radius
    Area = pi*(R^2 - r^2);      %m^2; area of single cross section
    V = Area*L_beam;            %m^3; volume of individual long strut
    m = rho * V;                %kg; mass of individual strut
    
    I = 4*Area*(L_b/2)^2;       %m^4
    
    %% Beam model [conservation of energy]
    R_a = P/2;                  % left support reaction
    syms x
    M1 = R_a*x;                 % internal moment left half of beam
    M2 = R_a*(x + (L_beam/2)) - P*(x) ; % internal moment right half of beam
    
    % internal energy calculations
    A = M1^2 / (E*I);
    B = M2^2 / (E*I);
    U = 0.5*int(A,x,0,L_beam/2) + 0.5*int(B,x,L_beam/2,L_beam);
    
    % displacement at center [m]
    disp = double((-2*U)/P);

end