%Dr. Johnson's FEM Code
%Author: Aaron Johnson
%
%FUNCTION INPUTS
%
%connecttable describes the ELEMENTS
%   1: element number
%   2: node A
%   3: node B
%   4: cross-sectional area
%   5: elastic modulus
%   6: coefficient of thermal expansion (put 0 for all members)
%   7: temperature change (put 0 for all members)
%
%nodetable describes the NODES
%   1: node number
%   2: x-coordinate
%   3: y-coordinate
%   4: z-coodinate
%   5: x-displacement (number if known, NaN if unknown)
%   6: y-displacement (number if known, NaN if unknown)
%   7: z-displacement (number if known, NaN if unknown)
%   8: external force in x-direction (number if known, NaN if unknown)
%   9: external force in y-direction (number if known, NaN if unknown)
%   10: external force in z-direction (number if known, NaN if unknown)
%
%dof is the number of degrees of freedom in the truss (put 3 for ASEN 3112)
%
%
%FUNCTION OUTPUTS
%u = matrix of displacements for each node (size n x dof, where n is the number of nodes and dof is the number of degress of freedom)
%   For row n:
%   Column 1: x-displacement of node n
%   Column 2: y-displacement of node n
%   Column 3: z-displacement of node n
%
%F = matrix of external force on each node (size n x dof)
%   For row n:
%   Column 1: x-component of external force on node n
%   Column 2: y-component of external force on node n
%   Column 3: z-component of external force on node n
%
%sigma = vector of axial stress in each element (size m x 1, where m is the number of elements)
%   For row m:
%   Column 1: axial stress in element m

%function [u, F, sigma] = johnson_FEM_code(nodetable, connecttable, dof)
tic
%%%%%%%%%%%%PREPROCESSING%%%%%%%%%%%
% coordinate matrix for each node with the specified dof
co = nodetable(:,2:dof+1);

% external force matrix for each node with the specified dof
F = nodetable(:,8:dof+7);

% displacement matrix for each node with the specified dof
u = nodetable(:,5:dof+4);

Nel = size(connecttable,1); % number of elements
Nnodes = size(co,1);        % number of nodes

%%% initialize global stiffness matrix 'K' and force vector 'F'
K = zeros(Nnodes*dof,Nnodes*dof);
%Ftemp = zeros(Nnodes*dof,1);    % Column vector
%Fthermal = Ftemp;               % Column vector
Ftemp = reshape(F',Nnodes*dof,1);   % Column vector
Fthermal = zeros(Nnodes*dof,1);     % Column vector
sigma = nan(Nel,1);

%%% Assemble global K
% For each ELEMENT
for A = 1:Nel                   % for each element
    n = (co(connecttable(A,3),:) - co(connecttable(A,2),:));  % vector from node 2 to node 1
    L = norm(n);                % length of element
    n = n./L;                   % unit vector        
    Area = connecttable(A,4);   % area of element
    E = connecttable(A,5);      % elastic modulus of element
    alpha = connecttable(A,6);  % coefficient of thermal expansion
    delT = connecttable(A,7);   % temperature change
 
    k11 = (E*Area/L)*(n'*n);    % khat matrix, always 3x3
     
    % local stiffness matrix and force vector
    % localstiffness = [k11 -k11;-k11 k11];
   
    n1end = connecttable(A,2)*dof;    % Node 1 index end
    n1start = n1end-dof+1;            % Node 1 index start
    n2end = connecttable(A,3)*dof;    % Node 2 index end
    n2start = n2end-dof+1;            % Node 2 index start

    K(n1start:n1end, n1start:n1end) = K(n1start:n1end, n1start:n1end) + k11;
    K(n2start:n2end, n2start:n2end) = K(n2start:n2end, n2start:n2end) + k11;
    K(n1start:n1end, n2start:n2end) = K(n1start:n1end, n2start:n2end) - k11;
    K(n2start:n2end, n1start:n1end) = K(n2start:n2end, n1start:n1end) - k11;
    
    Fthermal(n1start:n1end) = Fthermal(n1start:n1end) + alpha*delT*E*Area*n';
    Fthermal(n2start:n2end) = Fthermal(n2start:n2end) - alpha*delT*E*Area*n';
end

Kfull = K;

%%% Incorporate boundary conditions
% For each NODE
for A=Nnodes:-1:1               % Work backwards so as to not mess up indexing
    nend = A*dof;               % Last index for this node
    for B=dof:-1:1              % Check each dof
        if isnan(F(A,B))        % If the force is NaN (unknown) at this node in this dof
            K(nend-dof+B,:) = [];               % remove row of K
            Ftemp(nend-dof+B,:) = [];           % remove row of Ftemp
            Kcolumn = K(:,nend-dof+B).*u(A,B);  % store column of K multiplied by known displacement     W19
            K(:,nend-dof+B) = [];               % remove column of K
            Ftemp = Ftemp - Kcolumn;            % modify Ftemp by the column of K x displacement         W19
        else                    % if the force is known at this node in this dof
            Ftemp(nend-dof+B) = Ftemp(nend-dof+B) - Fthermal(nend-dof+B); % add thermal effects
        end
    end
end
%     
%     if isnan(F(A,1))                % Nodes where force is Nan (unknown)
%         K(nstart:nend,:) =[];           % remove rows of K
%         K(:,nstart:nend) =[];           % remove columns of K
%         Ftemp(nstart:nend,:) = [];      % remove rows of Ftemp
%     else                            % Nodes where force is known
%         Ftemp(nstart:nend) = F(A,:)';   % populate rows of Ftemp
%     end
% end

% Solve for unknown displacements u
% utemp = K\Ftemp;
utemp = K^(-1)*Ftemp;

% Reassemble u vector
index = 1;                          % Start an index to pick off elements of the utemp vector
for A=1:Nnodes                      % For each node
    for B=1:dof
        if isnan(u(A,B))            % If the displacement of this node is unknown
            u(A,B) = utemp(index);  % Pick off the next elements of the utemp vector
            index = index+1;        % Advance the index
        end
    end
end

% Calculate reaction forces
u2 = u'; u2 = u2(:);                        % Reshape the u vector back into a column vector
Ftemp2 = Kfull*u2 + Fthermal;               % Calculate the column vector of all external forces
F = reshape(Ftemp2,[dof,Nnodes]); F = F';   % Reshape the F vector

% Calculate internal forces
for A = 1:Nel                   % for each element
    n = (co(connecttable(A,3),:) - co(connecttable(A,2),:));    % vector from node 2 to node 1
    delta = u(connecttable(A,3),:) - u(connecttable(A,2),:);   % elongation of bar (node 2 - node 1)
    L = norm(n);                % length of element
    n = n./L;                   % unit vector        
    Area = connecttable(A,4);   % area of element
    E = connecttable(A,5);      % elastic modulus of element
    alpha = connecttable(A,6);  % coefficient of thermal expansion
    delT = connecttable(A,7);   % temperature change
    P = (E*Area/L)*dot(delta,n) - alpha*delT*E*Area;
    sigma(A) = P/Area;
end
toc