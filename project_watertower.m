%% Programming Project Water Tower Frame

      nnodes =  16;
      ndof = nnodes*6;

%coord in mm
coord = 1000*[0 0 0; 10 0 0; 10 10 0; 0 10 0; 0 0 10; 10 0 10; 10 10 10; 0 10 10; 0 0 20; 10 0 20; 10 10 20; 0 10 20; 0 0 30; 10 0 30; 10 10 30; 0 10 30];


%holds external nodal forces
concen = zeros(nnodes, 6);
concen(13, 3) = -2;
concen(14, 3) = -2;
concen(15, 3) = -2;
concen(16, 3) = -2;

%fixity shows dof where 0 is fixed dof and nan is free dof
%value would indicate settlement
fixity = nan(nnodes, 6);
fixity(1, 1:6) = [0 0 0 nan nan nan];
fixity(2, 1:6) = [0 0 0 nan nan nan];
fixity(3, 1:6) = [0 0 0 nan nan nan];
fixity(4, 1:6) = [0 0 0 nan nan nan];

      nele = 32;
      ends = zeros(nele,2);

      %end nodes of each element
      ends(1, 1:2) = [1 13];
      ends(2, 1:2) = [2 14];
      ends(3, 1:2) = [3 15];
      ends(4, 1:2) = [4 16];
      ends(5, 1:2) = [1 6];
      ends(6, 1:2) = [2 5];
      ends(7, 1:2) = [2 7];
      ends(8, 1:2) = [3 6];
      ends(9, 1:2) = [3 8];
      ends(10, 1:2) = [4 7];
      ends(11, 1:2) = [4 5];
      ends(12, 1:2) = [1 8];
      ends(13, 1:2) = [5 10];
      ends(14, 1:2) = [6 9];
      ends(15, 1:2) = [6 11];
      ends(16, 1:2) = [7 10];
      ends(17, 1:2) = [7 12];
      ends(18, 1:2) = [8 11];
      ends(19, 1:2) = [8 9];
      ends(20, 1:2) = [5 12];
      ends(21, 1:2) = [9 14];
      ends(22, 1:2) = [10 13];
      ends(23, 1:2) = [10 15];
      ends(24, 1:2) = [11 14];
      ends(25, 1:2) = [11 16];
      ends(26, 1:2) = [12 15];
      ends(27, 1:2) = [12 13];
      ends(28, 1:2) = [9 16];
      ends(29, 1:2) = [13 14];
      ends(30, 1:2) = [14 15];
      ends(31, 1:2) = [15 16];
      ends(32, 1:2) = [16 13];
      

%member properties
to_vector = ones(nele,1);
A = (1/4)*500^2*pi*to_vector; % mm^2  
A(1:4) = (1/4)*1000^2*pi; %area of four main members
I = ((pi*500^4)/64)*to_vector; %mm^4
I(1:4) = ((pi*1000^4)/64) %mm^4

J = (pi*500^4)/32*to_vector; %mm^4
J(1:4) = (pi*1000^4)/32; %mm^4

E = 200*to_vector; %Kn/mm
v = 0.3*to_vector; 
beta_ang = zeros(nele,1)';

%displacements
DEFL = zeros(nnodes, 6);

%dofs each node
node_id = zeros(nnodes,6);

for i = 1:nnodes
    for j = 1:6
        node_id(i,j) = (i-1)*6+j;
    end 
end

%dofs each member
mem_id = zeros(nele, 12);

for i = 1:nele
    mem_id(i,1:6) = node_id(ends(i,1),1:6);
    mem_id(i,7:12) = node_id(ends(i,2),1:6);
end


D = fixity';
D = D(:);

%find fixed dof and free dof
fixed_dof = find(D==0);
free_dof = find(isnan(D));

%Ktotal * deltatotal = Ptotal
P_total = concen';
P_total = P_total(:);

L = zeros(nele,1);

%get lenght of each member
for i = 1:nele
    x1 = coord(ends(i,1), 1);
    x2 = coord(ends(i,2), 1);
    y1 = coord(ends(i,1), 2);
    y2 = coord(ends(i,2), 2);
    z1 = coord(ends(i,1), 3);
    z2 = coord(ends(i,2), 3);

    L(i) = sqrt((x2-x1).^2 + (y2-y1).^2 + (z2-z1).^2);

end


%create xprime axis coordinates for each member
x_axis = zeros(nele,3);

for i = 1:nele
    x1 = coord(ends(i,1), 1);
    x2 = coord(ends(i,2), 1);
    y1 = coord(ends(i,1), 2);
    y2 = coord(ends(i,2), 2);
    z1 = coord(ends(i,1), 3);
    z2 = coord(ends(i,2), 3);

    x_axis(i,1) = (x2 - x1)/L(i);
    x_axis(i,2) = (y2 - y1)/L(i);
    x_axis(i,3) = (z2 - z1)/L(i);

end

%create three dimensional local stiffness matrices
k_stack_local = zeros(12,12,nele);

for i = 1:nele
    k_stack_local(1:12,1:12,i) = estiff(A(i), I(i), I(i), J(i), E(i), v(i), L(i));
end

k_trans_stack = zeros(12, 12, nele);

for i = 1:nele
    k_trans_stack(1:12,1:12,i) = etran(beta_ang(i),x_axis(i,:));
end

k_stack_global = zeros(12,12,nele);

for i = 1:nele
    k_stack_global(1:12,1:12,i) = (k_trans_stack(1:12,1:12,i))' * k_stack_local(1:12,1:12,i) * k_trans_stack(1:12,1:12,i);
end

k_total = zeros(ndof, ndof);

for i = 1:nele

    k_total(mem_id(i, 1:12), mem_id(i,1:12)) = k_total(mem_id(i, 1:12), mem_id(i,1:12)) + k_stack_global(1:12, 1:12, i);
end

kff = k_total(free_dof, free_dof);
ksf = k_total(fixed_dof, free_dof);

p_free = P_total(free_dof);

d_free = kff\p_free;

%insert dfree into dtotal
d_total = zeros(ndof, 1);

for i = 1:length(d_free)
    
    d_total(free_dof(i)) = d_free(i);
end

DEFL = reshape(d_total, [6, nnodes])';

react = ksf * d_free;

%add reactions to the P_total which has all the forces
for i = 1:size(react)
    P_total(fixed_dof(i)) = react(i);
end


%ELE_FOR matrix has element each row for each member and all the dofs for
%each element in the columns
ELE_FOR = zeros(nele, 12);

for i = 1:nele

    d_global = d_total(mem_id(i,1:12));

    d_local = k_trans_stack(1:12,1:12,i) * d_global;

    f_local = k_stack_local(1:12,1:12,i) * d_local;

    f_local = f_local';

    ELE_FOR(i, 1:12) = f_local;

end

%% plot struct
plotstructure(coord, ends)


%% dictionart of variables

%       node i's coordinates
%                           coord(i,1) = X coordinate
%                           coord(i,2) = Y coordinate
%                           coord(i,3) = Z coordinate
% concentrated loads for node i's 6 d.o.f.
%                            concen(i,1) = force in global X direction
%                            concen(i,2) = force in global Y direction
%                            concen(i,3) = force in global Z direction
%                            concen(i,4) = moment about global X axis
%                            concen(i,5) = moment about global Y axis
%                            concen(i,6) = moment about global Z axis
%   prescribed displacements for node i's 6 d.o.f.
%                          Note: A free d.o.f. will have a value of NaN
%                          and hence, you will find the Matlab function
%                          isnan very useful.
%                          Examples: If fixity(15,3) is set to NaN, then node 15's
%                                      Z-disp component is free;
%                                    If fixity(2,6) is set to 0.0, then node 2's
%                                      Z-rotation component is supported;
%                            fixity(i,1) = prescribed disp. in global X direction
%                            fixity(i,2) = prescribed disp. in global Y direction
%                            fixity(i,3) = prescribed disp. in global Z direction
%                            fixity(i,4) = prescribed rotation about global X axis
%                            fixity(i,5) = prescribed rotation about global Y axis
%                            fixity(i,6) = prescribed rotation about global Z axis


%       Izz(i)         ==  element i's moment of inertia about its local z-z axis
%       Iyy(i)         ==  element i's moment of inertia about its local y-y axis
%       J(i)           ==  element i's torsional constant
%       E(i)           ==  element i's material elastic modulus, Young's Modulus
%       v(i)           ==  element i's material Poisson's ratio
%       beta_ang(i)    ==  element i's web rotation angle.  These values are
%                          provided for those students who are required to calculate
%                          their own unit web vectors (see above).  It is based
%                          on the structure's undeformed geometry.
%                              Use the following convention for
%                                 defining a member's default web orientation:
%                                 A vector defing the element's local y-axis
%                                 with respect to the global coordinate system
%                                 will have a positive component in the global
%                                 Y direction.  If the element's local x-axis,
%                                 its length axis, is aligned with the global Y
%                                 axis, then element's local y-axis is aligned
%                                 with global negative X axis.  After this initial
%                                 orientation, element i may be rotated about
%                                 its local x-axis by the amount defined by
%                                 its web rotation angle, beta_ang(i).  The
%                                 angle is in radians and assumes a right-hand
%                                 convention about the local x-axis which runs from
%                                 the element's start node to its finish node.
% %

%     Output Information:
%       DEFL(i,1:6)      ==  node i's calculated 6 d.o.f. deflections
%                              DEFL(i,1) = displacement in X direction
%                              DEFL(i,2) = displacement in Y direction
%                              DEFL(i,3) = displacement in Z direction
%                              DEFL(i,4) = rotation about X direction
%                              DEFL(i,5) = rotation about Y direction
%                              DEFL(i,6) = rotation about Z direction
%       REACT(i,1:6)     ==  reactions for supported node i's 6 d.o.f.
%                              REACT(i,1) = force in X direction
%                              REACT(i,2) = force in Y direction
%                              REACT(i,3) = force in Z direction
%                              REACT(i,4) = moment about X direction
%                              REACT(i,5) = moment about Y direction
%                              REACT(i,6) = moment about Z direction
%       ELE_FOR(i,1:12)  ==  element i's internal forces and moments
%                            Note: All values reference the element's local
%                                  coordinate system.
%                              ELE_FOR(i,1)  = x-force at start node
%                              ELE_FOR(i,2)  = y-force at start node
%                              ELE_FOR(i,3)  = z-force at start node
%                              ELE_FOR(i,4)  = x-moment at start node
%                              ELE_FOR(i,5)  = y-moment at start node
%                              ELE_FOR(i,6)  = z-moment at start node
%                              ELE_FOR(i,7)  = x-force at end node
%                              ELE_FOR(i,8)  = y-force at end node
%                              ELE_FOR(i,9)  = z-force at end node
%                              ELE_FOR(i,10) = x-moment at end node
%                              ELE_FOR(i,11) = y-moment at end node
%                              ELE_FOR(i,12) = z-moment at end node
%       AFLAG            ==  logical flag to indicate if a successful
%                            analysis has been completed
%                              AFLAG = 1     Successful
%                              AFLAG = 0     Unstable Structure
%                              AFLAG = inf   No analysis code available
%
%


