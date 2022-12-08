function[elk] = estiff(A,Izz,Iyy,J,E,v,L)


%rows of element stiffness matrix rows 1 to 12
r1 = [A/L 0 0 0 0 0 -A/L 0 0 0 0 0];
r2 = [0 (12*Izz)/(L.^3) 0 0 0 (6*Izz)/(L.^2) 0 (-12*Izz)/(L.^3) 0 0 0 (6*Izz)/(L.^2)];
r3 = [0 0 (12*Iyy)/L.^3 0 (-6*Iyy)/(L.^2) 0 0 0 (-12*Iyy)/(L.^3) 0 (-6*Iyy)/(L.^2) 0];
r4 = [0 0 0 J/(2*(1+v)*L) 0 0 0 0 0 -J/(2*(1+v)*L) 0 0];
r5 = [0 0 -(6*Iyy)/L.^2 0 (4*Iyy)/L 0 0 0 (6*Iyy)/L.^2 0 (2*Iyy)/L 0];
r6 = [0 (6*Izz)/L.^2 0 0 0 (4*Izz)/L 0 -(6*Izz)/L.^2 0 0 0 (2*Izz)/L];
r7 = [-A/L 0 0 0 0 0 A/L 0 0 0 0 0];
r8 = [0 -(12*Izz)/(L.^3) 0 0 0 -(6*Izz)/(L.^2) 0 (12*Izz)/(L.^3) 0 0 0 (-6*Izz)/(L.^2)];
r9 = [0 0 (-12*Iyy)/L.^3 0 (6*Iyy)/(L.^2) 0 0 0 (12*Iyy)/(L.^3) 0 (6*Iyy)/(L.^2) 0];
r10 = [0 0 0 -J/(2*(1+v)*L) 0 0 0 0 0 J/(2*(1+v)*L) 0 0];
r11 = [0 0 -(6*Iyy)/L.^2 0 (2*Iyy)/L 0 0 0 (6*Iyy)/L.^2 0 (4*Iyy)/L 0];
r12 = [0 (6*Izz)/L.^2 0 0 0 (2*Izz)/L 0 -(6*Izz)/L.^2 0 0 0 (4*Izz)/L];

elk = E*[r1; r2; r3; r4; r5; r6; r7; r8; r9; r10; r11; r12];
end
