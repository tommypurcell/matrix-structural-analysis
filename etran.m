%--------------------------
%Thomas Purcell
%etran edited for project 
%--------------------------
% so xaxis should be edoubleprimex already


function[gamma] = etran(beta_ang, xaxis)

%---
%this is the test member from lecture 9 to confirm it works
% xaxis = [0.8639 0.4319 0.2529];
% beta_ang = 30;
%---

if xaxis == [0 1 0]
    ey = transpose([-1 0 0]);
elseif xaxis == [0 -1 0]
     ey = transpose([1 0 0]);
else 
    ey = transpose([0 1 0]);
end


ez_double_prime = cross(xaxis, ey)/(norm(cross(xaxis, ey)));

ey_double_prime = cross(ez_double_prime, xaxis); %/(norm(cross(xaxis, ez_double_prime)));


% R1 = [
%     dot(xaxis, ex) dot(xaxis, ey) dot(xaxis, ez)
%     dot(ey_double_prime, ex) dot(ey_double_prime, ey) dot(ey_double_prime, ez)
%     dot(ez_double_prime, ex) dot(ez_double_prime, ey) dot(ez_double_prime, ez)]

R1 = [xaxis; ey_double_prime; ez_double_prime];

R2 = [1 0 0; 0 cosd(beta_ang) sind(beta_ang); 0 -sind(beta_ang) cosd(beta_ang)];

small_gamma = R2*R1;

zero_matrix = zeros(3);

gamma = [small_gamma zero_matrix zero_matrix zero_matrix;
            zero_matrix small_gamma zero_matrix zero_matrix;
            zero_matrix zero_matrix small_gamma zero_matrix;
            zero_matrix zero_matrix zero_matrix small_gamma];

end
