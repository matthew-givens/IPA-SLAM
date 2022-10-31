%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matthew Givens
% ASEN 6014
% Homework 2, Problem 14.1
% This function takes position and velocity inputs in the Hill frame and
% maps them into the ECI frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x_ECI = Hill2Inertial(x_H,x_chief)

r_cN = x_chief(1:3,1);
rdot_cN = x_chief(4:6,1);

rho = x_H(1:3,1);
rhodot = x_H(4:6,1);

o_r = r_cN/norm(r_cN);
o_h = cross(r_cN,rdot_cN)/norm(cross(r_cN,rdot_cN));
o_theta = cross(o_h,o_r);
fdot = norm(cross(r_cN,rdot_cN))/norm(r_cN)^2;
omega_ON = fdot*[0;0;1];

ON = [o_r,o_theta,o_h]';
NO = ON';

x_ECI(1:3,1) = r_cN + NO*rho;
x_ECI(4:6,1) = rdot_cN + NO*(rhodot+cross(omega_ON,rho));

end