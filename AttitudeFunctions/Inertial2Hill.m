%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matthew Givens
% ASEN 6014
% Homework 2, Problem 14.1
% This function takes position and velocity inputs in the ECI frame and
% maps them into rectilinear Hill frame components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vec,ON] = Inertial2Hill(x_c,x_d)

r_cN = x_c(1:3,1);
rdot_cN = x_c(4:6,1);
r_dN = x_d(1:3,1);
rdot_dN = x_d(4:6,1);

o_r = r_cN/norm(r_cN);
o_h = cross(r_cN,rdot_cN)/norm(cross(r_cN,rdot_cN));
o_theta = cross(o_h,o_r);
fdot = norm(cross(r_cN,rdot_cN))/norm(r_cN)^2;
omega_ON = fdot*[0;0;1];

ON = [o_r,o_theta,o_h]';

rho = ON*(r_dN-r_cN);
o_rho = rho/norm(rho);
rhodot = ON*(rdot_dN-rdot_cN) - cross(omega_ON,rho);

vec = [rho;rhodot];

end