function [x_eci,R] = ECEFtoECI(x,t,params)

r_ecef = x(1:3,1);
v_ecef = x(4:6,1);

theta = params.theta_f + params.omega_e(3)*t;

R = [cos(theta),-sin(theta),0;  
    sin(theta),cos(theta),0;
    0,0,1];

r_eci = R*r_ecef;
v_eci = R*v_ecef + cross(params.omega_e,R*r_ecef);
x_eci = [r_eci;v_eci];
end