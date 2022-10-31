%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IPA SLAM State Propagation + Marginalization
% Author: Matthew Givens
% Date: August 24, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -This assumes MEKF quaternion attitude and discrete-time IMU inputs.
% -States are in ECEF-like frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%   mu_km1 - vehicle state vector prior to propagation
%   R_window - SR information matrix including active map states
%   u - input vector of delta velocity and delta angle [dv_k;da_k]
%   dt - time increment
%   params - structure containting R_Q (inverse of process noise Q),
%            R_NB (DCM from body to CON frame), g (magnitude of gravity),
%            and omega_E (rotation rate of planet)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mu_out,R_out] = propagate(mu_km1,R_window,u,dt,params)

m = size(R_window,1);
Phi_k = eye(m);
Gamma_k = zeros(m,params.q);

% Extract inputs
dv_k = u(1:3,1);
da_k = u(4:6,1);
da_kp = u(7:9,1);

% Extract States
x = mu_km1(end-9:end,1);
p = x(1:3,1);
r = params.pB_IMU;
R_EN = q2C(x(7:10,1));
R_EB = R_EN*params.R_NB;

% Process IMU outputs (sculling & coning)
dvN_k = params.R_NB*(dv_k-lever_arm_effect(da_k,da_kp,r,dt));
g_k = -params.mu/norm(p)^3*p;
dvg_k = dt*(g_k-2*cross(params.omega_E,x(4:6,1)));
dvc_k = R_EN*dvN_k + dvg_k;
dp_k = dt*x(4:6,1) + dt/2*dvc_k;

% Update states with IMU data
x_new(1:6,1) = x(1:6,1) + [dp_k;dvc_k];
dpsi = R_EN*params.R_NB*da_k;
dq = [dpsi./2;1];
x_new(7:10,1) = qXp(qinv(dq),x(7:10,1));
x_new(7:10,1) = x_new(7:10,1)/norm(x_new(7:10,1));

% Populate STMs
Phi_x = [eye(3),dt*(eye(3)-dt*X(params.omega_E)),dt/2*X(R_EN*dvN_k);
        zeros(3),eye(3)-2*dt*X(params.omega_E),X(R_EN*dvN_k);
        zeros(3),zeros(3),eye(3)];
Gamma_x = [-dt/2*R_EB,zeros(3);
           -R_EB,zeros(3);
           zeros(3),-params.R_NB];

Phi_k(end-params.m+1:end,end-params.m+1:end) = Phi_x;
Gamma_k(end-params.m+1:end,:) = Gamma_x;

% Build system and minimize
R_tilde = R_window/Phi_k;
R_pre = [params.R_Q,zeros(params.q,size(R_tilde,1));R_tilde*Gamma_k,R_tilde];
[~,R_post] = qr(R_pre);

R_out = R_post(params.q+1:end,params.q+1:end);
mu_out = [mu_km1(1:end-params.n,1);x_new];

end