%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IPA SLAM State Augmentation
% Author: Matthew Givens
% Date: August 24, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -Due to the Markov Property, only the prior dynamic states need be
%  included in x_km1 and R_km1 for state augmentation (disregard map states)
% -This assumes MEKF quaternion attitude and discrete-time IMU inputs.
% -States are in ECEF-like frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%   x_km1 - vehicle state vector prior to propagation
%   R_km1 - SR information matrix of x_km1
%   u - input vector of delta velocity and delta angle [dv_k;da_k]
%   dt - time increment
%   params - structure containting R_Q (inverse of process noise Q),
%            R_NB (DCM from body to CON frame), g (magnitude of gravity),
%            and omega_E (rotation rate of planet)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mu_out,R_out] = augment(x_km1,R_km1,u,dt,params)

% Extract inputs
dv_k = u(1:3,1);
da_k = u(4:6,1);
da_kp = u(7:9,1);

% Extract States
r = x_km1(1:3,1);
R_EN = q2C(x_km1(7:10,1));
R_EB = R_EN*params.R_NB;

% Process IMU outputs (sculling & coning)
dvN_k = params.R_NB*(dv_k-lever_arm_effect(da_k,da_kp,params.pB_IMU,dt));
g_k = -params.mu/norm(r)^3*r;
dvg_k = dt*(g_k-2*cross(params.omega_E,x_km1(4:6,1)));
dvc_k = R_EN*dvN_k + dvg_k;
dp_k = dt*x_km1(4:6,1) + dt/2*dvc_k;

% Update states with IMU data
x_new(1:6,1) = x_km1(1:6,1) + [dp_k;dvc_k];
dpsi = R_EN*params.R_NB*da_k;
dq = [dpsi./2;1];
x_new(7:10,1) = qXp(qinv(dq),x_km1(7:10,1));
x_new(7:10,1) = x_new(7:10,1)/norm(x_new(7:10,1));

% Populate STMs
Phi_x = [eye(3),dt*(eye(3)-dt*X(params.omega_E)),dt/2*X(R_EN*dvN_k);
        zeros(3),eye(3)-2*dt*X(params.omega_E),X(R_EN*dvN_k);
        zeros(3),zeros(3),eye(3)];
Gamma_x = [eye(3),-dt/2*R_EB,zeros(3);
           zeros(3),-R_EB,zeros(3);
           zeros(3),zeros(3),-params.R_NB];

% Build system and minimize
Q = [(1e-10)^2*eye(3),zeros(3,6);
    zeros(6,3),params.Q];
GQ = Gamma_x*Q*Gamma_x.';
R_Q = eye(params.m)/chol(GQ);
V1 = -R_Q*Phi_x;
R_pre = [R_km1,zeros(params.m);
         V1,R_Q];
[~,R_post] = qr(R_pre);

R_out = R_post;
mu_out(1:params.n,1) = x_km1(1:params.n,1);
mu_out(params.n+1:2*params.n,1) = x_new(:,1);

end