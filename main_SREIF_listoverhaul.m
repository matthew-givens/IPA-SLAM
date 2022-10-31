clear
clc
close all

% Main SREIF Script

% Load data
addpath('NewAttitudeLibrary')
addpath('Data')
addpath('FilterFunctions')
load('measurement_data_L=6_SURF.mat')
maxI = 410;
load('distortion_jacobians.mat')
params.H_dist = H_dist;
params.H_undist = H_undist;

VO_updates = 1;
params.n = 10; % Navigation state size
params.m = 9; % Navigation error state size
params.postfit = 1;
params.debug_image = 0;
params.n_max = 3;
params.num_iter = 4;
params.gating = 1;
params.beta = 0.5;

% Load IMU data
true_table = readtable('Data/truth.csv','VariableNamingRule','preserve');
true_array = table2array(true_table);
imu_table = readtable('Data/dlc.csv','VariableNamingRule','preserve');
imu_array = table2array(imu_table);

imu_array(:,1) = imu_array(:,1)/1e9;
t_start = imu_array(1,1);
imu_array(:,1) = imu_array(:,1) - t_start;
true_array(:,1) = true_array(:,1)/1e9;
true_array(:,1) = true_array(:,1) - t_start;

%%% Frames %%%
% E is ECEF frame
% B is IMU frame
% N is Center of Navigation (CON) frame

%%% Simulation Parameters %%%
params.omega_E = [0;0;2*pi/(24*3600)];
params.g = norm((sum(imu_array(1:10351,2:4))/10351)/0.02);
params.r_E = norm(true_array(1,2:4)); % Launch Earth radius [m]
params.mu = params.g*params.r_E^2;
params.theta_f = 0; % Earth initial sidereal time [rad] (arbitrary)

params.pN_IMU = [-0.08035, 0.28390, -1.42333 ].';
params.R_NB = [-0.2477 -0.1673  0.9543;
               -0.0478  0.9859  0.1604;
               -0.9677 -0.0059 -0.2522];
params.pB_IMU = params.R_NB.'*params.pN_IMU;

i_start = 1;
t = imu_array(i_start:end,1);
tspan = length(t);
dv = imu_array(i_start:end,1);

xE = zeros(10,tspan);
xE(:,1) = true_array(1,2:11).';
xE_true = zeros(10,tspan);
xE_true(:,1) = xE(:,1);

Pnav = zeros(9,9,tspan);
sigma_p = 10*ones(1,3);
sigma_v = 2*ones(1,3);
sigma_a = 0.05*ones(1,3);
Pnav(:,:,1) = diag([sigma_p.^2,sigma_v.^2,sigma_a.^2]);
sigma_dv = 0.0005;
sigma_theta = 0.00003;
params.Q = [sigma_dv^2*eye(3),zeros(3);
    zeros(3),sigma_theta^2*eye(3)];
params.R_Q = inv(chol(params.Q));

% Camera Sensor
sigma_m = 0.5;
params.q = size(params.Q,1);
params.V = sigma_m^2*eye(3);
params.R_NC = [0.1742 -0.9500 0.2590
              -0.9847 -0.1681 0.0457
               0.0001 -0.2630 -0.9648];
params.pN_C = [-0.09649, 0.27776, -1.54774].';
params.R_CB = params.R_NC.'*params.R_NB;
params.K = [4794.313768 0.0         720.027670;
            0.0         4800.628517 539.999614;
            0.0         0.0         1.0];
params.f = [params.K(1,1) params.K(2,2)];
params.c = [params.K(1,3) params.K(2,3)];
params.k1 = -7.577065e-02;
params.k2 = 2.939401e-02;
params.k3 = -2.308007e-01;
k = [params.k1 params.k2 params.k3];
params.p1 = -6.624078e-03;
params.p2 = 3.305442e-04;
p = [params.p1 params.p2];
params.intrinsics = cameraIntrinsics(params.f,params.c,[1080 1440],'RadialDistortion',k,...
    'TangentialDistortion',p);
params.Kinv = inv(params.K);
params.V = sigma_m^2*eye(2);
params.R_C = inv(chol(params.V));

% Load images
imageNames = dir(fullfile('Images','image_*.000000.png'));
imageNames = {imageNames.name}';
pat = digitsPattern;
imageData = extract(imageNames,pat);
imageTimes = str2double(imageData(:,1))./1e9 - t_start;
j = 1;
tol = 5e-3;
prevI = [];
points = cell(1,maxI);

rho_true = zeros(1,tspan);
rho = zeros(1,tspan);
rho_true(1) = norm(true_array(1,2:4));
rho(1) = norm(xE(1:3,1));

%%%% SREIF Preliminaries %%%%
% muSize is the size of the map for the full vector, Rsize is the size of
% the map for the error vector
next_idx = 1;
Rnext_idx = 1;

mu = zeros(params.n+3000,tspan); % Vehicle + landmark states
R = zeros(params.n-1+3000,params.n-1+3000); % All covariances (sliding window)
mu(1:params.n,1) = xE(:,1);
R(1:params.n-1,1:params.n-1,1) = chol(eye(params.n-1)/Pnav(:,:,1));
r = cell(1,maxI);

% 'lists' is a nested structure that contains cell arrays for both features
% and parents as well as the current indices of the navigation state and
% the identifiers of all currently active landmarks. Features and parents
% are organized by the order of their instantiation and thus their
% identifiers correspond to their indicies in the cell array
lists.nav_state.idx = 1:params.n;
lists.nav_state.Ridx = 1:params.m;
lists.features = {};
lists.parent = {};
lists.parent_update = [];
lists.n_parent = 0;
lists.active = [];
lists.mu_prop_window = lists.nav_state.idx;
lists.R_prop_window = lists.nav_state.Ridx;
lists.mu_update_window = lists.mu_prop_window;
lists.R_update_window = lists.R_prop_window;
lists.new = [];
lists.deactivate = [];

I = 0;
new_feature = 0;

for k = 2:tspan

    if t(k) < 207 + 25

        xE_true(:,k) = true_array(2*k,2:11).';
        xE(:,k) = xE_true(:,k);
        mu(1:params.n,k) = xE(:,k);
        rho(k) = norm(xE(1:3,k));
        rho_true(k) = norm(true_array(2*k,2:4));
        Pnav(:,:,k) = Pnav(:,:,k-1);

    else

        % SREIF Propagation
        dt = t(k) - t(k-1);
        u(1:3,1) = imu_array(k-1,2:4).';
        u(4:6,1) = imu_array(k-1,5:7).';
        u(7:9,1) = imu_array(k,5:7).';

        mu(1:lists.nav_state.idx(end),k) = mu(1:lists.nav_state.idx(end),k-1);

        if new_feature == 0 % Propagate via Marginalization

            lists.mu_prop_window = lists.mu_update_window;
            lists.R_prop_window = lists.R_update_window;
            mu_in = mu(lists.mu_prop_window,k-1);
            R_in = R(lists.R_prop_window,lists.R_prop_window);

            [mu_out,R_out] = IPASLAM_propagate(mu_in,R_in,u,dt,params);

        else % Propagate via Augmentation

            lists.mu_prop_window = lists.nav_state.idx;
            lists.R_prop_window = lists.nav_state.Ridx;

            lists.n_parent = lists.n_parent + 1;
            lists.parent{lists.n_parent}.idx = lists.nav_state.idx;
            lists.parent{lists.n_parent}.Ridx = lists.nav_state.Ridx;
            lists.parent_update = [lists.parent_update;lists.parent{lists.n_parent}.idx];
            
            lists.nav_state.idx = lists.nav_state.idx + params.n;
            lists.nav_state.Ridx = lists.nav_state.Ridx + params.m;

            mu_in = mu(lists.mu_prop_window,k-1);
            R_in = R(lists.R_prop_window,lists.R_prop_window);

            [mu_out,R_out] = IPASLAM_augment(mu_in,R_in,u,dt,params);

            lists.mu_prop_window = [lists.mu_prop_window,lists.mu_prop_window + params.n];
            lists.mu_update_window = lists.mu_update_window(1):lists.mu_update_window(end) + params.n;
            lists.R_prop_window = [lists.R_prop_window,lists.R_prop_window+params.m];
            lists.R_update_window = lists.R_update_window(1):lists.R_update_window(end)+params.m;

            next_idx = lists.nav_state.idx(1);
            Rnext_idx = lists.nav_state.Ridx(1);
            new_feature = 0;

        end

        mu(lists.mu_prop_window,k) = mu_out;
        R(lists.R_prop_window,lists.R_prop_window) = R_out;

        xE(:,k) = mu_out(end-params.n+1:end,1);
        Snav = eye(params.m)/R_out(end-params.m+1:end,end-params.m+1:end);
        Pnav(:,:,k) = Snav*Snav.';
        rho(k) = norm(xE(1:3,k));

        % Update
        % Time check for update
        t_diff = abs(imageTimes(j) - t(k));
        if t_diff <= tol && VO_updates == 1 && j < 350

            if sum(z_ID(:,j)) ~= 0

                % Data Association
                lists = data_association_listoverhaul(lists,z_ID(:,j),next_idx,Rnext_idx,params);
                new_feature = any(lists.new);
                
                min_update = min(lists.active);
                if isempty(lists.deactivate)
                    lists.mu_update_window = lists.features{min_update}.idx(1):lists.nav_state.idx(end);
                    lists.R_update_window = lists.features{min_update}.Ridx(1):lists.nav_state.Ridx(end);
                    lists.min_parent = lists.features{min_update}.parent;
                else
                    min_deactivate = min(lists.deactivate);
                    if min_update < min_deactivate
                        lists.mu_update_window = lists.features{min_update}.idx(1):lists.nav_state.idx(end);
                        lists.R_update_window = lists.features{min_update}.Ridx(1):lists.nav_state.Ridx(end);
                        lists.min_parent = lists.features{min_update}.parent;
                    else
                        lists.mu_update_window = lists.features{min_deactivate}.idx(1):lists.nav_state.idx(end);
                        lists.R_update_window = lists.features{min_deactivate}.Ridx(1):lists.nav_state.Ridx(end);
                        lists.min_parent = lists.features{min_deactivate}.parent;
                    end
                end
                
                mu_in = mu(lists.mu_update_window,k);
                R_in = R(lists.R_update_window,lists.R_update_window);

                lists.mu_offset = lists.mu_update_window(1) - 1;
                lists.R_offset = lists.R_update_window(1) - 1;

                if params.debug_image == 1
                    I = imread(fullfile('Images',imageNames{j}));
                end

                [x_plus,R_plus,r{j},lists,points{j}] = IPASLAM_update_listoverhaul(mu_in,R_in,z(:,j),...
                    lists,params,I);
                mu(lists.mu_update_window,k) = x_plus;
                xE(:,k) = x_plus(end-params.n+1:end,1);
                R(lists.R_update_window,lists.R_update_window) = R_plus;
                Snav = eye(params.m)/R_plus(end-params.m+1:end,end-params.m+1:end);
                Pnav(:,:,k) = Snav*Snav.';

            end

            % Iterate
            j = j + 1;
            k_final = k;

        end

    end

    if 2*k <= 79998

        xE_true(1:6,k) = true_array(2*k,2:7).';
        xE_true(7:10,k) = true_array(2*k,8:11).';
        pE_true(:,k) = ECEFtoENU(xE(1:3,1),xE_true(1:3,k),params,'vector');

        rho_true(k) = norm(true_array(2*k,2:4));

    else

        xE_true(1:6,k) = xE_true(1:6,k-1);
        xE_true(7:10,k) = xE_true(7:10,k-1);
        pE_true(:,k) = pE_true(:,k-1);
        rho_true(k) = rho_true(k-1);

    end

    percent = k*100/tspan;
    clc
    fprintf('Percent Complete: %3.2f\n',percent)

end
%%
t_elapsed = t - 207;
dx = xE_true(1:6,:) - xE(1:6,:);
for k = 1:tspan
    dq = qXp(xE_true(7:10,k),qinv(xE(7:10,k)));
    dx(7:9,k) = 2*dq(1:3)*180/pi;
end

pxstd = extract_nSigma(Pnav,1,3);
pystd = extract_nSigma(Pnav,2,3);
pzstd = extract_nSigma(Pnav,3,3);
vxstd = extract_nSigma(Pnav,4,3);
vystd = extract_nSigma(Pnav,5,3);
vzstd = extract_nSigma(Pnav,6,3);
a1std = extract_nSigma(Pnav,7,1)*180/pi;
a2std = extract_nSigma(Pnav,8,1)*180/pi;
a3std = extract_nSigma(Pnav,9,1)*180/pi;

%%
figure
plot(t_elapsed,(rho_true-rho(1))./1000,'--','Linewidth',2)
hold on
plot(t_elapsed,(rho-rho(1))./1000)
xlabel('Time [s]')
ylabel('Altitude [km]')
legend('Integrated','Truth','Location','Northwest')
grid on
if VO_updates == 1
    xline(t_elapsed(k_final))
end

%%
figure
subplot(3,1,1)
plot(t_elapsed,dx(1,:).'); hold on;
plot(t_elapsed,pxstd,'r')
plot(t_elapsed,-pxstd,'r')
ylabel('\Delta x [m]')
grid on
if VO_updates == 1
    xline(t_elapsed(k_final))
end
subplot(3,1,2)
plot(t_elapsed,dx(2,:).'); hold on;
plot(t_elapsed,pystd,'r')
plot(t_elapsed,-pystd,'r')
ylabel('\Delta y [m]')
grid on
if VO_updates == 1
    xline(t_elapsed(k_final))
end
subplot(3,1,3)
plot(t_elapsed,dx(3,:).'); hold on;
plot(t_elapsed,pzstd,'r')
plot(t_elapsed,-pzstd,'r')
ylabel('\Delta z [m]')
grid on
if VO_updates == 1
    xline(t_elapsed(k_final))
end

figure
subplot(3,1,1)
plot(t_elapsed,dx(4,:).'); hold on;
plot(t_elapsed,vxstd,'r')
plot(t_elapsed,-vxstd,'r')
ylabel('\Delta v_x [m]')
grid on
if VO_updates == 1
    xline(t_elapsed(k_final))
end
subplot(3,1,2)
plot(t_elapsed,dx(5,:).'); hold on;
plot(t_elapsed,vystd,'r')
plot(t_elapsed,-vystd,'r')
ylabel('\Delta v_y [m]')
grid on
if VO_updates == 1
    xline(t_elapsed(k_final))
end
subplot(3,1,3)
plot(t_elapsed,dx(6,:).'); hold on;
plot(t_elapsed,vzstd,'r')
plot(t_elapsed,-vzstd,'r')
ylabel('\Delta v_z [m]')
grid on
if VO_updates == 1
    xline(t_elapsed(k_final))
end

figure
subplot(3,1,1)
plot(t_elapsed,dx(7,:).'); hold on;
plot(t_elapsed,a1std,'r')
plot(t_elapsed,-a1std,'r')
ylabel('\delta \psi_1 [deg]')
grid on
if VO_updates == 1
    xline(t_elapsed(k_final))
end
legend('Error','\pm 3 \sigma')
subplot(3,1,2)
plot(t_elapsed,dx(8,:).'); hold on;
plot(t_elapsed,a2std,'r')
plot(t_elapsed,-a2std,'r')
ylabel('\delta \psi_2 [deg]')
grid on
if VO_updates == 1
    xline(t_elapsed(k_final))
end
subplot(3,1,3)
plot(t_elapsed,dx(9,:).'); hold on;
plot(t_elapsed,a3std,'r')
plot(t_elapsed,-a3std,'r')
ylabel('\delta \psi_3 [deg]')
grid on
if VO_updates == 1
    xline(t_elapsed(k_final))
end

figure
plot(t_elapsed,xE(7:10,:))
hold on
plot(t_elapsed,xE_true(7:10,:))
if VO_updates == 1
    xline(t_elapsed(k_final))
end

%% Plot residuals
if VO_updates == 1

    c = 1;
    for k = 1:length(r)

        if ~isempty(r{k})

            for j = 1:length(r{k})
                u_res(c) = r{k}{j}(1);
                v_res(c) = r{k}{j}(2);
                t_res(c) = k;
                c = c + 1;
            end

        end

    end

    figure
    subplot(2,1,1)
    plot(t_res,u_res,'r.')
    hold on
    plot(t_res,1*sigma_m*ones(1,length(t_res)),'b')
    plot(t_res,-1*sigma_m*ones(1,length(t_res)),'b')
    axis([0 t_res(end) -6*sigma_m 6*sigma_m])
    legend('Post-fit Residuals','\pm 1\sigma')
    ylabel('u [pix]')
    subplot(2,1,2)
    plot(t_res,v_res,'r.')
    hold on
    plot(t_res,1*sigma_m*ones(1,length(t_res)),'b')
    plot(t_res,-1*sigma_m*ones(1,length(t_res)),'b')
    axis([0 t_res(end) -6*sigma_m 6*sigma_m])
    xlabel('Time [s]')
    ylabel('v [pix]')

    % Recover features
    mu_final = mu(1:lists.nav_state.idx(end),k_final);
    R_final = R(1:lists.nav_state.Ridx(end),1:lists.nav_state.Ridx(end));
    %l = zeros(length(lists.features),3);
    %mU = zeros(length(lists.features),3);
    j = 1;
    for i = 1:length(lists.features)

        if ~isempty(lists.features{i})

            idx = lists.features{i}.idx;
            p = lists.features{i}.parent;
            pidx = lists.parent{p}.idx;

            x_i = mu_final(pidx(1:3));
            R_EN_i = q2C(mu_final(pidx(7:10)));
            c = mu_final(idx(1));
            m = lists.features{i}.m.';
            l(j,:) = x_i + R_EN_i*params.R_NC*m/c;

            mU(j,:) = ECEFtoENU(xE_true(1:3,1),l(j,:).',params,'vector').';

            j = j + 1;

        end

    end

    for k = 1:k_final

        xU(:,k) = ECEFtoENU(xE_true(1:3,1),xE(1:3,k),params,'vector');
        xU_true(:,k) = ECEFtoENU(xE_true(1:3,1),xE_true(1:3,k),params,'vector');

    end

    figure
    mcheck = mU(:,3)>-10000;
    plot3(mU(mcheck,1)./1000,mU(mcheck,2)./1000,mU(mcheck,3)./1000,'*')
    hold on
    plot3(xU_true(1,:)./1000,xU_true(2,:)./1000,xU_true(3,:)./1000)
    plot3(xU(1,:)./1000,xU(2,:)./1000,xU(3,:)./1000)
    xlabel('East [km]')
    ylabel('North [km]')
    zlabel('Up [km]')
    grid on
    axis equal
    legend('Landmarks','Truth','VO Solution')

else

    for k = 1:tspan

        xU(:,k) = ECEFtoENU(xE_true(1:3,1),xE(1:3,k),params,'vector');
        xU_true(:,k) = ECEFtoENU(xE_true(1:3,1),xE_true(1:3,k),params,'vector');

    end

    figure
    plot3(xU_true(1,:)./1000,xU_true(2,:)./1000,xU_true(3,:)./1000)
    hold on
    plot3(xU(1,:)./1000,xU(2,:)./1000,xU(3,:)./1000)
    xlabel('East [km]')
    ylabel('North [km]')
    zlabel('Up [km]')
    grid on
    axis equal
    legend('Truth','Integrated')

end