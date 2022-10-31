%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IPA SLAM
% Author: Matthew Givens
% Date: October 31, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%   x0 - Initial state
%   P0 - Initial covariance of states
%   z  - Measurements, formatted as cell array
%   params
%       -params.sizes - n (nav state size), m (error state size)
%       -params.noises - process and measurement noise parameters
%       -params.camera - camera calibration parameters
%       -params.tforms - relationships between frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IPA_SLAM(x0,P0,z,params)

%%%% Preliminaries %%%%
n = params.sizes.n;
m = params.sizes.m;
next_idx = 1;
Rnext_idx = 1;
new_feature = 0;

%%%% Pre-allocations %%%%
x = zeros(n,tspan); % Nav states
P = zeros(n,n,tspan); % Nav covariances
mu = zeros(n+3000,tspan); % Nav + landmark states
R = zeros(n-1+3000,n-1+3000); % All covariances (sliding window)
r = cell(1,maxI);

%%%% Initializations %%%%
x(:,1) = x0(:,1);
P(n,n,1) = P0(:,:,1);
mu(1:n,1) = x0(:,1);
R(1:n-1,1:n-1,1) = chol(eye(n-1)/P0(:,:,1));

% 'lists' is a nested structure that contains cell arrays for both features
% and parents as well as the current indices of the navigation state and
% the identifiers of all currently active landmarks. Features and parents
% are organized by the order of their instantiation and thus their
% identifiers correspond to their indicies in the cell array
lists.nav_state.idx = 1:n;
lists.nav_state.Ridx = 1:m;
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

for k = 2:tspan

    % IPA-SLAM Propagation
    dt = t(k) - t(k-1);
    u(1:3,1) = imu_array(k-1,2:4).';
    u(4:6,1) = imu_array(k-1,5:7).';

    mu(1:lists.nav_state.idx(end),k) = mu(1:lists.nav_state.idx(end),k-1);

    if new_feature == 0 % Propagate via Marginalization

        lists.mu_prop_window = lists.mu_update_window;
        lists.R_prop_window = lists.R_update_window;
        mu_in = mu(lists.mu_prop_window,k-1);
        R_in = R(lists.R_prop_window,lists.R_prop_window);

        [mu_out,R_out] = propagate(mu_in,R_in,u,dt,params);

    else % Propagate via Augmentation

        lists.mu_prop_window = lists.nav_state.idx;
        lists.R_prop_window = lists.nav_state.Ridx;

        lists.n_parent = lists.n_parent + 1;
        lists.parent{lists.n_parent}.idx = lists.nav_state.idx;
        lists.parent{lists.n_parent}.Ridx = lists.nav_state.Ridx;
        lists.parent_update = [lists.parent_update;lists.parent{lists.n_parent}.idx];

        lists.nav_state.idx = lists.nav_state.idx + n;
        lists.nav_state.Ridx = lists.nav_state.Ridx + m;

        mu_in = mu(lists.mu_prop_window,k-1);
        R_in = R(lists.R_prop_window,lists.R_prop_window);

        [mu_out,R_out] = augment(mu_in,R_in,u,dt,params);

        lists.mu_prop_window = [lists.mu_prop_window,lists.mu_prop_window + n];
        lists.mu_update_window = lists.mu_update_window(1):lists.mu_update_window(end) + n;
        lists.R_prop_window = [lists.R_prop_window,lists.R_prop_window+m];
        lists.R_update_window = lists.R_update_window(1):lists.R_update_window(end)+m;

        next_idx = lists.nav_state.idx(1);
        Rnext_idx = lists.nav_state.Ridx(1);
        new_feature = 0;

    end

    mu(lists.mu_prop_window,k) = mu_out;
    R(lists.R_prop_window,lists.R_prop_window) = R_out;

    x(:,k) = mu_out(end-n+1:end,1);
    S = eye(m)/R_out(end-m+1:end,end-params.m+1:end);
    P(:,:,k) = S*S.';

    % Update
    % Time check for update
    t_diff = abs(imageTimes(j) - t(k));
    if t_diff <= tol && VO_updates == 1 && j < 350

        if sum(z_ID(:,j)) ~= 0

            % Data Association
            lists = data_association(lists,z_ID(:,j),next_idx,Rnext_idx,params);
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

            [x_plus,R_plus,r{j},lists] = measurement_update(mu_in,R_in,z(:,j),...
                lists,params,I);
            mu(lists.mu_update_window,k) = x_plus;
            x(:,k) = x_plus(end-params.n+1:end,1);
            R(lists.R_update_window,lists.R_update_window) = R_plus;
            S = eye(params.m)/R_plus(end-params.m+1:end,end-params.m+1:end);
            P(:,:,k) = S*S.';

        end

        % Image counter
        j = j + 1;

    end

    percent = k*100/tspan;
    clc
    fprintf('Percent Complete: %3.2f\n',percent)

end

end