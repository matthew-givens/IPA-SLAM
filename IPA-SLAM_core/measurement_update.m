%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IPA SLAM Update
% Author: Matthew Givens
% Date: August 24, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_plus,R_plus,r,lists,points] = measurement_update(x_minus,R_minus,z,lists,...
    params,I)

n = params.sizes.n;
m = params.sizes.m;
r = [];
points = [];

update_list = setdiff(lists.active,lists.new);
init_list = lists.new;

%%%%%% UPDATING PERSISTENT FEATURES %%%%%%
if ~isempty(update_list)

    [x_plus,R_plus,r,points] = update_features(x_minus,R_minus,...
        update_list,z(update_list),params,lists,I);

else

    x_plus = x_minus;
    R_plus = R_minus;

end

%%%%%% INITIALIZING NEW FEATURES %%%%%%
% For each feature to be initialized, estimated position and attitude are the same
if ~isempty(init_list)

    [x_init,R_init,lists] = init_features(x_plus(end-params.n+1:end,1),...
        R_plus(end-params.m+1:end,end-params.m+1:end),init_list,z(init_list),params,lists,I);

    nNew = size(init_list,1);
    for j = 1:nNew
        lists.features{lists.new(j)}.status = 1;
    end
    lists.nav_state.idx = lists.nav_state.idx + nNew;
    lists.nav_state.Ridx = lists.nav_state.Ridx + nNew;

    % Update matrix
    R_new = zeros(size(R_plus,1)+nNew);
    R_new(end-m-nNew+1:end,end-m-nNew+1:end) = R_init;
    R_new(1:end-m-nNew,1:end-m-nNew) = R_plus(1:end-m,1:end-m);
    R_new(1:end-m-nNew,end-m+1:end) = R_plus(1:end-m,end-m+1:end);
    R_plus = R_new;

    x_new = zeros(size(x_plus,1)+nNew,1);
    x_new(1:size(x_plus,1),1) = x_plus;
    x_new(end-n-nNew+1:end,1) = x_init;
    x_plus = x_new;

    lists.mu_update_window = lists.mu_update_window(1):lists.mu_update_window(end)+nNew;
    lists.R_update_window = lists.R_update_window(1):lists.R_update_window(end)+nNew;
    
end

end