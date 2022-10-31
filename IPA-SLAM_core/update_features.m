%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IPA SLAM Update Features
% Author: Matthew Givens
% Date: August 24, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_plus,R_plus,r,points] = update_features(x_minus,R_minus,...
    update_list,z,params,lists,I)

nUpdate = size(update_list,1);
break_flag = 0;
points = [];
x_iter = x_minus;
S = eye(size(R_minus,1))/R_minus;
for q = 1:params.num_iter

% Pre-allocate Jacobian
H = zeros(2*nUpdate,size(R_minus,2));
db_iter = zeros(size(R_minus,1),1);
res = zeros(2*nUpdate,1);

    for i = 1:nUpdate

        x_nav = x_iter(end-9:end,1);

        R_EN = q2C(x_nav(7:10,1));
        R_EC = R_EN*params.R_NC;

        current_feature = update_list(i,:);
        current_parent = lists.features{current_feature}.parent;

        idx = lists.features{current_feature}.idx - lists.mu_offset;
        R_idx = lists.features{current_feature}.Ridx - lists.R_offset;

        parent_idx = lists.parent{current_parent}.idx - lists.mu_offset;
        R_parent_idx = lists.parent{current_parent}.Ridx - lists.R_offset;

        Hxyz = zeros(1,size(R_minus,2));

        x_j = x_iter(parent_idx(1,1:3),1);
        R_EN_j = q2C(x_iter(parent_idx(1,7:10),1));
        c = x_iter(idx(1),1);

        % Update
        m = lists.features{current_feature}.m.';
        pc_k = x_nav(1:3,1) + R_EN*params.pN_C;
        pc_j = x_j(1:3,1) + R_EN_j*params.pN_C;
        delta1 = pc_j-pc_k;
        delta3 = params.pN_C+1/c*params.R_NC*m;
        delta2 = x_j(1:3,1) - x_nav(1:3,1) + R_EN_j*delta3;
        hE = c*delta1 + R_EN_j*params.R_NC*m;
        hC = R_EC.'*hE;
        hC = hC/hC(3);
        z_hat(:,i) = params.K*hC;
        zu(:,i) = undistortPoints(z{i},params.intrinsics).';

        % Compute Jacobians
        dhdx = -c*R_EN.'*[eye(3),zeros(3),X(delta2)];
        dhdc = R_EN.'*delta1;
        dhdmx = c*R_EN.'*[eye(3),zeros(3),X(R_EN_j*delta3)];

        % Populate H matrix
        Hxyz(1:3,end-params.m+1:end) = dhdx;
        Hxyz(1:3,R_idx) = dhdc;
        Hxyz(1:3,R_parent_idx) = dhdmx;
        Huv = params.K*1/hC(3)*(eye(3)-hC/hC(3)*[0,0,1])*params.R_NC.'*Hxyz;
        
        % Underweighting
        HS = H(2*i-1:2*i,:)*S;
        V = params.V + params.beta*(HS*HS.');
        R_C = eye(2)/V;

        H(2*i-1:2*i,:) = R_C*Huv(1:2,:);
        r{i} = zu(:,i) - z_hat(1:2,i);
        r_tilde = R_C*r{i};
        res(2*i-1:2*i,1) = r_tilde;

        % Mahalanobis Distance Test for Outlier Rejection
        if params.gating == 1
            U = H(2*i-1:2*i,:)/R_minus;
            Sinv = eye(2)/(U*U.' + params.V);
            gamma(i) = sqrt(r{i}.'*Sinv*r{i});
            thresh = chi2inv(0.95,size(params.m+3*nUpdate,2));
            if gamma(i) > thresh
                H(2*i-1:2*i,:) = zeros(2,size(R_minus,1));
                res(2*i-1:2*i,1) = [0;0];
                r{i} = NaN(2,1);
                break_flag = 1;
                break
            end
        end

    end

    if params.debug_image == 1
        currI = undistortImage(I,params.intrinsics,'OutputView','full');
        imshow(I); hold on;
        points1 = plot(z_hat(1,:),z_hat(2,:),'ro');
        points2 = plot(zu(1,:),zu(2,:),'b+');

        delete(points1);delete(points2);
    end

    R_pre = [R_minus,db_iter;
             H,res];
    [~,R_post] = qr(R_pre);

    R_iter = R_post(1:end-2*nUpdate,1:end-1);
    db_iter = R_post(1:end-2*nUpdate,end);

    dx(:,q) = R_iter\db_iter;

    % dx is indexed in terms of the reduced (error) state which makes
    % it problematic when updating the full state. The error state will
    % always have dimension 1+numParent less than the full state. The
    % challenge then is to map dx to dx_full and then reconcile the error
    % state to the full state

    parent_list_update = lists.parent_update - lists.mu_offset;
    x_iter = update_state(x_iter,dx(:,q),parent_list_update(parent_list_update(:,1)>0,:),params);

end

x_plus = x_iter;
R_plus = R_iter;

%%%%%% POSTFIT RESIDUALS %%%%%%
if params.postfit == 1 && break_flag ~= 1

    x_nav = x_plus(end-9:end,1);
    R_EN = q2C(x_nav(7:10,1));
    R_EC = R_EN*params.R_NC;

    for i = 1:nUpdate

        current_feature = update_list(i,:);
        current_parent = lists.features{current_feature}.parent;
        idx = lists.features{current_feature}.idx - lists.mu_offset;
        parent_idx = lists.parent{current_parent}.idx - lists.mu_offset;

        x_j = x_iter(parent_idx(1,1:3),1);
        R_EN_j = q2C(x_iter(parent_idx(1,7:10),1));
        c = x_iter(idx(1),1);

        m = lists.features{current_feature}.m.';
        pc_k = x_nav(1:3,1) + R_EN*params.pN_C;
        pc_j = x_j(1:3,1) + R_EN_j*params.pN_C;
        delta1 = pc_j-pc_k;
        hE = c*delta1 + R_EN_j*params.R_NC*m;
        hC = R_EC.'*hE;
        hC = hC/hC(3);
        z_hat(:,i) = params.K*hC;
        zu(:,i) = undistortPoints(z{i},params.intrinsics).';

        r{i} = zu(:,i) - z_hat(1:2,i);

    end

    points.z_hat = z_hat;
    points.z = zu;

    if params.debug_image == 1
        currI = undistortImage(I,params.intrinsics,'OutputView','full');
        imshow(I); hold on;
        points1 = plot(z_hat(1,:),z_hat(2,:),'ro','MarkerSize',10);
        points2 = plot(zu(1,:),zu(2,:),'b+','MarkerSize',10);
        
        delete(points1);delete(points2);
    end

end

end