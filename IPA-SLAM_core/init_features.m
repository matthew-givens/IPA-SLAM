%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IPA SLAM Initialize Feature
% Author: Matthew Givens
% Date: August 24, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_init,R_init,lists] = init_features(x_minus,R_minus,init_list,z,params,lists,I)

nNew = size(init_list,1);

x_init = [zeros(nNew,1);x_minus];

Ptemp = zeros(nNew+params.m);

p = x_minus(1:3,1);
c0 = 1/(norm(p)-params.r_E);
sigma_c0 = c0*100;

idxNew = 1;
for i = 1:nNew

    % Project and rotate measurement to earth-fixed frame
    % undistort to use pinhole model
    hC_d = [z{i}.';1];
    hC_u(:,i) = undistortPoints(hC_d(1:2,1).',params.intrinsics).';
    hC = params.Kinv*[hC_u(1:2,i);1];
    lists.features{init_list(i)}.m = (hC./norm(hC)).';

    % Compute first-estimate Jacobian
    

    % Stack new ID states
    x_init(idxNew,1) = c0;

    Ptemp(idxNew,idxNew) = sigma_c0^2;
    idxNew = idxNew + 1;

end

if params.debug_image == 1
    currI = undistortImage(I,params.intrinsics,'OutputView','full');
    imshow(I); hold on;
    points = plot(hC_u(1,:),hC_u(2,:),'b+');

    delete(points)
end

Rtemp = [inv(chol(Ptemp(1:nNew,1:nNew))),zeros(nNew,params.m);...
    zeros(params.m,nNew),R_minus(end-params.m+1:end,end-params.m+1:end)];
R_init = Rtemp;

end