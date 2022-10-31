% Jan 2022 - new quaternion functions (JPL/Breckenridge convention: [x;y;z;w])
% ijk = 1, i*i = j*j = k*k = -1
% q2C outputs the rotation matrix C corresponding to quaternion q

function C = q2C(q)

% Same equation, different forms
% C = eye(3) - 2*q(4,1)*X(q(1:3,1)) + 2*X(q(1:3,1))^2;
C = (2*q(4)^2-1)*eye(3) - 2*q(4)*X(q(1:3,1)) + 2*q(1:3,1)*q(1:3,1).';

end