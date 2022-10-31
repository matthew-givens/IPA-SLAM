% Jan 2022 - new quaternion functions (JPL convention: [x;y;z;w])
% q2C outputs the quaternion q corresponding to rotation matrix C

function q = C2q(C)

[e,theta] = C2e(C);

q = e2q(e,theta);

end