% Jan 2022 - new quaternion functions (JPL/Breckenridge convention: [x;y;z;w])
% qX generates the skew symmetric matrix qx

function qx = qX(q)

qx = [q(4,1)*eye(3) - X(q(1:3,1)),q(1:3,1);
      -q(1:3,1).',q(4,1)];

end