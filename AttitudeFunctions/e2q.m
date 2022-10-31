% Jan 2022 - new quaternion functions (JPL convention: [x;y;z;w])
% q2C outputs the quaternion q corresponding to euler axis/angle e,theta

function q = e2q(e,theta)

q = [e*sin(theta/2);
     cos(theta/2)];

end