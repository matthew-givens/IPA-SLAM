% Jan 2022 - new quaternion functions (JPL/Breckenridge convention: [x;y;z;w])
% qinv outputs the congugate quaternion of q

function qc = qconj(q)

qc = [-q(1:3,1);q(4)];

end