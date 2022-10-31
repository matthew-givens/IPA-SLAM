% Jan 2022 - new quaternion functions (JPL/Breckenridge convention: [x;y;z;w])
% qXp multiplies q (otimes) p

function qp = qXp(q,p)

qp = qX(q)*p;

end