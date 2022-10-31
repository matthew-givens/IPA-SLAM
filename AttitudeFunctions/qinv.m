% Jan 2022 - new quaternion functions (JPL convention: [x;y;z;w])
% qinv outputs the inverse quaternion of q

function qi = qinv(q)

qi = qconj(q)/norm(q)^2;

end