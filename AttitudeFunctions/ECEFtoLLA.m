function out = ECEFtoLLA(x,r_eq)

% Altitude
r = norm(x);

% Latitude
p_vec = [x(1:2,1);0];
p = norm(p_vec);
phi = acos(dot(x,p_vec)/(r*p));

% Longitude
xdir = [1;0;0];
lambda = -acos(dot(xdir,p_vec)/p);

out = [phi;lambda;r-r_eq];

end