function output = ECEFtoENU(origin,p,params,outflag)

%%% Inputs %%%
% Use column vectors
% 'origin' is the position on Earth's surface in ECEF that is to be the
% origin of the new ENU frame
% 'p' is another position vector of interest in the ECEF frame
% 'params' must include params.omega_E, the rotation rate in rad/s of the 
% planet, e.g. params.omega_E = 2*pi/T where T is the rotation period in s
% 'outflag' can be either 'vector' or 'matrix', the latter of which only
% depends on the 'origin' input

%%% Outputs %%%
% Depending on 'outflag', the output is either a 'vector' in ENU pointing
% from the 'origin' to 'p' in ENU or it's a 'matrix' that represents the 
% rotation from ECEF to ECI 

r = origin(1:3,1);

xhat = cross(params.omega_E,r)/norm(cross(params.omega_E,r));
zhat = r/norm(r);
yhat = cross(zhat,xhat);

R = [xhat yhat zhat].';

position_ENU = R*p;

switch outflag
    
    case 'matrix'
        
        output = R;
        
    case 'vector'
        
        output = position_ENU - R*origin;
        
end

end