function x_out = convert_orbit_param(inputType,x_in,t0,outputType,mu)

% This function is designed to convert the vector x_in that is of type
% inputType to a vector of type outputType. The standard gravitational
% parameter mu for the relevant body must be provided. The function should
% work for all non-degenerate eccentricities. For eccentricities greater
% than one, simply replace M or M_0 with N or N_0 but use the same
% inputType argument as for eccentricities less than one. If the invariant
% set is used (M0), the initial time can be passed in using the t0 input,
% otherwise this input is unused.

% Classical element sets should be formatted as
% x_in = [a;e;inc,RA,w,M/f/M0]

% Angle inputs should be in radians, angle outputs are in radians
% Possible types for inputType and outputFlag are:
% - 'states' - 6x1 position/velocity vector in ECI frame
% - 'classical_elements_M' - classical elements w/mean anomaly in x_in(6)
% - 'classical_elements_f' - classical elements w/true anomaly in x_in(6)
% - 'classical_elements_M0' - classical elements with invariant initial
%   mean anomaly in x_in(6)
% More element sets to be added in the future

switch inputType
    
    case 'classical_elements_M'
        
        switch outputType
            
            case 'states'
                
                a = x_in(1);
                e = x_in(2);
                inc = x_in(3);
                RA = x_in(4);
                w = x_in(5);
                
                % Eccentricity check
                if e > 1
                    
                    N = x_in(6);
                    
                    % Solving for hyperbolic anomaly
                    H_guess = N;
                    H = newton_solve_H(N,e,H_guess);
                    f = 2*atan(sqrt((e+1)/(e-1))*tanh(H/2));
                    
                elseif e < 1
                    
                    M = x_in(6);
                    
                    % Solving for eccentric anomaly
                    E_guess = M;
                    E = newton_solve_E(M,e,E_guess);
                    f = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
                    
                end
                
                % Universal equations
                p = a*(1-e^2);
                r_norm = p/(1+e*cos(f));
                v_norm = sqrt(mu*(2/r_norm-1/a));
                h = sqrt(mu*p);
                theta = w + f;
                
                % Vectors of position and velocity
                r = r_norm*[cos(RA)*cos(theta)-sin(RA)*sin(theta)*cos(inc);
                    sin(RA)*cos(theta)+cos(RA)*sin(theta)*cos(inc);
                    sin(theta)*sin(inc)];
                v = -mu/h*[cos(RA)*(sin(theta)+e*sin(w))+sin(RA)*(cos(theta)+e*cos(w))*cos(inc);
                    sin(RA)*(sin(theta)+e*sin(w))-cos(RA)*(cos(theta)+e*cos(w))*cos(inc);
                    -(cos(theta)+e*cos(w))*sin(inc)];
                
                x_out = [r;v];
                
            case 'classical_elements_f'
                
                % Solving for eccentric anomaly, true anomaly
                e = x_in(2);
                M = x_in(6);
                E_guess = M;
                E = newton_solve_E(M,e,E_guess);
                f = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
                
                x_out = [x_in(1:5,1);f];
                
            case 'qnonsingular'
                
                a = x_in(1);
                e = x_in(2);
                inc = x_in(3);
                RA = x_in(4);
                w = x_in(5);
                M = x_in(6);
                
                x = convert_orbit_param('classical_elements_M',x_in,0,'states',mu);
                r = x(1:3,1);
                v = x(4:6,1);
                h = cross(r,v);
                i_e = cross(v,h)/mu-r/norm(r);
                e_x = i_e(1);
                e_y = i_e(2);
                
                % Solving for eccentric anomaly, true anomaly
                E_guess = M;
                E = newton_solve_E(M,e,E_guess);
                f = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
                u = w + f;
                
                x_out = [a;e_x;e_y;inc;RA;u];
                
            otherwise
                
                error('Not a valid entry for [outputFlag]')
                
        end
        
    case 'classical_elements_f'
        
        a = x_in(1);
        e = x_in(2);
        inc = x_in(3);
        RA = x_in(4);
        w = x_in(5);
        f = x_in(6);
                        
        switch outputType
            
            case 'states'

                % Universal equations
                p = a*(1-e^2);
                r_norm = p/(1+e*cos(f));
                v_norm = sqrt(mu*(2/r_norm-1/a));
                h = sqrt(mu*p);
                theta = w + f;
                
                % Vectors of position and velocity
                r = r_norm*[cos(RA)*cos(theta)-sin(RA)*sin(theta)*cos(inc);
                    sin(RA)*cos(theta)+cos(RA)*sin(theta)*cos(inc);
                    sin(theta)*sin(inc)];
                v = -mu/h*[cos(RA)*(sin(theta)+e*sin(w))+sin(RA)*(cos(theta)+e*cos(w))*cos(inc);
                    sin(RA)*(sin(theta)+e*sin(w))-cos(RA)*(cos(theta)+e*cos(w))*cos(inc);
                    -(cos(theta)+e*cos(w))*sin(inc)];
                
                x_out = [r;v];
                
            case 'qnonsingular'
                
                x = convert_orbit_param('classical_elements_f',x_in,0,'states',mu);
                r = x(1:3,1);
                v = x(4:6,1);
                h = cross(r,v);
                i_e = cross(v,h)/mu-r/norm(r);
                e_x = i_e(1);
                e_y = i_e(2);
                
                u = w + f;
                x_out = [a;e_x;e_y;inc;RA;u];
                
            otherwise
                
                error('Not a valid entry for [outputFlag]')
                
        end
        
    case 'classical_elements_M0'
        
        switch outputType
            
            case 'states'
                
                a = x_in(1);
                e = x_in(2);
                inc = x_in(3);
                RA = x_in(4);
                w = x_in(5);
                
                % Eccentricity check
                if e > 1
                    
                    % Change in N since t_0
                    N_0 = x_in(6);
                    n = sqrt(mu/(-a)^3);
                    N = N_0 + n*t0;
                    
                    % Solving for hyperbolic anomaly
                    H_guess = N;
                    H = newton_solve_H(N,e,H_guess);
                    f = 2*atan(sqrt((e+1)/(e-1))*tanh(H/2));
                    
                elseif e < 1
                    
                    % Change in M since t_0
                    M_0 = x_in(6);
                    n = sqrt(mu/a^3);
                    M = M_0 + n*t0;
                    
                    % Solving for eccentric anomaly
                    E_guess = M;
                    E = newton_solve_E(M,e,E_guess);
                    f = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
                    
                end
                
                % Universal equations
                p = a*(1-e^2);
                r_norm = p/(1+e*cos(f));
                v_norm = sqrt(mu*(2/r_norm-1/a));
                h = sqrt(mu*p);
                theta = w + f;
                
                % Vectors of position and velocity
                r = r_norm*[cos(RA)*cos(theta)-sin(RA)*sin(theta)*cos(inc);
                    sin(RA)*cos(theta)+cos(RA)*sin(theta)*cos(inc);
                    sin(theta)*sin(inc)];
                v = -mu/h*[cos(RA)*(sin(theta)+e*sin(w))+sin(RA)*(cos(theta)+e*cos(w))*cos(inc);
                    sin(RA)*(sin(theta)+e*sin(w))-cos(RA)*(cos(theta)+e*cos(w))*cos(inc);
                    -(cos(theta)+e*cos(w))*sin(inc)];
                
                x_out = [r;v];
                
            otherwise
                
                error('Not a valid entry for [outputFlag]')
                
        end
        
        
    case 'states'
        
        switch outputType
            
            case 'classical_elements_M'
                
                r = x_in(1:3);
                v = x_in(4:6);
                h = cross(r,v);
                
                % Perifocal frame
                i_e = cross(v,h)/mu-r/norm(r);
                i_h = h/norm(h);
                i_p = cross(i_h,i_e);
                i_r = r/norm(r);
                
                % Semi-major axis
                a = -mu/2*(norm(v)^2/2-mu/norm(r))^(-1);
                
                % Eccentricity
                e = norm(i_e);
                
                % Inclination (inc)
                inc = acos(i_h(3));
                
                % Right Ascension of the Ascending Node (RAAN)
                RA = atan2(i_h(1),-i_h(2));
                
                % Argument of periapsis (arg)
                w = atan2(i_e(3),i_p(3));
                
                if e > 1
                    
                    % Hyperbolic anomaly
                    n = sqrt(mu/(-a)^3);
                    p = a*(1-e^2);
                    f = acos((1/e)*((p/norm(r))-1));
                    H = 2*atanh(tan(f/2)/sqrt((e+1)/(e-1)));
                    N = e*sinh(H) - H;
                    
                    % Orbital elements
                    x_out = [a,e,inc,RA,w,N];
                    
                elseif e < 1
                    
                    % Mean anomaly
                    n = sqrt(mu/a^3);
                    sigma = dot(r,v)/sqrt(mu);
                    E = atan2(sigma/sqrt(a),1-norm(r)/a);
                    M = E - e*sin(E);
                    
                    % Orbital elements
                    x_out = [a;e;inc;RA;w;M];
                    
                end
                
            case 'classical_elements_f'
                
                r = x_in(1:3);
                v = x_in(4:6);
                h = cross(r,v);
                
                % Perifocal frame
                i_e = cross(v,h)/mu-r/norm(r);
                i_h = h/norm(h);
                i_p = cross(i_h,i_e);
                i_r = r/norm(r);
                PN = [i_e,i_p,i_h].';
                
                % Semi-major axis
                a = 1/(2/norm(r)-norm(v)^2/mu);
                
                % Eccentricity
                e = norm(i_e);
                
                % Inclination (inc)
                inc = acos(PN(3,3));
                
                % Right Ascension of the Ascending Node (RAAN)
                RA = atan2(PN(3,1),-PN(3,2));
                
                % Argument of periapsis (arg)
                w = atan2(PN(1,3),PN(2,3));
                
                if e > 1
                    
                    % Hyperbolic anomaly
                    n = sqrt(mu/(-a)^3);
                    p = a*(1-e^2);
                    f = acos((1/e)*((p/norm(r))-1));
                    
                    % Orbital elements
                    x_out = [a,e,inc,RA,w,f];
                    
                elseif e < 1
                    
                    % True anomaly
                    f = atan2(dot(cross(i_e,i_r),i_h),dot(i_e,i_r));
                    
                    % Orbital elements
                    x_out = [a;e;inc;RA;w;f];
                    
                end
                
            case 'classical_elements_M0'
                
                r = x_in(1:3);
                v = x_in(4:6);
                h = cross(r,v);
                
                % Perifocal frame
                i_e = cross(v,h)/mu-r/norm(r);
                i_h = h/norm(h);
                i_p = cross(i_h,i_e);
                i_r = r/norm(r);
                
                % Semi-major axis
                a = -mu/2*(norm(v)^2/2-mu/norm(r))^(-1);
                
                % Eccentricity
                e = norm(i_e);
                
                % Inclination (inc)
                inc = acos(i_h(3));
                
                % Right Ascension of the Ascending Node (RAAN)
                RA = atan2(i_h(1),-i_h(2));
                
                % Argument of periapsis (arg)
                w = atan2(i_e(3),i_p(3));
                
                if e > 1
                    
                    % Hyperbolic anomaly
                    n = sqrt(mu/(-a)^3);
                    p = a*(1-e^2);
                    f = acos((1/e)*((p/norm(r))-1));
                    H = 2*atanh(tan(f/2)/sqrt((e+1)/(e-1)));
                    N = e*sinh(H) - H;
                    N_0 = N - n*t0;
                    
                    % Orbital elements
                    x_out = [a,e,inc,RA,w,N_0];
                    
                elseif e < 1
                    
                    % Mean anomaly
                    n = sqrt(mu/a^3);
                    sigma = dot(r,v)/sqrt(mu);
                    E = atan2(sigma/sqrt(a),1-norm(r)/a);
                    M = E - e*sin(E);
                    M_0 = mod(M - n*t0,2*pi);
                    
                    % Orbital elements
                    x_out = [a;e;inc;RA;w;M_0];
                    
                end
                
            case 'qnonsingular'
                
                oe = convert_orbit_param('states',x_in,0,'classical_elements_f',mu);
                x_out = convert_orbit_param('classical_elements_f',oe,0,'qnonsingular',mu);
                                               
            otherwise
                
                error('Not a valid entry for [outputFlag]')
                
        end
        
    case 'qnonsingular'
        
        a = x_in(1);
        e_x = x_in(2);
        e_y = x_in(3);
        inc = x_in(4);
        RA = x_in(5);
        u = x_in(6);
        e = sqrt(e_x^2+e_y^2);
        
        w = acos(e_x/e);
        f = u - w;
        
        switch outputType
            
            case 'states'

                oe = [a;e;inc;RA;w;f];
                
                x_out = convert_orbit_param('classical_elements_f',oe,0,'states',mu);
                                
            case 'classical_elements_f'
                
                x_out = [a;e;inc;RA;w;f];
                
        end
        
    otherwise
        
        error('Not a valid entry for [inputType]')
        
end

end

function E = newton_solve_E(M,e,E_guess)
% Solves for the eccentric anomaly given a mean anomaly, eccentricity, and
% an initial guess.

% Initialization
r = 1;
k = 1;
E = E_guess;

tol = 1e-10; % tolerance

while r > tol
    
    ratio = (E-e*sin(E)-M)/(1-e*cos(E));
    E_new = E - ratio;
    
    r = abs(ratio);
    E = E_new;
    
    k = k + 1; % Counter
    
end

end

function H = newton_solve_H(N,e,H_guess)
% Solves for the hyperbolic eccentric anomaly given a mean anomaly,
% eccentricity, and an initial guess.

% Initialization
r = 1;
k = 1;
H = H_guess;

tol = 1e-10; % tolerance

while r > tol
    
    ratio = (e*sinh(H)-H-N)/(e*cosh(H)-1);
    H_new = H - ratio;
    
    r = abs(ratio);
    H = H_new;
    
    k = k + 1; % Counter
    
end

end