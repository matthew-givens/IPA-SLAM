function output = ECItoECEF(input,t,params)

theta = params.theta_f + params.omega_E(3)*t;

R = [cos(theta),-sin(theta),0;  
    sin(theta),cos(theta),0;
    0,0,1].';

if size(input,2) == 3

    output = R;
    
else
    
    if size(input,1) == 3
    
        output = R*input;
        
    elseif size(input,1) == 6
        
        r = R*input(1:3,1);
        v = R*input(4:6,1) - cross(params.omega_e,R*r);
        output = [r;v];
        
    end

end

end