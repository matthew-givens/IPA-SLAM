function output = ECEFtoNED(input,params)

r = input(1:3,1);

yhat = cross(params.omega_e,r)/norm(cross(params.omega_e,r));
zhat = -r/norm(r);
xhat = cross(yhat,zhat);

R = [xhat yhat zhat].';

if size(input,2) == 3
    
    output = R;
    
else
    
    output = R*input;
   
end

end