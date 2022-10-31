% Jan 2022 - new quaternion functions (JPL convention: [x;y;z;w])
% C2e outputs the euler axis/angle corresponding to rotation matrix C

function [e,theta] = C2e(C)

theta = acos((trace(C)-1)/2);

if theta == 0
    
    e = [0;0;1];
    
elseif theta == pi
        
        % Eq 2.115
        eeT = (C + eye(3))/2;
        
        [val, i] = max([norm(eeT(:,1)),norm(eeT(:,2)),norm(eeT(:,3))]);
        
        if val == 0
            
            e = [0;0;0];
            
        else
            
            e = eeT(:,i)/val;
            
        end
        
else
    
    % Eq 2.114
    e = 1/(2*sin(theta))*[C(2,3)-C(3,2);C(3,1)-C(1,3);C(1,2)-C(2,1)];
        
end

end

