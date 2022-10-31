function v_lever = lever_arm_effect(dtheta,dthetak,r,dt)

v_cent = -1/dt*[(dtheta(2)^2+dtheta(3)^2)*r(1);
                 (dtheta(1)^2+dtheta(3)^2)*r(2);
                 (dtheta(1)^2+dtheta(2)^2)*r(3)];
v_eul = zeros(3,1);%cross((dthetak-dtheta)/dt,r);

v_lever = v_cent + v_eul;

end