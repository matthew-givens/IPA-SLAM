function [q_new,A] = updateQuaternion(da,q_old)

% Update quaternion based on output of coning integration (dtheta = da)
% Scalar component is last

gamma = norm(da)/2;

if gamma >= 1e-5

    Psi = sin(gamma)/2/gamma*da;

else

    Psi = da./2;

end

A = [cos(gamma)*eye(3)-X(Psi),Psi;
         -Psi.',cos(gamma)];
q_new = A*q_old;

q_new = q_new/norm(q_new);

end