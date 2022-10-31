function x_plus = update_state(x_minus,dx,parent_list,params)
% First entry will always be a feature
% Last entry will always be attitude state
x_size = size(x_minus,1);
parent_list = [parent_list;zeros(1,10)];

% dsize will be the number of parents + 1
x_plus = zeros(x_size,1);
i = 1;
j = 1;
k = 1;
for q = 1:x_size - params.n

    if i == parent_list(k,7)

        % Compute quaternion update
        dq = qinv([dx(j:j+2,1)./2;sqrt(1-0.25*dx(j:j+2,1).'*dx(j:j+2,1))]);
        x_plus(i:i+3,1) = qXp(qinv(dq),x_minus(i:i+3,1));

        x_plus(i:i+3,1) = x_plus(i:i+3,1)/norm(x_plus(i:i+3,1));
        i = i + 4;
        j = j + 3;
        k = k + 1;

        if k == size(parent_list,1)
            break
        end
        
    else

        x_plus(i) = x_minus(i) + dx(j);
        i = i + 1;
        j = j + 1;

    end

end

x_plus(end-9:end-4) = x_minus(end-9:end-4) + dx(end-8:end-3);
%dq = [0;0;0;1];
dq = qinv([dx(end-2:end,1)./2;sqrt(1-0.25*dx(end-2:end,1).'*dx(end-2:end,1))]);
x_plus(end-3:end,1) = qXp(qinv(dq),x_minus(end-3:end,1));
x_plus(end-3:end,1) = x_plus(end-3:end,1)/norm(x_plus(end-3:end,1));

end