% Acceleration: a> 0.1 m/s2 
% Deceleration: a< -0.1 m/s2 
% Cruising: -0.1 m/s2 ≤ a ≤0.1 m/s2

% T, V, P, D
function cur_Traj = identify_DrMo(cur_Traj, threshold, t_dim, v_dim, p_dim)

format long

[m,~]=size(cur_Traj);
driving_mode=zeros(m,1); % [cruise, dec, idle, acc, stop] ---> [0, 1, 2, 3, 4]

T=cur_Traj(:,t_dim);
V=cur_Traj(:,v_dim); % km/h
P=cur_Traj(:,p_dim); % SA/MA
for i = 3:m
    if P(i) == 1
        dV = V(i) - V(i-1);
        a = dV / 3.6 / T(i); % avg acceleration: m/s2
        if dV>120 || dV<-120 || V(i)>150
            driving_mode(i) = 4;
        elseif a > threshold
            driving_mode(i) = 3;
        elseif a < -threshold
            driving_mode(i) = 1;
        else
            driving_mode(i) = 0;
        end
    elseif P(i) == 2 
        driving_mode(i) = 2;
    elseif P(i) == 3 
        driving_mode(i) = 4;
    end
end

% the end point of a cur_Trajectory
if P(m) == 1
    if V(m) > 5 
        driving_mode(m) = 0;
    elseif V(m) <= 5 
        driving_mode(m) = 2;
    end
elseif P(m) == 2 
    driving_mode(m) = 2;
elseif P(m) == 3 
    driving_mode(m) = 4;
end

if size(cur_Traj,2) <= p_dim+1
    cur_Traj = [cur_Traj driving_mode];
else
    cur_Traj(:, p_dim+2) = driving_mode;
end

end

