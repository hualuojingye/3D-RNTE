% Acceleration: a> 0.1 m/s2 
% Deceleration: a< -0.1 m/s2 
% Cruising: -0.1 m/s2 ≤ a ≤0.1 m/s2


% T, V, P, D
function Identify_DrMo(threshold, base_path, param_idxs)

format long

t_idx = param_idxs.t_idx;
v_idx = param_idxs.v_idx;
p_idx = param_idxs.p_idx;
d_idx = param_idxs.d_idx;

filename_list = dir(base_path+"*.txt");
for filename_idx = 1:length(filename_list)
    Traj = load(base_path+filename_list(filename_idx).name);
    [m,~]=size(Traj);
    driving_mode=zeros(m,1); % [cruise, dec, idle, acc, stop] ---> [0, 1, 2, 3, 4]
    
    T=Traj(:,t_idx);
    V=Traj(:,v_idx); % km/h
    P=Traj(:,p_idx); % SA/MA
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

    % the end point of a trajectory
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


    if size(Traj,2) < d_idx
        Traj = [Traj driving_mode];
    else
        Traj(:, d_idx) = driving_mode;
    end
    writematrix(Traj, base_path+filename_list(filename_idx).name);
end
