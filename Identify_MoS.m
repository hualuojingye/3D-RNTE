% Identification of MA (1) SA with engine-on (2) SA with engine-off (3)
% T, L, V, P
function Identify_MoS(threshold, base_path, param_idxs)

t_idx = param_idxs.t_idx;
l_idx = param_idxs.l_idx;
v_idx = param_idxs.v_idx;
p_idx = param_idxs.p_idx;

format long

filename_list = dir(base_path+"*.txt");
for filename_idx = 1:length(filename_list)
    Traj = load(base_path+filename_list(filename_idx).name);
    [m,~]=size(Traj);
    Cat=zeros(m,1); % Cat: 1:MA; 2:SA with engine on; 3:SA with engine-off
    
    v23=Traj(:,v_idx); % speed: km/h
    l23=Traj(:,l_idx);  % distance: m
    t23=Traj(:,t_idx);  % time interval: s
    
    % Define two variable-length matrices SA[] (to store SA of uncertain length);
    SA=[];
    i=2;
    while i<=m
        s=l23(i);
        tsa=0;
        if s>=threshold  % MA
            if (v23(i)-v23(i-1))<=120
                Cat(i)=1;
            else
                Cat(i)=3;
            end
            i=i+1;
        end
        
        if s<threshold  % SA
            SA=[SA,s];
            for j=i+1:m
                s=l23(j);
                
                if s<threshold
                    SA=[SA,s];
                else
                    break;
                end
            end % Obtain the SA matrices with consecutive distances less than 2 meters, with indices from i to j-1
            
            SA=[];
            % Check if the stay time is greater than 3 minutes
            for k=i:j-1
                tsa=tsa+t23(k);
            end
            if tsa>180  %SA with engine off cat=3
                for p= i:j-1
                    Cat(p)=3;
                end
            else
                for q= i:j-1
                    Cat(q)=2;
                end
            end
            i=j;
        end
    end


    if size(Traj, 2) < p_idx
        Traj = [Traj Cat];
    else
        Traj(:, p_idx) = Cat;
    end
    writematrix(Traj, base_path+filename_list(filename_idx).name);
end
