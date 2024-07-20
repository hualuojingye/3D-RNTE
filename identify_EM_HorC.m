function cur_traj = identify_EM_HorC(cur_traj, l_dim, p_dim, HorC_dim)

format long;

[m,~]=size(cur_traj);

EM_HorC=zeros(m,1); % EM_HorC: 0:Hot phase; 1:Cold phase;

L=cur_traj(:,l_dim);  % distance: m
Pattern=cur_traj(:,p_dim);

L_all=sum(L(:,1)); % total distance
beta=0.389;
Lcold=L_all.*beta; % total distance of cold-start emissions

%% Determine the number of flameouts, pattern==3
N_p3=0;
for i=1:m-1
    cur_p=Pattern(i);
    next_p=Pattern(i+1);
    if cur_p==3 && next_p~=3
        N_p3=N_p3+1;
    end
end

Lcold_Stp=Lcold/N_p3;

%% Determine emission type
i=2;
while i<=m
    Lstps=0;
    [prev_p, cur_p] = deal(Pattern(i-1), Pattern(i));
    if prev_p==3 && cur_p~=3  
        while cur_p~=3 && i<=m
            Lstps=Lstps+L(i);
            if Lstps<=Lcold_Stp
                EM_HorC(i)=1;
            else
                break;
            end
            
            if i<=m
                i=i+1;
            end
        end
    end

    i=i+1;
end

if size(cur_traj,2) < HorC_dim
    cur_traj = [cur_traj EM_HorC];
else
    cur_traj(:, HorC_dim) = EM_HorC;
end

end