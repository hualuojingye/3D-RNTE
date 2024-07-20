function cur_Traj=identify_MoS(cur_Traj, threshold, t_dim, l_dim, v_dim, p_dim)

format long

[m,~]=size(cur_Traj);
Cat=zeros(m,1); % Cat: 1:MA; 2:SA with engine on; 3:SA with engine-off
    
v23=cur_Traj(:,v_dim); % speed: km/h
l23=cur_Traj(:,l_dim);  % distance: m
t23=cur_Traj(:,t_dim);  % time interval: s

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
        end 
        
        SA=[];
        
        for k=i:j-1
            tsa=tsa+t23(k);
        end
        if tsa>180  
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


if size(cur_Traj, 2) < p_dim
    cur_Traj = [cur_Traj Cat];
else
    cur_Traj(:, p_dim) = Cat;
end

end
