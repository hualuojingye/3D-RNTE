tic

% Folder name for the trajectory dataset from DiDi
% [order_id, road_id, utctime, lon, lat, speed, angle, 
% x, y, elev, t, ed, ev, d_3d, v_3d]
base_path = "Folder names for trajectory files";
prefix_result = 'Result folder name';

[MoS_param, DrMo_param] = deal(15, 0.1);

idle_rate = 0.725;
M = 1456; % kg, vehicle mass
tem = 27.5; % temperature

[elev_dim, p_dim, HorC_dim, drmo_dim] = deal(10, 16, 17, 18);
[t_dim, l_dim, v_dim] = deal(11, 14, 15); 

filename_list = dir(base_path+"*.txt");
for filename_idx = 1:length(filename_list)
    cur_fn = filename_list(filename_idx).name;
    disp(cur_fn);

    Traj = load(base_path+cur_fn);

    % Data grouping
    vid = unique(Traj(:,1));
    [r,~] = size(vid);
    
    [all_FC_mat, all_CO_mat, all_NOx_mat]  = deal([], [], []);
    for k=1:r
        disp(['k= ', num2str(k)]);
        L = (Traj(:,1)==vid(k));
        cur_traj = Traj(L,:);
        
        % Recognition pattern
        cur_traj = identify_MoS(cur_traj, MoS_param, t_dim, l_dim, v_dim, p_dim);
        cur_traj = identify_EM_HorC(cur_traj, l_dim, p_dim, HorC_dim);
        cur_traj = identify_DrMo(cur_traj, DrMo_param, t_dim, v_dim, p_dim);

        % 3D-RNTE method
        [T, L, V, Pattern, HorC, DBE] = deal(cur_traj(:,t_dim), cur_traj(:,l_dim), ...
            cur_traj(:, v_dim), cur_traj(:,p_dim), cur_traj(:,HorC_dim), cur_traj(:,drmo_dim));
        
        [m, ~] = size(cur_traj);
        [FC_mat, CO_mat, NOx_mat] = deal(zeros(m,1), zeros(m,1), zeros(m,1));

        for i = 2:m
            [p, horc, dbe] = deal(Pattern(i), HorC(i), DBE(i));
            [pV, cV, cL, cT] = deal(V(i-1), V(i), L(i)/1000, T(i)); % km/h, km, s

            FC_dH = 0.0;
            % mobile activities
            if p == 1 
                dH = cur_traj(i,elev_dim) - cur_traj(i-1,elev_dim);
                FC_dH = (M*9.8*dH)/43070; % ((kg * N/kg * m) ==> J) / (J/g) 
			    dCO_H = 0.132*FC_dH;
			    dNOx_H= 0.0145*FC_dH;
                
                % hot phase
                [CO_hot, NOx_hot] = deal(0.0, 0.0);
                if dbe == 0 % cruise
                    FC_mat(i,1) = (217+0.253*cV+0.00965*cV^2)/(1+0.096*cV-0.000421*cV^2)*cL;
                    FC_mat(i,1) = FC_mat(i,1)+FC_dH;
                    CO_hot = (71.7+11.4*cV)/(1+35.4*cV-0.248*cV^2)*cL+dCO_H;
                    NOx_hot = (0.0929-0.00149*cV+0.00000653*cV^2)/(1-0.0122*cV+0.0000397*cV^2)*cL+dNOx_H;
                elseif dbe == 1 % deceleration
                    l = cL;
                    Ek = 0.3858/10000*(cV^2-pV^2)/l;
                    k1 = 0.621+0.000777*pV-0.0189*sqrt(cV);
                    [A, B, B1] = deal(30, 0.0075, 0.09);
    
                    kx = 0.046 + 100/M + 0.00421*pV+0.0026*cV;
                    ky = kx^(0.75);
                    ka = kx^(3.81)*(2-kx^(3.81));
    
                    Idle_FC = 0.361*idle_rate*cT;
                    Dec_FC = Idle_FC + (kx*A+ky*k1*B*(cV^2+pV^2)+ka*B1*M*Ek)*l*idle_rate;
                    FC_mat(i,1) = max(Idle_FC, Dec_FC);
                    FC_mat(i,1) = FC_mat(i,1)+FC_dH;
                    CO_hot = 0.132*FC_mat(i,1);  
			        NOx_hot= 0.0145*FC_mat(i,1); 
                elseif dbe == 3 % acceleration
                    l = cL;
                    Ek = 0.3858/10000*(cV^2-pV^2)/l;
                    k1 = 0.616+0.000544*cV-0.0171*sqrt(pV);
                    k2 = 1.376+0.00205*cV-0.00538*pV;
    
                    [A, B, B1, B2] = deal(30, 0.0075, 0.09, 0.045);
                    Idle_FC = 0.361*idle_rate*cT;
                    Acc_FC = Idle_FC + (A+k1*B*(cV^2+pV^2)+B1*M*Ek+k2*B2*M*Ek^2)*l*idle_rate;
                    FC_mat(i,1) = max(Idle_FC, Acc_FC);
                    FC_mat(i,1) = FC_mat(i,1)+FC_dH;
                    CO_hot = 0.132*FC_mat(i,1);  
			        NOx_hot= 0.0145*FC_mat(i,1); 
                end
                
                % cold phase
                [CO_cold, NOx_cold] = deal(0.0, 0.0);
                if horc == 1 
                    if cV>=5 && cV<25
                       CO_ratio = 5.03e-2*cV-0.363*tem+8.604;
                       NOx_ratio=4.58e-2*cV+7.47e-3*tem+0.764;
                    elseif cV>=25 && cV<=45
                       CO_ratio = 5.03e-2*cV-0.363*tem+8.604;
                       NOx_ratio=4.84e-2*cV+2.28e-3*tem+0.685;
                    else
                       [CO_ratio, NOx_ratio] = deal(2, 2);
                    end
                    CO_cold = (CO_ratio-1)*CO_hot;
                    NOx_cold = (NOx_ratio-1)*NOx_hot;
                end
                CO_mat(i,1) = CO_mat(i,1)+CO_hot+CO_cold;
                NOx_mat(i,1) = NOx_mat(i,1)+NOx_hot+NOx_cold;

            %SA with engine-on
            elseif p == 2 
                dH = cur_traj(i,elev_dim) - cur_traj(i-1,elev_dim);
                FC_dH = (M*9.8*dH)/43070; % ((kg * N/kg * m) ==> J) / (J/g)
			    dCO_H = 0.132*FC_dH;  
			    dNOx_H= 0.0145*FC_dH;  
			    
                FC_mat(i,1) = 0.361*idle_rate*cT; % ml/s * g/ml * s
			    FC_mat(i,1) = FC_mat(i,1) + FC_dH;
			    CO_hot = 13.889*idle_rate*cT;
                NOx_hot = 0.566*idle_rate*cT;

                % cold phase
                [CO_cold, NOx_cold] = deal(0.0, 0.0);
                if horc == 1 
                    if cV>=5 && cV<25
                       CO_ratio = 5.03e-2*cV-0.363*tem+8.604;
                       NOx_ratio=4.58e-2*cV+7.47e-3*tem+0.764;
                    elseif cV>=25 && cV<=45
                       CO_ratio = 5.03e-2*cV-0.363*tem+8.604;
                       NOx_ratio=4.84e-2*cV+2.28e-3*tem+0.685;
                    else
                       [CO_ratio, NOx_ratio] = deal(2, 2);
                    end
                    CO_cold = (CO_ratio-1)*CO_hot;
                    NOx_cold = (NOx_ratio-1)*NOx_hot;
                end
			    CO_mat(i,1) = CO_mat(i,1)+CO_hot+CO_cold;
                NOx_mat(i,1) = NOx_mat(i,1)+NOx_hot+NOx_cold;
            end
        end
        all_FC_mat = cat(1, all_FC_mat, FC_mat); 
		all_CO_mat = cat(1, all_CO_mat, CO_mat); 
		all_NOx_mat = cat(1, all_NOx_mat, NOx_mat); 
    end

    if size(Traj,2)<=(p_dim-1)
        Traj=cat(2, Traj, all_FC_mat);
        Traj=cat(2, Traj, all_CO_mat);
		Traj=cat(2, Traj, all_NOx_mat); 
    else
        Traj(:,drmo_dim+1) = all_FC_mat;
		Traj(:,drmo_dim+2) = all_CO_mat;
		Traj(:,drmo_dim+3) = all_NOx_mat;
    end

    out_fp = base_path + prefix_result + cur_fn;
    writematrix(Traj, out_fp);
end

toc