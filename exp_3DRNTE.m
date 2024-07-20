% clc
% tic

% Folder names for all trajectory files from January 9-13, 2020
% [VID, utctime, mmX, mmY, Roadid, Trajid, X_proj, Y_proj, T, L, V, elev, Pattern, HorC, DrMo]
mm_3d_bfp = "Folder names for trajectory files";

% trajectory files
[me_idx, mt_idx, ml_idx, mv_idx] = deal(12, 9, 10, 11);
[mp_idx, mhc_idx, md_idx] = deal(13, 14, 15);

% parameters of activities and driving modes
[MoS_param, DrMo_param] = deal(15, 0.1);

% real fuel consumption of trajectory, kg
real_FC = 60.46;

idle_rate = 0.725; % g/ml
M = 1456; % vehicle mass: kg

%% baseline + road shape + slope + 3D speed + driving modes, correspongding to mm_3d_bfp

format long

FC = 0.0;
FC_dH = 0.0; % fuel consumption of road slope
filename_list = dir(mm_3d_bfp+"*.txt");
for filename_idx = 1:length(filename_list)
    Traj = load(mm_3d_bfp+filename_list(filename_idx).name);
    
    % mm_3d_bfp
    Traj = identify_MoS(Traj, MoS_param, mt_idx, ml_idx, mv_idx, mp_idx);
    Traj = identify_EM_HorC(Traj, ml_idx, mp_idx, mhc_idx);
    Traj = identify_DrMo(Traj, DrMo_param, mt_idx, mv_idx, mp_idx);
    
    [T, L, V, Pattern, DBE] = deal(Traj(:,mt_idx), Traj(:,ml_idx), Traj(:,mv_idx), Traj(:,mp_idx),Traj(:,md_idx));
    [m, ~] = size(Traj);
    for i = 1:m
        [p, dbe] = deal(Pattern(i), DBE(i));
        if p == 1
            % road slope
            dH = Traj(i,me_idx) - Traj(i-1,me_idx);
            FC_dH = FC_dH + (M*9.8*dH)/43070; % ((kg * N/kg * m) ==> J) / (J/g)

            if dbe == 0 % cruise
                Cruise_FC = (217+0.253*V(i)+0.00965*V(i)^2)/(1+0.096*V(i)-0.000421*V(i)^2)*L(i)/1000;
                FC = FC + Cruise_FC;  
            elseif dbe == 1 % decceleration
                l = L(i)/1000;
                Ek = 0.3858/10000*(V(i)^2-V(i-1)^2)/l;
                k1 = 0.621+0.000777*V(i-1)-0.0189*sqrt(V(i));
                [A, B, B1] = deal(30, 0.0075, 0.09);

                kx = 0.046 + 100/M + 0.00421*V(i-1)+0.0026*V(i);
                ky = kx^(0.75);
                ka = kx^(3.81)*(2-kx^(3.81));

                Idle_FC = 0.361*idle_rate*T(i);
                Dec_FC = Idle_FC + (kx*A+ky*k1*B*(V(i)^2+V(i-1)^2)+ka*B1*M*Ek)*l*idle_rate;
                FC = FC+max(Idle_FC, Dec_FC);
            elseif dbe == 2 % idle
                Idle_FC = 0.361*idle_rate*T(i);
                FC = FC + Idle_FC;
            elseif dbe == 3 % acceleration
                l = L(i)/1000;
                Ek = 0.3858/10000*(V(i)^2-V(i-1)^2)/l;
                k1 = 0.616+0.000544*V(i)-0.0171*sqrt(V(i-1));
                k2 = 1.376+0.00205*V(i)-0.00538*V(i-1);

                [A, B, B1, B2] = deal(30, 0.0075, 0.09, 0.045);
                Idle_FC = 0.361*idle_rate*T(i);
                Acc_FC = Idle_FC + (A+k1*B*(V(i)^2+V(i-1)^2)+B1*M*Ek+k2*B2*M*Ek^2)*l*idle_rate;
                FC = FC+max(Idle_FC, Acc_FC);
            end 
        elseif p == 2 % SA with engine-on
            Idle_FC = 0.361*idle_rate*T(i); % ml/s * g/ml * s
            FC = FC + Idle_FC; 

            % road slope
            dH = Traj(i,me_idx) - Traj(i-1,me_idx);
            FC_dH = FC_dH + (M*9.8*dH)/43070; % ((kg * N/kg * m) ==> J) / (J/g)
        end 
    end
end

FC_all = FC + FC_dH;
Acc = roundn((1-abs(FC_all/1000-real_FC)/real_FC)*100, -2);
disp(['3D-RNTE: FC_all = ', num2str(FC_all/1000), ', ', 'Acc = ', num2str(Acc), ' %, ']);

