% =========================================================================
% calculate_receiver_position.m (全系统版: G/C/E/J + R)
% 功能: 
%   - 使用 G(GPS), C(BDS), E(GAL), J(QZSS) 和 R(GLONASS) 进行混合定位。
%   - 自动适配不同系统的伪距信号 (C1C/C2I/C1P等)。
%   - 确保不同轨道模型 (开普勒 vs 状态矢量) 的卫星位置均被正确引入。
% =========================================================================

function [receiver_pos, receiver_clk_err_sec, available_sat_states, lat_deg, lon_deg, alt_m] = calculate_receiver_position(obs_data, nav_data, epoch_idx)
    
    % 物理常数
    c = 299792458.0; 
    omega_e = 7.2921151467e-5; % 地球自转角速度 (WGS-84)
    
    % 初始化
    max_iterations = 10; 
    tolerance = 1e-4;
    receiver_pos_approx = [0; 0; 0]; 
    receiver_clk_m_approx = 0;
    
    available_sat_states = struct(); 
    current_epoch = obs_data(epoch_idx);
    t_obs = current_epoch.time;
    observed_sats_all = fieldnames(current_epoch.data);
    
    % --- 第一次遍历：计算所有可用卫星的状态 ---
    sats_for_positioning = {};
    
    for i = 1:length(observed_sats_all)
        sat_id = observed_sats_all{i};
        sys = upper(sat_id(1));

        % [关键修改] 增加对 'R' (GLONASS) 的支持
        if ~(sys == 'G' || sys == 'C' || sys == 'E' || sys == 'J' || sys == 'R')
            continue;
        end
        
        sat_obs = current_epoch.data.(sat_id);
        
        % --- 伪距信号选择逻辑 (Priority Selector) ---
        pr_code = '';
        
        if (sys == 'G' || sys == 'E' || sys == 'J')
            % GPS/Galileo/QZSS: 优先 L1 C/A (C1C)
            if isfield(sat_obs.pseudorange, 'C1C')
                pr_code = 'C1C';
            elseif isfield(sat_obs.pseudorange, 'C1X'), pr_code = 'C1X'; 
            end
            
        elseif sys == 'C'
            % BDS: 优先 B1I (C2I)
            if isfield(sat_obs.pseudorange, 'C2I')
                pr_code = 'C2I';
            elseif isfield(sat_obs.pseudorange, 'C1I'), pr_code = 'C1I'; % 有些旧文件可能是 C1I
            end
            
        elseif sys == 'R'
            % [新增] GLONASS: 优先 C1C (L1 C/A) 或 C1P (L1 P)
            if isfield(sat_obs.pseudorange, 'C1C')
                pr_code = 'C1C';
            elseif isfield(sat_obs.pseudorange, 'C1P')
                pr_code = 'C1P';
            end
        end
        
        % 兜底: 如果首选都没找到，取第一个非空的伪距
        if isempty(pr_code)
            avail_codes = fieldnames(sat_obs.pseudorange);
            for k=1:length(avail_codes)
                val = sat_obs.pseudorange.(avail_codes{k});
                if ~isnan(val) && val > 0
                    pr_code = avail_codes{k};
                    break;
                end
            end
        end
        
        if isempty(pr_code), continue; end
        pseudorange = sat_obs.pseudorange.(pr_code);
        if isnan(pseudorange) || pseudorange <= 0, continue; end
        
        try
            % 调用 calculate_satellite_state
            % 注意: 该函数内部已经封装了 GLONASS 的特殊轨道算法 (Runge-Kutta/外推)
            [sat_pos, sat_vel, sat_clk_err] = ...
                calculate_satellite_state(t_obs, pseudorange, sat_id, nav_data);
            
            available_sat_states.(sat_id).position = sat_pos;
            available_sat_states.(sat_id).velocity = sat_vel;
            available_sat_states.(sat_id).clock_error = sat_clk_err;
            available_sat_states.(sat_id).pseudorange_raw = pseudorange;
            sats_for_positioning{end+1} = sat_id;
        catch ME
            % 静默跳过无法计算轨道的卫星
            continue;
        end
    end
    
    % --- 第二次遍历：最小二乘解算 (Least Squares) ---
    num_sats = length(sats_for_positioning);
    if num_sats < 4
        receiver_pos = NaN(3,1); receiver_clk_err_sec = NaN;
        lat_deg = NaN; lon_deg = NaN; alt_m = NaN;
        return;
    end
    
    for iter = 1:max_iterations
        H = zeros(num_sats, 4); 
        L = zeros(num_sats, 1);
        
        for i = 1:num_sats
            sat_id = sats_for_positioning{i};
            sat_state = available_sat_states.(sat_id);
            
            sat_pos = sat_state.position;       % ECEF 坐标 (已包含轨道模型差异的处理)
            sat_clk_err = sat_state.clock_error;
            pr_raw = sat_state.pseudorange_raw;
            
            % 1. 几何距离计算
            geom_range_raw = norm(sat_pos - receiver_pos_approx);
            
            % 2. 地球自转校正 (Sagnac Effect Correction)
            % 这一步对于所有系统的卫星（包括 GLONASS）都是通用的，因为 sat_pos 都是 ECEF
            travel_time = geom_range_raw / c; 
            theta = omega_e * travel_time;
            
            R_z = [cos(theta), sin(theta), 0;
                  -sin(theta), cos(theta), 0;
                   0,          0,          1];
            
            sat_pos_rot = R_z * sat_pos;
            
            % 3. 更新后的几何距离与方向矢量
            geom_range = norm(sat_pos_rot - receiver_pos_approx);
            los_vec = (sat_pos_rot - receiver_pos_approx) / geom_range;
            
            % 4. 组建矩阵 (观测方程线性化)
            % P = rho + c*(dt_r - dt_s) + I + T + epsilon
            % l_i = P_observed - P_computed
            % P_computed = geom_range + clk_bias_r - c*clk_err_s
            l_i = pr_raw - (geom_range + receiver_clk_m_approx - c * sat_clk_err);
            
            H(i, :) = [-los_vec', 1];
            L(i, 1) = l_i;
        end
        
        % 求解增量 dx
        if rank(H) < 4
            break; % 矩阵奇异
        end
        dx = H \ L;
        
        receiver_pos_approx = receiver_pos_approx + dx(1:3);
        receiver_clk_m_approx = receiver_clk_m_approx + dx(4);
        
        if norm(dx(1:3)) < tolerance
            break;
        end
    end
    
    receiver_pos = receiver_pos_approx;
    receiver_clk_err_sec = receiver_clk_m_approx / c;
    
    % --- 坐标转换 (ECEF -> Geodetic) ---
    [lat_rad, lon_rad, alt_m] = ecef2geodetic_custom(receiver_pos(1), receiver_pos(2), receiver_pos(3));
    lat_deg = rad2deg(lat_rad);
    lon_deg = rad2deg(lon_rad);
end

% 内部辅助函数: ECEF转大地坐标 (防止依赖外部文件缺失)
function [phi, lambda, h] = ecef2geodetic_custom(x, y, z)
    a = 6378137.0; 
    f = 1/298.257223563; 
    b = a * (1-f);
    e2 = f * (2-f);
    ep2 = (a^2 - b^2) / b^2;
    p = sqrt(x^2 + y^2);
    theta = atan2(z*a, p*b);
    lon = atan2(y, x);
    lat = atan2(z + ep2*b*sin(theta)^3, p - e2*a*cos(theta)^3);
    N = a / sqrt(1 - e2*sin(lat)^2);
    alt = p / cos(lat) - N;
    phi = lat; lambda = lon; h = alt;
end











% % ============== calculate_receiver_position.m (新版: 支持 G/C/E/J 定位) ==============
% function [receiver_pos, receiver_clk_err_sec, available_sat_states, lat_deg, lon_deg, alt_m] = calculate_receiver_position(obs_data, nav_data, epoch_idx)
%     % 功能:
%     %   - 使用 G, C, E, J 四系统的观测数据进行最小二乘定位。
%     %   - 依赖: calculate_satellite_state.m (必须是支持 G/C/E/J 的新版本)
%     
%     c = 299792458.0; max_iterations = 10; tolerance = 1e-4;
%     receiver_pos_approx = [0; 0; 0]; receiver_clk_m_approx = 0;
%     
%     available_sat_states = struct(); 
%     current_epoch = obs_data(epoch_idx);
%     t_obs = current_epoch.time;
%     observed_sats_all = fieldnames(current_epoch.data);
%     
%     % --- 第一次遍历：计算所有可用卫星的状态 (G/C/E/J) ---
%     sats_for_positioning = {};
%     for i = 1:length(observed_sats_all)
%         sat_id = observed_sats_all{i};
%         sys = upper(sat_id(1));
% 
%         % 检查是否是G, C, E, J 系统
%         if ~(sys == 'G' || sys == 'C' || sys == 'E' || sys == 'J')
%             continue;
%         end
%         
%         sat_obs = current_epoch.data.(sat_id);
%         
%         % 根据不同系统选择合适的伪距码
%         % G (GPS L1)     -> C1C
%         % E (Galileo E1) -> C1C
%         % J (QZSS L1)    -> C1C
%         % C (BDS B1I)    -> C2I
%         pr_code = '';
%         if (sys == 'G' || sys == 'E' || sys == 'J')
%             if isfield(sat_obs.pseudorange, 'C1C')
%                 pr_code = 'C1C';
%             % (如果您的 E/J 文件使用 C1X，可以在这里添加 elseif)
%             end
%         elseif sys == 'C'
%             if isfield(sat_obs.pseudorange, 'C2I')
%                 pr_code = 'C2I';
%             % (如果您的 BDS 文件使用 C1X，可以在这里添加 elseif)
%             end
%         end
%         
%         if isempty(pr_code), continue; end
% 
%         pseudorange = sat_obs.pseudorange.(pr_code);
%         if isnan(pseudorange), continue; end
%         
%         try
%             % !! 关键: 调用我们刚刚替换的、支持多系统的 'calculate_satellite_state.m'
%             [sat_pos, sat_vel, sat_clk_err] = ...
%                 calculate_satellite_state(t_obs, pseudorange, sat_id, nav_data);
%             
%             available_sat_states.(sat_id).position = sat_pos;
%             available_sat_states.(sat_id).velocity = sat_vel;
%             available_sat_states.(sat_id).clock_error = sat_clk_err;
%             available_sat_states.(sat_id).pseudorange_raw = pseudorange;
%             sats_for_positioning{end+1} = sat_id;
%         catch ME
%             % fprintf('⚠️  在定位中跳过卫星 %s 的状态计算: %s\n', sat_id, ME.message);
%             continue;
%         end
%     end
%     
%     % --- 第二次遍历：使用状态可用的卫星进行定位 ---
%     sats_for_positioning = fieldnames(available_sat_states);
%     if length(sats_for_positioning) < 4
%         receiver_pos = NaN(3,1); receiver_clk_err_sec = NaN;
%         lat_deg = NaN; lon_deg = NaN; alt_m = NaN;
%         return;
%     end
%     
%     for iter = 1:max_iterations
%         H = []; L = []; 
%         for i = 1:length(sats_for_positioning)
%             sat_id = sats_for_positioning{i};
%             sat_state = available_sat_states.(sat_id);
%             sat_pos = sat_state.position;
%             sat_clk_err = sat_state.clock_error;
%             pseudorange = sat_state.pseudorange_raw;
%             
%             omega_e = 7.2921151467e-5;
%             travel_time = pseudorange / c;
%             R_e = [cos(omega_e * travel_time), sin(omega_e * travel_time), 0;
%                   -sin(omega_e * travel_time), cos(omega_e * travel_time), 0;
%                    0, 0, 1];
%             sat_pos_rotated = R_e * sat_pos;
% 
%             geom_range = norm(sat_pos_rotated - receiver_pos_approx);
%             l_i = pseudorange - (geom_range + receiver_clk_m_approx - c * sat_clk_err);
%             los_vec = (sat_pos_rotated - receiver_pos_approx) / geom_range;
% 
%             h_i = [-los_vec', 1];
%             H(i, :) = h_i;
%             L(i, 1) = l_i;
%         end
%         
%         try
%             dx = H \ L;
%         catch
%             % 矩阵奇异，解算失败
%             receiver_pos = NaN(3,1); receiver_clk_err_sec = NaN;
%             lat_deg = NaN; lon_deg = NaN; alt_m = NaN;
%             return;
%         end
%         
%         receiver_pos_approx = receiver_pos_approx + dx(1:3);
%         receiver_clk_m_approx = receiver_clk_m_approx + dx(4);
%         
%         if norm(dx(1:3)) < tolerance
%             break;
%         end
%     end
%     receiver_pos = receiver_pos_approx;
%     receiver_clk_err_sec = receiver_clk_m_approx / c;
%     
%     % --- 坐标转换 (ECEF -> Geodetic) ---
%     a = 6378137.0; f = 1 / 298.257223563; e_sq = f * (2 - f);
%     X = receiver_pos(1); Y = receiver_pos(2); Z = receiver_pos(3);
%     lon_rad = atan2(Y, X);
%     p = sqrt(X^2 + Y^2);
%     lat_rad_old = atan2(Z, p * (1 - e_sq));
%     tolerance_lat = 1e-12; delta_lat = tolerance_lat + 1;
%     while delta_lat > tolerance_lat
%         N = a / sqrt(1 - e_sq * sin(lat_rad_old)^2);
%         lat_rad_new = atan2(Z + N * e_sq * sin(lat_rad_old), p);
%         delta_lat = abs(lat_rad_new - lat_rad_old);
%         lat_rad_old = lat_rad_new;
%     end
%     lat_rad = lat_rad_new;
%     N = a / sqrt(1 - e_sq * sin(lat_rad)^2);
%     alt_m = (p / cos(lat_rad)) - N;
%     lat_deg = rad2deg(lat_rad);
%     lon_deg = rad2deg(lon_rad);
% end