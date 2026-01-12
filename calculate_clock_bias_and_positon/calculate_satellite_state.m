% ============== calculate_satellite_state.m (v2.0: 增加 GLONASS 支持) ==============
function [sat_pos, sat_vel, sat_clk_err] = calculate_satellite_state(t_obs, pseudorange, satellite_id, nav_data)
% CALCULATE_SATELLITE_STATE - 计算卫星位置、速度和钟差。
%
% 更新: 增加了对 GLONASS (R) 的运动学外推支持，用于天空图绘制。
% 适配: parse_rinex_nav_multi_gnss (G=1, C=2, R=3, E=4, J=5)

% --- 0. 根据卫星ID判断系统并获取星历 ---
sys_char = upper(satellite_id(1));
prn = str2double(satellite_id(2:end));

% 系统索引映射
switch sys_char
    case 'G', sys_idx = 1;
    case 'C', sys_idx = 2;
    case 'R', sys_idx = 3;
    case 'E', sys_idx = 4;
    case 'J', sys_idx = 5;
    otherwise
        error('不支持的卫星系统: %s', satellite_id);
end

if prn < 1 || prn > size(nav_data, 1) || sys_idx > size(nav_data, 2) || isempty(nav_data{prn, sys_idx})
    % 静默失败，由上层处理
    error('未找到卫星 %s 的星历数据。', satellite_id);
end

eph_sets = nav_data{prn, sys_idx};

% --- 物理常数 ---
c = 299792458.0; 

% 计算观测时刻 (Second of Week / Second of Day)
[~,~,~,H,M,S] = datevec(t_obs);
% 注意：GLONASS 的参考时间通常是当天的 UTC 时间，而 GPS/BDS 是周内秒。
% 为了统一，我们将 t_obs 转换为该星历参考系下的相对时间 dt。

% =========================================================================
%  分支 A: GLONASS (R) 处理逻辑 (运动学外推)
% =========================================================================
if sys_char == 'R'
    % 1. 寻找最近的星历 (基于 Toc)
    % GLONASS 的 Toc 在 parse_rinex_nav_multi_gnss 中已解析为 struct
    min_diff = inf; best_eph_index = -1;
    
    % 将 t_obs 转为 serial date number 以便比较
    t_obs_num = datenum(t_obs);
    
    for i = 1:length(eph_sets)
        % 构建星历参考时间 (Toc)
        toc_vec = [eph_sets(i).Toc.Year, eph_sets(i).Toc.Month, eph_sets(i).Toc.Day, ...
                   eph_sets(i).Toc.Hour, eph_sets(i).Toc.Minute, eph_sets(i).Toc.Second];
        t_toc_num = datenum(toc_vec);
        
        diff_sec = abs((t_obs_num - t_toc_num) * 86400);
        if diff_sec < min_diff
            min_diff = diff_sec;
            best_eph_index = i;
        end
    end
    
    if best_eph_index == -1
         error('未找到 GLONASS 卫星 %s 的合适星历。', satellite_id);
    end
    eph = eph_sets(best_eph_index);
    
    % 2. 计算时间差 dt (考虑信号传播时间)
    transit_time = pseudorange / c;
    t_signal_arrival = datenum(t_obs);
    t_signal_transmit = t_signal_arrival - (transit_time / 86400);
    
    toc_vec = [eph.Toc.Year, eph.Toc.Month, eph.Toc.Day, ...
               eph.Toc.Hour, eph.Toc.Minute, eph.Toc.Second];
    t_ref = datenum(toc_vec);
    
    % dt = 发射时刻 - 参考时刻 (秒)
    dt = (t_signal_transmit - t_ref) * 86400;
    
    % 3. 检查 dt 有效性 (GLONASS 星历有效期较短，通常 +-15分钟)
    % 但为了画图，我们放宽到 +- 1小时
    if abs(dt) > 3600
        % warning('GLONASS %s 时间差过大 (%.1f 秒)，计算可能不准。', satellite_id, dt);
    end
    
    % 4. 运动学外推位置 (P = P0 + V*t + 0.5*A*t^2)
    % 注意：RINEX 导航文件中的坐标是 PZ-90，与 WGS-84 差异很小(米级)，画图可忽略转换。
    pos_0 = [eph.pos_x; eph.pos_y; eph.pos_z] * 1000; % km -> m
    vel_0 = [eph.vel_x; eph.vel_y; eph.vel_z] * 1000; % km/s -> m/s
    acc_0 = [eph.acc_x; eph.acc_y; eph.acc_z] * 1000; % km/s^2 -> m/s^2
    
    sat_pos = pos_0 + vel_0 * dt + 0.5 * acc_0 * (dt^2);
    sat_vel = vel_0 + acc_0 * dt; % 简单速度外推
    
    % 5. 钟差计算 (Clk = -tau_n + gamma_n * dt)
    % 注意符号：RINEX 中 tau_n 是卫星钟相对于系统时的偏差
    sat_clk_err = -eph.tau_n + eph.gamma_n * dt;
    
    return; % GLONASS 计算结束
end

% =========================================================================
%  分支 B: GPS / BDS / Galileo / QZSS (开普勒轨道)
% =========================================================================

% --- 物理常数 ---
omega_e_dot = 7.2921151467e-5;
F_rel = -4.442807633e-10;
if sys_char == 'C'
    mu_E = 3.986004418e14; % CGCS2000
else
    mu_E = 3.986005e14;   % WGS-84
end

% --- 寻找最佳星历 (基于 Toe / SOW) ---
t_obs_sow = (weekday(t_obs)-1)*86400 + H*3600 + M*60 + S;

min_diff = inf; best_eph_index = -1;
for i = 1:length(eph_sets)
    time_diff = abs(t_obs_sow - eph_sets(i).Toe);
    if time_diff > 302400, time_diff = 604800 - time_diff; end
    if time_diff < min_diff
        min_diff = time_diff;
        best_eph_index = i;
    end
end

if best_eph_index == -1
    error('未能为卫星 %s 找到合适的星历。', satellite_id);
end
eph = eph_sets(best_eph_index);

% --- 计算信号发射时刻 ---
transit_time = pseudorange / c;
t_transmit_sow = t_obs_sow - transit_time;

% --- 轨道参数计算 ---
A = eph.sqrtA^2; e = eph.e;
tk = t_transmit_sow - eph.Toe;
if tk > 302400, tk = tk - 604800; end
if tk < -302400, tk = tk + 604800; end

n0 = sqrt(mu_E / A^3); n = n0 + eph.Delta_n; Mk = eph.M0 + n * tk;
Ek = Mk; for i=1:10, Ek = Mk + e * sin(Ek); end
vk = atan2(sqrt(1-e^2)*sin(Ek), cos(Ek)-e); phik = vk + eph.omega;

delta_uk = eph.Cus*sin(2*phik) + eph.Cuc*cos(2*phik); 
delta_rk = eph.Crs*sin(2*phik) + eph.Crc*cos(2*phik); 
delta_ik = eph.Cis*sin(2*phik) + eph.Cic*cos(2*phik);

uk = phik + delta_uk; 
rk = A*(1-e*cos(Ek)) + delta_rk; 
ik = eph.i0 + delta_ik + eph.IDOT * tk;

xk_prime = rk*cos(uk); yk_prime = rk*sin(uk); 
Omegak = eph.OMEGA0 + (eph.OMEGA_DOT - omega_e_dot)*tk - omega_e_dot*eph.Toe;

Xk = xk_prime*cos(Omegak) - yk_prime*cos(ik)*sin(Omegak); 
Yk = xk_prime*sin(Omegak) + yk_prime*cos(ik)*cos(Omegak); 
Zk = yk_prime*sin(ik);

sat_pos = [Xk; Yk; Zk];

% --- 速度计算 (近似) ---
Ek_dot = n/(1-e*cos(Ek)); 
vk_dot = Ek_dot*sqrt(1-e^2)/(1-e*cos(Ek));
uk_dot = vk_dot + 2*(eph.Cus*cos(2*phik)-eph.Cuc*sin(2*phik))*vk_dot; 
rk_dot = A*e*sin(Ek)*Ek_dot + 2*(eph.Crs*cos(2*phik)-eph.Crc*sin(2*phik))*vk_dot; 
ik_dot = eph.IDOT + 2*(eph.Cis*cos(2*phik)-eph.Cic*sin(2*phik))*vk_dot; 
Omegak_dot = eph.OMEGA_DOT - omega_e_dot;
xk_prime_dot = rk_dot*cos(uk)-yk_prime*uk_dot; 
yk_prime_dot = rk_dot*sin(uk)+xk_prime*uk_dot;
Vx=xk_prime_dot*cos(Omegak)-yk_prime_dot*cos(ik)*sin(Omegak)+yk_prime*sin(ik)*sin(Omegak)*ik_dot-Yk*Omegak_dot; 
Vy=xk_prime_dot*sin(Omegak)+yk_prime_dot*cos(ik)*cos(Omegak)-yk_prime*sin(ik)*cos(Omegak)*ik_dot+Xk*Omegak_dot; 
Vz=yk_prime_dot*sin(ik)+yk_prime*cos(ik)*ik_dot;
sat_vel = [Vx; Vy; Vz];

% --- 钟差计算 ---
dtr = F_rel*e*eph.sqrtA*sin(Ek); 
dt_clock = t_transmit_sow - eph.Toe;
if dt_clock > 302400, dt_clock = dt_clock - 604800; end
if dt_clock < -302400, dt_clock = dt_clock + 604800; end
        
switch sys_char
    case {'G', 'J'} 
        dtsv_poly = eph.af0 + eph.af1*dt_clock + eph.af2*(dt_clock^2);
        sat_clk_err = dtsv_poly + dtr - eph.TGD;
    case 'E' 
        dtsv_poly = eph.af0 + eph.af1*dt_clock + eph.af2*(dt_clock^2);
        sat_clk_err = dtsv_poly + dtr - eph.TGD_E1E5a; 
    case 'C' 
        toc_datetime = datetime(eph.Toc.Year,eph.Toc.Month,eph.Toc.Day,eph.Toc.Hour,eph.Toc.Minute,eph.Toc.Second);
        [~,~,~,H_toc,M_toc,S_toc] = datevec(toc_datetime);
        t_toc_sow = (weekday(toc_datetime)-1)*86400 + H_toc*3600 + M_toc*60 + S_toc;
        dt_clock_bds = t_transmit_sow - t_toc_sow;
        if dt_clock_bds > 302400, dt_clock_bds = dt_clock_bds - 604800; end
        if dt_clock_bds < -302400, dt_clock_bds = dt_clock_bds + 604800; end
        dtsv_poly = eph.A0 + eph.A1*dt_clock_bds + eph.A2*(dt_clock_bds^2);
        sat_clk_err = dtsv_poly + dtr - eph.TGD1;
end
end










% % ============== calculate_satellite_state.m (新版: 适配 G/C/R/E/J) ==============
% function [sat_pos, sat_vel, sat_clk_err] = calculate_satellite_state(t_obs, pseudorange, satellite_id, nav_data)
% % CALCULATE_SATELLITE_STATE - 计算卫星位置、速度和钟差。
% %
% % 适配 `parse_rinex_nav_gps_bds` (G=1, C=2, R=3, E=4, J=5) 的数据结构
% 
% % --- 0. 根据卫星ID判断系统并获取星历 ---
% sys_char = upper(satellite_id(1));
% prn = str2double(satellite_id(2:end));
% 
% % 系统索引映射 (必须与解析器 parse_rinex_nav_gps_bds 严格一致!)
% switch sys_char
%     case 'G', sys_idx = 1;
%     case 'C', sys_idx = 2;
%     case 'R', sys_idx = 3;
%     case 'E', sys_idx = 4;
%     case 'J', sys_idx = 5;
%     otherwise
%         error('不支持的卫星系统: %s', satellite_id);
% end
% 
% if prn < 1 || prn > size(nav_data, 1) || sys_idx > size(nav_data, 2) || isempty(nav_data{prn, sys_idx})
%     error('在导航数据中未找到卫星 %s (PRN %d, SysIdx %d) 的星历。', satellite_id, prn, sys_idx);
% end
% 
% eph_sets = nav_data{prn, sys_idx};
% 
% % --- GLONASS (R) 的特殊处理 ---
% if sys_char == 'R'
%     error('GLONASS (R) 系统的轨道计算逻辑 (Runge-Kutta 积分) 未实现。');
% end
% 
% % --- G/E/C/J 的 "最佳星历" 选择 (基于 Toe) ---
% [~,~,~,H,M,S] = datevec(t_obs);
% t_obs_sow = (weekday(t_obs)-1)*86400 + H*3600 + M*60 + S;
% 
% min_diff = inf; best_eph_index = -1;
% for i = 1:length(eph_sets)
%     time_diff = abs(t_obs_sow - eph_sets(i).Toe);
%     if time_diff > 302400, time_diff = 604800 - time_diff; end
%     if time_diff < min_diff
%         min_diff = time_diff;
%         best_eph_index = i;
%     end
% end
% 
% if best_eph_index == -1
%     error('未能为卫星 %s 在时刻 %s 找到合适的星历。', satellite_id, datestr(t_obs));
% end
% eph = eph_sets(best_eph_index);
% 
% % --- 物理常数 ---
% c = 299792458.0; 
% omega_e_dot = 7.2921151467e-5;
% F_rel = -4.442807633e-10;
% 
% if sys_char == 'C'
%     mu_E = 3.986004418e14; % CGCS2000 for BeiDou
% else
%     mu_E = 3.986005e14;   % WGS-84 for GPS/Galileo/QZSS
% end
% 
% % --- 计算信号发射时刻 ---
% transit_time = pseudorange / c;
% t_transmit_sow = t_obs_sow - transit_time;
% 
% % --- 轨道参数 (G/E/C/J 通用) ---
% A = eph.sqrtA^2; e = eph.e;
% tk = t_transmit_sow - eph.Toe;
% if tk > 302400, tk = tk - 604800; end
% if tk < -302400, tk = tk + 604800; end
% 
% % --- 1. 计算卫星位置 (G/E/C/J 通用) ---
% n0 = sqrt(mu_E / A^3); n = n0 + eph.Delta_n; Mk = eph.M0 + n * tk;
% Ek = Mk; for i=1:10, Ek = Mk + e * sin(Ek); end
% vk = atan2(sqrt(1-e^2)*sin(Ek), cos(Ek)-e); phik = vk + eph.omega;
% delta_uk = eph.Cus*sin(2*phik) + eph.Cuc*cos(2*phik); 
% delta_rk = eph.Crs*sin(2*phik) + eph.Crc*cos(2*phik); 
% delta_ik = eph.Cis*sin(2*phik) + eph.Cic*cos(2*phik);
% uk = phik + delta_uk; 
% rk = A*(1-e*cos(Ek)) + delta_rk; 
% ik = eph.i0 + delta_ik + eph.IDOT * tk;
% xk_prime = rk*cos(uk); yk_prime = rk*sin(uk); 
% Omegak = eph.OMEGA0 + (eph.OMEGA_DOT - omega_e_dot)*tk - omega_e_dot*eph.Toe;
% Xk = xk_prime*cos(Omegak) - yk_prime*cos(ik)*sin(Omegak); 
% Yk = xk_prime*sin(Omegak) + yk_prime*cos(ik)*cos(Omegak); 
% Zk = yk_prime*sin(ik);
% sat_pos = [Xk; Yk; Zk];
% 
% % --- 2. 计算卫星速度 (G/E/C/J 通用) ---
% Ek_dot = n/(1-e*cos(Ek)); 
% vk_dot = Ek_dot*sqrt(1-e^2)/(1-e*cos(Ek));
% uk_dot = vk_dot + 2*(eph.Cus*cos(2*phik)-eph.Cuc*sin(2*phik))*vk_dot; 
% rk_dot = A*e*sin(Ek)*Ek_dot + 2*(eph.Crs*cos(2*phik)-eph.Crc*sin(2*phik))*vk_dot; 
% ik_dot = eph.IDOT + 2*(eph.Cis*cos(2*phik)-eph.Cic*sin(2*phik))*vk_dot; 
% Omegak_dot = eph.OMEGA_DOT - omega_e_dot;
% xk_prime_dot = rk_dot*cos(uk)-yk_prime*uk_dot; 
% yk_prime_dot = rk_dot*sin(uk)+xk_prime*uk_dot;
% Vx=xk_prime_dot*cos(Omegak)-yk_prime_dot*cos(ik)*sin(Omegak)+yk_prime*sin(ik)*sin(Omegak)*ik_dot-Yk*Omegak_dot; 
% Vy=xk_prime_dot*sin(Omegak)+yk_prime_dot*cos(ik)*cos(Omegak)-yk_prime*sin(ik)*cos(Omegak)*ik_dot+Xk*Omegak_dot; 
% Vz=yk_prime_dot*sin(ik)+yk_prime*cos(ik)*ik_dot;
% sat_vel = [Vx; Vy; Vz];
% 
% % --- 3. 计算卫星钟差 (G/E/C/J) ---
% dtr = F_rel*e*eph.sqrtA*sin(Ek); % 相对论效应 (通用)
% 
% dt_clock = t_transmit_sow - eph.Toe; % G, E, J 使用 Toe
% if dt_clock > 302400, dt_clock = dt_clock - 604800; end
% if dt_clock < -302400, dt_clock = dt_clock + 604800; end
%         
% switch sys_char
%     case {'G', 'J'} % GPS 和 QZSS (使用 af0, af1, af2 和 TGD)
%         dtsv_poly = eph.af0 + eph.af1*dt_clock + eph.af2*(dt_clock^2);
%         sat_clk_err = dtsv_poly + dtr - eph.TGD;
%         
%     case 'E' % Galileo (使用 af0, af1, af2 和 TGD_E1E5a/b)
%         dtsv_poly = eph.af0 + eph.af1*dt_clock + eph.af2*(dt_clock^2);
%         sat_clk_err = dtsv_poly + dtr - eph.TGD_E1E5a; % 假设使用 E1/E5a (C1C)
%         
%     case 'C' % BeiDou (使用 A0, A1, A2 和 TGD1/2)
%         toc_datetime = datetime(eph.Toc.Year,eph.Toc.Month,eph.Toc.Day,eph.Toc.Hour,eph.Toc.Minute,eph.Toc.Second);
%         [~,~,~,H_toc,M_toc,S_toc] = datevec(toc_datetime);
%         t_toc_sow = (weekday(toc_datetime)-1)*86400 + H_toc*3600 + M_toc*60 + S_toc;
%         dt_clock_bds = t_transmit_sow - t_toc_sow;
%         if dt_clock_bds > 302400, dt_clock_bds = dt_clock_bds - 604800; end
%         if dt_clock_bds < -302400, dt_clock_bds = dt_clock_bds + 604800; end
%         dtsv_poly = eph.A0 + eph.A1*dt_clock_bds + eph.A2*(dt_clock_bds^2);
%         sat_clk_err = dtsv_poly + dtr - eph.TGD1; % B1I频点(C2I)使用TGD1
% end
% end