% =========================================================================
% simulate_gnss_spoofing.m (函数版)
% 功能: 基于"物理特征缺失"的攻击仿真 (v2.0)
% 核心思想: 
%   1. SDR攻击: 模拟单源欺骗。只有在SDR发射机波束范围内的卫星保留波动(模拟遮挡)，
%      其余方向卫星因物理缺失，信号被替换为平滑直线(环境噪声)。
%   2. 重放攻击: 模拟无动作重放。所有卫星信号均为平滑直线，完全没有当前手势的物理遮挡特征。
%
% [调用格式]:
%   obs_spoofed = simulate_gnss_spoofing(obs_data, nav_data, ATTACK_TYPE);
%
% [输入参数]:
%   1. obs_data: 原始观测数据。
%   2. nav_data: 导航星历数据(用于计算卫星几何分布)。
%   3. ATTACK_TYPE: 攻击模式字符串。
%      - 'SDR'   : 单源物理遮挡攻击 (保留特定方向波动)。
%      - 'REPLAY': 全星座重放攻击 (全部抹平)。
%
% [返回值]:
%   obs_spoofed: 被注入欺骗特征后的观测数据。
% =========================================================================

function obs_spoofed = simulate_gnss_spoofing(obs_data, nav_data, ATTACK_TYPE)

addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('sky_plot'));

fprintf('--> 启动攻击仿真 (Function版): Mode = [%s]...\n', ATTACK_TYPE);

%% 1. 基础设置
% 初始化输出数据 (复制一份，避免修改原数据)
obs_spoofed = obs_data;

% --- 内部参数配置 ---
SDR_AZIMUTH    = 120;   % [SDR参数] 发射机方位角 (度)
SDR_BEAM_WIDTH = 30;    % [SDR参数] 物理遮挡的有效波束宽度 (度)

%% 2. 提取卫星几何信息
% (用于 SDR 模式：判断哪些卫星位于发射机方向)
fprintf('   正在计算卫星几何分布...\n');

% 选取中间时刻计算几何 (近似认为手势期间卫星相对位置不变)
mid_idx = round(length(obs_data)/2);
try
    [rec_pos, ~, sat_states] = calculate_receiver_position(obs_data, nav_data, mid_idx);
    if isempty(rec_pos), error('接收机位置计算失败'); end
    [lat, lon, alt] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));
catch ME
    warning('无法计算几何分布，将降级为全星座处理: %s', ME.message);
    sat_states = struct();
end

% 构建卫星方位角映射表
sat_azimuths = containers.Map;
if ~isempty(fieldnames(sat_states))
    sat_list = fieldnames(sat_states);
    for k = 1:length(sat_list)
        sid = sat_list{k};
        sat_pos = sat_states.(sid).position;
        [e, n, ~] = ecef2enu(sat_pos(1)-rec_pos(1), sat_pos(2)-rec_pos(2), sat_pos(3)-rec_pos(3), lat, lon, alt);
        az = atan2d(e, n); if az < 0, az = az + 360; end
        sat_azimuths(sid) = az;
    end
end

%% 3. 执行攻击 (数据篡改)
num_samples = length(obs_spoofed);
flatten_count = 0;
keep_count = 0;

for i = 1:num_samples
    if isempty(obs_spoofed(i).data), continue; end
    sats = fieldnames(obs_spoofed(i).data);
    
    for k = 1:length(sats)
        sid = sats{k};
        
        % 智能查找 SNR 字段 (S1C, S2I, S1I 等)
        target_field = '';
        snr_struct = obs_spoofed(i).data.(sid).snr;
        if isfield(snr_struct, 'S1C'), target_field = 'S1C';
        elseif isfield(snr_struct, 'S2I'), target_field = 'S2I';
        elseif isfield(snr_struct, 'S1I'), target_field = 'S1I';
        elseif ~isempty(fieldnames(snr_struct))
             % 兜底: 取第一个字段
             fns = fieldnames(snr_struct); target_field = fns{1};
        end
        
        if isempty(target_field), continue; end
        original_val = obs_spoofed(i).data.(sid).snr.(target_field);
        if isnan(original_val) || original_val == 0, continue; end
        
        % --- 核心逻辑: 决定是保留波动还是抹平 ---
        should_flatten = false;
        
        if strcmp(ATTACK_TYPE, 'REPLAY')
            % [重放攻击]: 没有任何动作特征 -> 全部抹平
            should_flatten = true;
            
        elseif strcmp(ATTACK_TYPE, 'SDR')
            % [SDR攻击]: 只有 SDR 方向的卫星保留波动，其他方向抹平
            if isKey(sat_azimuths, sid)
                sat_az = sat_azimuths(sid);
                % 计算角度差 (处理 0/360 跨越)
                diff_az = abs(sat_az - SDR_AZIMUTH);
                if diff_az > 180, diff_az = 360 - diff_az; end
                
                if diff_az > SDR_BEAM_WIDTH / 2
                    % 卫星不在 SDR 物理波束内 -> 无法被手遮挡 -> 应该是平的
                    should_flatten = true;
                else
                    % 卫星在 SDR 方向 -> 能被手遮挡 -> 保留真实波动
                    should_flatten = false;
                end
            else
                % 几何未知的卫星，默认抹平以防漏网
                should_flatten = true; 
            end
        else
            error('未知的攻击类型: %s (支持 SDR / REPLAY)', ATTACK_TYPE);
        end
        
        % --- 执行抹平操作 ---
        if should_flatten
            % 使用平滑基线替代原始信号
            % 这里简单取 42 dB 作为基准，加微量白噪模拟接收机热噪
            noise = randn(1) * 0.5; 
            obs_spoofed(i).data.(sid).snr.(target_field) = 42 + noise; 
            flatten_count = flatten_count + 1;
        else
            keep_count = keep_count + 1;
        end
    end
end

%% 4. 结果摘要
fprintf('   仿真完成: 抹平点数 %d, 保留点数 %d\n', flatten_count, keep_count);
if strcmp(ATTACK_TYPE, 'SDR')
    fprintf('   [SDR Info] Azimuth: %.1f°, BeamWidth: %.1f°\n', SDR_AZIMUTH, SDR_BEAM_WIDTH);
end
fprintf('✅ 已返回欺骗后的数据结构 obs_spoofed。\n');
end