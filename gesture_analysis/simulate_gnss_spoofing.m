% =========================================================================
% simulate_gnss_spoofing_v2.m
% 功能: 基于"物理特征缺失"的攻击仿真
% 核心思想: 
%   1. SDR攻击: 只有SDR物理方向有波动，其余方向因物理缺失而呈平滑直线。
%   2. 重放攻击: 数据中完全没有当前手势的物理特征，全星座均为平滑直线。
% =========================================================================

%% 1. 基础设置与数据备份
clearvars -except obs_data nav_data;
if ~exist('obs_data', 'var'), error('请先加载正常数据 obs_data!'); end

% 备份原始数据 (以便恢复)
if ~exist('obs_data_backup', 'var')
    obs_data_backup = obs_data;
    fprintf('📦 已备份原始数据至 obs_data_backup\n');
else
    obs_data = obs_data_backup; % 每次运行前先恢复
    fprintf('🔄 已从备份恢复原始数据\n');
end

% --- 攻击配置 ---
ATTACK_TYPE = 'REPLAY';  % 可选: 'SDR' (单源物理遮挡) 或 'REPLAY' (无动作重放)
SDR_AZIMUTH = 120;    % 假设 SDR 发射机位于方位角 120 度
SDR_BEAM_WIDTH = 30;  % SDR 物理遮挡的波束宽度 (度)

fprintf('⚠️  正在构建 [%s] 攻击仿真环境...\n', ATTACK_TYPE);

%% 2. 提取卫星几何信息 (用于判断哪些卫星在 SDR 方向)
% 我们需要知道每颗卫星的方位角，以便决定是保留还是抹平
fprintf('   正在计算卫星几何分布...\n');
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('sky_plot'));

% 选取一个中间时刻计算几何 (近似认为手势期间卫星位置不变)
mid_idx = round(length(obs_data)/2);
[rec_pos, ~, sat_states] = calculate_receiver_position(obs_data, nav_data, mid_idx);
if isempty(rec_pos), error('无法计算接收机位置，仿真终止'); end
[lat, lon, alt] = ecef2geodetic(rec_pos(1), rec_pos(2), rec_pos(3));

sat_azimuths = containers.Map;
sat_list = fieldnames(sat_states);
for k = 1:length(sat_list)
    sid = sat_list{k};
    sat_pos = sat_states.(sid).position;
    [e, n, ~] = ecef2enu(sat_pos(1)-rec_pos(1), sat_pos(2)-rec_pos(2), sat_pos(3)-rec_pos(3), lat, lon, alt);
    az = atan2d(e, n); if az < 0, az = az + 360; end
    sat_azimuths(sid) = az;
end

%% 3. 执行攻击 (数据篡改)

num_samples = length(obs_data);
flatten_count = 0;
keep_count = 0;

for i = 1:num_samples
    if isempty(obs_data(i).data), continue; end
    sats = fieldnames(obs_data(i).data);
    
    for k = 1:length(sats)
        sid = sats{k};
        
        % 获取原始 C/N0 (假设 S1C 或 S2I)
        target_field = '';
        if isfield(obs_data(i).data.(sid).snr, 'S1C'), target_field = 'S1C';
        elseif isfield(obs_data(i).data.(sid).snr, 'S2I'), target_field = 'S2I';
        elseif isfield(obs_data(i).data.(sid).snr, 'S1I'), target_field = 'S1I';
        end
        
        if isempty(target_field), continue; end
        original_val = obs_data(i).data.(sid).snr.(target_field);
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
                    % 卫星不在 SDR 物理波束内 -> 手挡不到 -> 应该是平的
                    should_flatten = true;
                else
                    % 卫星在 SDR 方向 -> 手挡住了 SDR -> 保留波动 (或注入波动)
                    should_flatten = false;
                end
            else
                should_flatten = true; % 未知位置的卫星默认抹平
            end
        end
        
        % --- 执行抹平操作 ---
        if should_flatten
            % 使用一个带噪声的常数来模拟"环境声曲线"
            % 这里简单取 40 dB 作为基准，加一点点高斯白噪
            % (更高级的做法是取该卫星前段时间的均值)
            noise = randn(1) * 0.2; 
            obs_data(i).data.(sid).snr.(target_field) = 42 + noise; 
            flatten_count = flatten_count + 1;
        else
            keep_count = keep_count + 1;
        end
    end
end

%% 4. 结果摘要
fprintf('\n=== 仿真结果 [%s] ===\n', ATTACK_TYPE);
if strcmp(ATTACK_TYPE, 'SDR')
    fprintf('   SDR物理方位: %.1f° (波束宽 %.1f°)\n', SDR_AZIMUTH, SDR_BEAM_WIDTH);
    fprintf('   -> 位于 SDR 波束内的卫星保留了真实波动 (模拟遮挡SDR)。\n');
    fprintf('   -> 其他方向卫星已被替换为平滑环境噪声 (模拟物理缺失)。\n');
else
    fprintf('   -> 全星座数据已被替换为平滑环境噪声 (模拟无动作重放)。\n');
end
fprintf('   (受影响采样点: 抹平 %d 个, 保留 %d 个)\n', flatten_count, keep_count);

fprintf('\n👉 下一步: 请运行 run_gesture_analysis_boundary_trackV3.m 查看防御效果。\n');
fprintf('   预期结果: \n');
fprintf('   1. SDR模式: 可能检测到少量 Hit 卫星，但因数量不足或分布过于集中，无法解算出有效轨迹。\n');
fprintf('   2. Replay模式: GVI 能量极低，直接提示 "未检测到有效手势片段"。\n');