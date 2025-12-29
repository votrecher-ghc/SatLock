% =========================================================================
% test_snr_baseline_algorithm.m (函数版)
% 功能: 基准线趋势算法验证 (Algorithm Validation)
% 描述:
%   这是一个调试与验证工具，用于测试 "全局基准线 + 窗口趋势判定" 预处理算法
%   的有效性。它不进行手势定位，而是专注于展示 SNR 信号是如何被清洗的。
%   通过对比原始信号和处理后信号，您可以直观地调整噪声阈值和趋势窗口参数。
%
% [调用格式]:
%   [cn0_flat, valid_sats, t_grid] = test_snr_baseline_algorithm(obs_data, nav_data);
%
% [输入参数]:
%   1. obs_data (struct数组): 
%      原始观测数据。
%   2. nav_data (struct结构体): 
%      (可选) 导航星历数据。本函数主要处理 SNR，不涉及定位，
%      保留此参数是为了与其他分析函数接口保持一致。
%
% [返回值说明]:
%   1. cn0_flat (double矩阵 [Samples x Sats]): 
%      [处理后] 的信噪比矩阵。
%      - 每一列对应一颗卫星。
%      - 震荡已被归零到基准线，仅保留了显著的单边趋势 (手势动作)。
%   2. valid_sats (cell数组): 
%      卫星 ID 列表 (如 {'G01', 'C02', ...})，与矩阵的列一一对应。
%   3. t_grid (datetime列向量): 
%      插值对齐后的统一时间轴。
%
% [核心算法]:
%   1. Global Baseline: 利用众数 (Mode) 锁定环境底噪基准线。
%   2. Spike Rejection: 剔除短时 (<=5点) 的瞬时跳变毛刺。
%   3. Trend Window Check: 
%      - 若窗口内波动方向一致 (全正或全负) -> 判定为趋势，保留原始值。
%      - 若窗口内波动方向混合 (震荡) -> 判定为噪声，强制压平至基准线。
% =========================================================================

function [cn0_flat, valid_sats, t_grid] = test_snr_baseline_algorithm(obs_data, nav_data)

addpath(genpath('sky_plot')); 
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

fprintf('--> 启动基准线趋势算法验证 (Function版: Baseline Trend Algorithm)...\n');

%% 1. 参数设置
% --- [核心参数 - 请根据实际调整] ---
PARA.diff_lag_N         = 5;     % 趋势窗口 (N点): 用于判断窗口内的单调性/震荡
PARA.noise_cutoff_db    = 1;     % 噪声阈值 (dB): 与基准线偏差 <= 此值，视为平稳/噪声
PARA.spike_th_db        = 1;     % 毛刺幅度阈值 (dB): 偏离基准线 > 此值才检查是否为毛刺
PARA.spike_max_duration = 5;     % 毛刺最大持续点数 (1-5): 超过此宽度的跳变不视作毛刺

% --- 自动计算采样率 ---
raw_times = [obs_data.time];
mean_dt = mean(diff(raw_times));
PARA.sampling_rate = round(1 / seconds(mean_dt));
fprintf('   自动检测采样率: %d Hz\n', PARA.sampling_rate);

%% 2. 数据提取
fprintf('--> [计算] 正在提取全星座数据...\n');

all_sat_ids = {};
for i = 1:min(100, length(obs_data)), if ~isempty(obs_data(i).data), all_sat_ids = [all_sat_ids, fieldnames(obs_data(i).data)']; end; end
unique_sat_ids = unique(all_sat_ids);
valid_sats = {};
for i = 1:length(unique_sat_ids)
    sid = unique_sat_ids{i};
    if ismember(sid(1), ['G','C','E','J']), valid_sats{end+1} = sid; end
end

t_grid = (min(raw_times) : seconds(1/PARA.sampling_rate) : max(raw_times))';
t_grid_plot = t_grid + hours(8) - seconds(18);
num_samples = length(t_grid);
num_sats = length(valid_sats);

cn0_raw = NaN(num_samples, num_sats);
for s_idx = 1:num_sats
    sat_id = valid_sats{s_idx};
    target_snr_code = '';
    for k = 1:min(50, length(obs_data))
        if isfield(obs_data(k).data, sat_id) && isfield(obs_data(k).data.(sat_id), 'snr')
            fields = fieldnames(obs_data(k).data.(sat_id).snr);
            if ismember('S1C', fields), target_snr_code = 'S1C'; elseif ismember('S2I', fields), target_snr_code = 'S2I'; else, if ~isempty(fields), target_snr_code = fields{1}; end; end
            if ~isempty(target_snr_code), break; end
        end
    end
    if isempty(target_snr_code), continue; end
    s_times = []; s_cn0 = [];
    for k = 1:length(obs_data)
        if isfield(obs_data(k).data, sat_id) && isfield(obs_data(k).data.(sat_id).snr, target_snr_code)
            val = obs_data(k).data.(sat_id).snr.(target_snr_code);
            if ~isnan(val) && val > 10, s_times = [s_times; obs_data(k).time]; s_cn0 = [s_cn0; val]; end
        end
    end
    if length(s_times) > 20
        [u_times, u_idx] = unique(s_times);
        cn0_raw(:, s_idx) = interp1(u_times, s_cn0(u_idx), t_grid, 'linear', NaN);
    end
end

%% 3. 核心算法: 全局基准线 + 趋势清洗
fprintf('--> [计算] 执行新算法 (Baseline Mode + Window Check)...\n');

cn0_flat = cn0_raw;
N = PARA.diff_lag_N;
NoiseTh = PARA.noise_cutoff_db;
SpikeTh = PARA.spike_th_db;
SpikeDur = PARA.spike_max_duration;

for s = 1:num_sats
    raw_col = cn0_raw(:, s);
    
    % [Step 1] 计算全局基准线 (众数)
    % 只要数据不全是 NaN
    valid_data = raw_col(~isnan(raw_col));
    if isempty(valid_data), continue; end
    
    % 取整后计算众数，确保基准线是一个整数水位
    baseline = mode(round(valid_data)); 
    
    % 初始化输出列
    col = raw_col;
    
    % [Step 2] 逐点遍历与处理
    for t = 1:num_samples
        curr_val = raw_col(t);
        if isnan(curr_val), continue; end
        
        % -----------------------------------------------------
        % A. 毛刺检测 (Spike Check) - 针对基准线的瞬时大跳变
        % -----------------------------------------------------
        if abs(curr_val - baseline) > SpikeTh
            % 向后看 SpikeDur 个点，看是否回到基准线附近
            is_spike = false;
            for k = 1:SpikeDur
                if t + k > num_samples, break; end
                future_val = raw_col(t+k);
                if abs(future_val - baseline) <= NoiseTh
                    is_spike = true;
                    break;
                end
            end
            if is_spike
                col(t) = baseline; % 确认为毛刺，压平
                continue; % 处理下一点
            end
        end
        
        % -----------------------------------------------------
        % B. 窗口趋势判定 (Trend vs Oscillation)
        % -----------------------------------------------------
        % 提取当前窗口 [t, t+N]
        win_end = min(t + N - 1, num_samples);
        window_vals = raw_col(t : win_end);
        
        % 计算窗口内每个点相对于基准线的偏差
        diffs = window_vals - baseline;
        
        % 找出"显著偏差"的点 (绝对值 > 噪声阈值)
        sig_indices = abs(diffs) > NoiseTh;
        sig_diffs = diffs(sig_indices);
        
        if isempty(sig_diffs)
            % 情况1: 窗口内全是微小波动 -> 强制设为基准线
            col(t) = baseline;
        else
            % 情况2: 存在显著偏差，检查方向一致性
            all_pos = all(sig_diffs > 0); % 全在基准线上方
            all_neg = all(sig_diffs < 0); % 全在基准线下方
            
            if all_pos || all_neg
                % 单边趋势 (可能是 44 42 42 ... 或 45 39 39 ...) -> 保留原始值
                col(t) = curr_val;
            else
                % 混合方向 (震荡/穿过基准线) -> 强制设为基准线
                col(t) = baseline;
            end
        end
    end
    cn0_flat(:, s) = col;
end

fprintf('✅ 计算完成 (已返回处理后的 SNR 矩阵)。\n');

%% 4. 绘图 (默认配色)
fprintf('--> [绘图] 生成对比图...\n');
plot_cnt = 0;
for s = 1:num_sats
    if all(isnan(cn0_raw(:, s))), continue; end
    sat_id = valid_sats{s};
    plot_cnt = plot_cnt + 1;

    figure('Name', sprintf('[%s] Baseline Algorithm Check', sat_id), 'Position', [100+(plot_cnt*20), 100, 1000, 600], 'Color', 'w');

    % 子图1: 原始数据
    subplot(2, 1, 1);
    plot(t_grid_plot, cn0_raw(:, s), 'LineWidth', 1);
    
    % 在原图上画出计算出的基准线，方便核对
    base_val = mode(round(cn0_raw(~isnan(cn0_raw(:,s)), s)));
    yline(base_val, 'r--', 'Baseline', 'LineWidth', 1);
    
    title(sprintf('卫星 %s: 原始信号 (Global Baseline = %d)', sat_id, base_val));
    ylabel('SNR'); grid on; axis tight;
    ylim_val = [min(cn0_raw(:,s))-1, max(cn0_raw(:,s))+1];
    ylim(ylim_val);
    datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');

    % 子图2: 处理后数据
    subplot(2, 1, 2);
    plot(t_grid_plot, cn0_flat(:, s), 'LineWidth', 1.5);
    title(sprintf('新算法处理后 (N=%d, NoiseTh=%d) - 震荡归零，趋势保留', PARA.diff_lag_N, PARA.noise_cutoff_db));
    ylabel('SNR'); xlabel('Time'); grid on; axis tight;
    ylim(ylim_val);
    datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');

    if plot_cnt >= 50, fprintf('⚠️ 已显示前 10 颗卫星...\n'); break; end
end
fprintf('✅ 绘图完毕。\n');
end