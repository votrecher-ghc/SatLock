% =========================================================================
% test_snr_flattening.m (函数版)
% 功能: 验证 "智能削平" 算法 (v11.1 - 整数逻辑优化版)
% 描述:
%   这是一个用于验证信号预处理效果的调试工具。
%   它实现了一种"智能削平"逻辑，旨在保留显著的信号跳变(台阶)，同时
%   强力抹平微小的观测噪声和瞬间的宽毛刺。
%   该算法模拟了实时处理过程(逐点推进)，适用于流式应用。
%
% [调用格式]:
%   [cn0_flat, valid_sats, t_grid] = test_snr_flattening(obs_data, nav_data);
%
% [输入参数]:
%   1. obs_data (struct数组): 原始观测数据。
%   2. nav_data (struct结构体): (可选) 导航星历数据(仅作占位，保持接口统一)。
%
% [返回值说明]:
%   1. cn0_flat (double矩阵): 
%      [处理后] 的信噪比矩阵。震荡已被削平，呈现干净的台阶状。
%   2. valid_sats (cell数组): 
%      对应的卫星ID列表。
%   3. t_grid (datetime列向量): 
%      对应的时间轴。
%
% [核心算法]:
%   1. 宽毛刺剔除 (Wide Spike Removal): 
%      检测大幅跳变，若在 N 点内迅速回弹至原值，则视为干扰并剔除。
%   2. 智能稳态锁死 (Steady-State Locking):
%      若当前值与上一时刻值的偏差 <= 阈值，强制保持上一时刻值不变(削平)。
%   3. 趋势更新 (Trend Update):
%      仅当偏差 > 阈值且非毛刺时，才更新信号电平，形成台阶。
% =========================================================================

function [cn0_flat, valid_sats, t_grid] = test_snr_flattening(obs_data, nav_data)

addpath(genpath('sky_plot')); 
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));

fprintf('--> 启动智能削平验证 (Function版 v11.1: Integer Logic Fix)...\n');

%% 1. 参数设置

% --- [核心参数] ---
PARA.diff_lag_N         = 5;     % [趋势] 趋势窗口 (N点): 用于辅助判断长时趋势
PARA.noise_cutoff_db    = 1;     % [去噪] 噪声阈值 (dB): 波动幅度 <= 此值时视为噪声被削平
PARA.spike_th_db        = 2;     % [去噪] 毛刺幅度阈值 (dB): 跳变 > 此值才检查毛刺
PARA.spike_max_duration = 5;     % [去噪] 毛刺最大持续点数 (1-5): 允许剔除较宽的毛刺
PARA.sampling_rate      = 25;    % [采样] 默认采样率，后续会自动校正

% 自动计算采样率
raw_times = [obs_data.time];
mean_dt = mean(diff(raw_times));
if ~isnan(mean_dt)
    PARA.sampling_rate = round(1 / seconds(mean_dt));
end
fprintf('   当前采样率: %d Hz\n', PARA.sampling_rate);


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

%% 3. 核心算法: 智能削平 + 宽毛刺剔除
fprintf('--> [计算] 执行去噪 (Cutoff<=%d, SpikeDur=%d)...\n', PARA.noise_cutoff_db, PARA.spike_max_duration);

cn0_flat = cn0_raw;
N = PARA.diff_lag_N;
Th = PARA.noise_cutoff_db;
SpikeTh = PARA.spike_th_db;
MaxDur = PARA.spike_max_duration;

for s = 1:num_sats
    col = cn0_flat(:, s);
    raw_col = cn0_raw(:, s);

    last_val = NaN;
    skip_counter = 0; % 用于跳过被判定为毛刺的点

    for t = 1:num_samples
        curr_raw = raw_col(t);

        % 处理跳过逻辑
        if skip_counter > 0
            col(t) = last_val;
            skip_counter = skip_counter - 1;
            continue;
        end

        if isnan(curr_raw), last_val = NaN; continue; end
        if isnan(last_val), last_val = curr_raw; col(t) = last_val; continue; end

        % --- [宽毛刺检测] ---
        is_spike = false;
        spike_width = 0;

        delta_curr = curr_raw - last_val; 
        % 只有当发生大幅跳变时 (> SpikeTh) 才检查
        if abs(delta_curr) > SpikeTh
            % 向后看 MaxDur 个点
            for k = 1:MaxDur
                if t + k > num_samples, break; end
                future_val = raw_col(t+k);

                % 检查是否回弹 (回到 last_val 附近的 Th 范围内)
                if abs(future_val - last_val) <= Th
                    is_spike = true;
                    spike_width = k;
                    break;
                end
            end
        end

        if is_spike
            % 发现毛刺 -> 削平
            col(t) = last_val;
            skip_counter = spike_width - 1; 
            continue;
        end
        
        % --- [常规削平] ---
        if t > N, trend_diff = abs(curr_raw - raw_col(t-N)); else, trend_diff = 0; end

        if trend_diff > Th
            % 趋势显著 -> 更新
            col(t) = curr_raw;
            last_val = curr_raw;
        else
            % 趋势平稳 -> 检查瞬时偏差
            deviation = abs(curr_raw - last_val);

            % [修正] 使用 <= 符号，确保设为1时包含1
            if deviation <= Th
                % 偏差在阈值内 -> 强制削平，视为同一水平面
                col(t) = last_val;
            else
                % 偏差大 (且不是毛刺) -> 确认为台阶
                col(t) = curr_raw;
                last_val = curr_raw;
            end
        end
    end
    cn0_flat(:, s) = col;
end

fprintf('✅ 计算完成 (已返回处理后的 SNR 矩阵)。\n');

%% 4. 绘图
fprintf('--> [绘图] 生成对比图...\n');
plot_cnt = 0;
for s = 1:num_sats
    if all(isnan(cn0_raw(:, s))), continue; end
    sat_id = valid_sats{s};
    plot_cnt = plot_cnt + 1;

    figure('Name', sprintf('[%s] Smart Flattening v11.1', sat_id), 'Position', [100+(plot_cnt*20), 100, 1000, 600], 'Color', 'w');

    subplot(2, 1, 1);
    plot(t_grid_plot, cn0_raw(:, s), 'LineWidth', 1);
    title(sprintf('卫星 %s: 原始信号', sat_id));
    ylabel('SNR'); grid on; axis tight;
    ylim_val = [min(cn0_raw(:,s))-1, max(cn0_raw(:,s))+1];
    ylim(ylim_val);
    datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');

    subplot(2, 1, 2);
    plot(t_grid_plot, cn0_flat(:, s), 'LineWidth', 1.5);
    title(sprintf('处理效果 (Th<=%d, SpikeDur=%d) - 消除深坑毛刺', PARA.noise_cutoff_db, PARA.spike_max_duration));
    ylabel('SNR'); xlabel('Time'); grid on; axis tight;
    ylim(ylim_val);
    datetick('x', 'HH:MM:ss', 'keepticks', 'keeplimits');

    if plot_cnt >= 10, fprintf('⚠️ 已显示前 10 颗卫星...\n'); break; end
end
fprintf('✅ 绘图完毕。\n');
end