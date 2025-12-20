% =========================================================================
% simulate_satellite_visibility.m
% 功能: 模拟绘制一天中可见卫星数量的变化曲线
% 核心逻辑:
%   1. 生成 0-24 小时的时间轴。
%   2. 使用正弦波叠加随机噪声模拟卫星数量的缓变特性。
%   3. 将数值约束在 [25, 35] 区间内并取整。
%   4. 绘制曲线 (使用默认蓝色)。
% =========================================================================

%% 1. 参数设置
clear; clc;

% 模拟参数
sim_duration_hours = 24;        % 模拟时长 (小时)
time_step_min = 5;              % 时间步长 (分钟)
min_sats = 25;                  % 最小卫星数 (基准下界)
max_sats = 35;                  % 最大卫星数 (基准上界)

%% 2. 生成模拟数据
% 创建时间轴 (单位: 小时)
t_hours = 0 : (time_step_min/60) : sim_duration_hours;
num_points = length(t_hours);

% 生成基础波动曲线 (模拟星座几何变化的周期性)
% 使用两个不同周期的正弦波叠加，制造非单调的波动感
base_wave = sin(t_hours / 3.8) + 0.6 * cos(t_hours / 1.5);

% 归一化并映射到目标区间 [min_sats, max_sats]
% 1. 将波形归一化到 [0, 1]
wave_norm = (base_wave - min(base_wave)) / (max(base_wave) - min(base_wave));

% 2. 映射到 [25, 35]
sat_counts_smooth = min_sats + wave_norm * (max_sats - min_sats);

% 3. 添加微小的随机扰动 (模拟个别卫星的进出)
noise = (rand(1, num_points) - 0.5) * 1.5; 
sat_counts_noisy = sat_counts_smooth + noise;

% 4. 取整 (卫星数量必须是整数)
sat_counts_final = round(sat_counts_noisy);

% 5. 强制边界约束 (防止噪声导致超出范围)
sat_counts_final(sat_counts_final < min_sats) = min_sats;
sat_counts_final(sat_counts_final > max_sats) = max_sats;

%% 3. 绘图 (Visualization)
f = figure('Name', 'Simulated Satellite Visibility', 'Color', 'w');
ax = axes('Parent', f);
hold(ax, 'on'); 
grid(ax, 'on');

% --- 绘制曲线 ---
% 严格遵守要求: 不指定颜色，使用Matlab默认蓝色
plot(t_hours, sat_counts_final, 'LineWidth', 1.5);

% --- 坐标轴设置 ---
xlabel('Time (Hours)');
ylabel('Number of Visible Satellites');
xlim([0, 24]);
ylim([20, 40]); % Y轴范围稍微留白，好看一点
set(ax, 'XTick', 0:4:24); % X轴每4小时一个刻度

% --- 统计信息与标题 ---
curr_max = max(sat_counts_final);
curr_min = min(sat_counts_final);
curr_avg = mean(sat_counts_final);

title_str = {
    '\bf Number of Visible Satellites (24H Simulation)',
    sprintf('\\rm Range: [%d, %d] | Average: %.1f', curr_min, curr_max, curr_avg)
};
title(title_str);

% 添加图例 (可选)
legend('Total Visible Satellites (GPS+BDS+GAL+GLO)', 'Location', 'best');

fprintf('✅ 模拟完成。平均卫星数: %.1f\n', curr_avg);