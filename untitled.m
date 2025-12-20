% =========================================================================
% plot_gvi_principle_v5_final.m
% 功能: GVI 原理示意图 (V5 - 终极修正版)
% 修正: 
%   1. [Baseline]: 彻底修复"下凹"问题。Baseline 像桥一样直直地跨过深坑。
%   2. [Raw Signal]: 模拟真实的信号陡降 (Deep Dip)。
%   3. [逻辑]: 清晰展示 GVI = |桥面 - 坑底|。
% =========================================================================

clear; clc; close all;

% 1. 设置
N = 500;                % 采样点数
x = 1:N;                % 横坐标

% 2. 生成"完美的"基准线 (The Bridge)
% 它是缓慢变化的，或者是直的，绝不会因为手势而下凹
% 模拟卫星仰角带来的极缓慢趋势
true_baseline = 45 + 0.5 * sin(2 * pi * x / 1000); 

% 3. 生成手势动作 (The Dip)
center_x = 250;         
sigma = 8;              % 极窄、极陡
% 这是一个纯粹的干扰项
gesture_dip = -15 * exp(-((x - center_x).^2) / (2 * sigma^2)); 

% 4. 生成原始信号 (Raw = Baseline + Dip + Noise)
noise = 0.2 * randn(size(x)); 
raw_cn0 = true_baseline + gesture_dip + noise;

% 5. 算法模拟
% 注意：为了绘图展示原理，我们直接画"理想基线" (true_baseline)
% 这样能最清晰地告诉评委：我们要减去的是背景，保留的是动作
volatility = abs(raw_cn0 - true_baseline);

% GVI 聚合放大
gvi = volatility * 4; 
gvi_display = movmean(gvi, 5); 

% ================= [开始绘图] =================
figure('Name', 'GVI Principle Final', 'Position', [100, 100, 800, 700]);

% --- 子图 1: 理想的分离 (Ideal Separation) ---
subplot(3, 1, 1);
hold on; grid on; box on;

% 1. 先画基线 (它是稳定的背景)
% 使用较粗的虚线，表示这是参考系
plot(x, true_baseline, '--', 'LineWidth', 2, 'DisplayName', 'Baseline (Trend)'); 

% 2. 后画原始信号
% 我们可以清楚地看到信号掉下去了，但基线还在上面
plot(x, raw_cn0, '-', 'LineWidth', 1.2, 'DisplayName', 'Raw Signal'); 

ylabel('C/N_0 (dB-Hz)');
title('(a) Signal Phenomenon: Fast Drop vs. Stable Baseline');
legend('Location', 'southeast'); 
ylim([28, 50]); 
xlim([0, N]);

% 加个箭头强调差值
text(center_x, 35, 'Difference', 'HorizontalAlignment','center', 'VerticalAlignment','bottom');
annotation('doublearrow', [0.52 0.52], [0.75 0.85]); % 手动调整位置的大概箭头

% --- 子图 2: 波动提取 ---
subplot(3, 1, 2);
hold on; grid on; box on;
plot(x, volatility, 'LineWidth', 1.2);
area_obj = area(x, volatility);
area_obj.FaceAlpha = 0.2;
area_obj.EdgeColor = 'none';

ylabel('Volatility (dB)');
title('(b) Volatility Extraction |Raw - Baseline|');
xlim([0, N]);

% --- 子图 3: GVI 指数 ---
subplot(3, 1, 3);
hold on; grid on; box on;
plot(x, gvi_display, 'LineWidth', 2);
yline(8, '--', 'LineWidth', 1.5, 'Label', 'Threshold');

xlabel('Sampling Points');
ylabel('GVI Amplitude');
title('(c) GVI Index: Energy Aggregation');
xlim([0, N]);

% 自动调整
set(gcf, 'PaperPositionMode', 'auto');



% %% =================================================================
% %  Paper_Figure_Outlier_Removal_v3_3.m
% %  功能: 信号预处理对比图 - 深色阴影背景，移除绿色描述文字
% % =================================================================
% clc; clear; close all;
% 
% % --- 1. 信号构建 ---
% N = 1000;                       % 采样点
% x_axis = 1:N;                   
% baseline = 48;                  % 基准线 (dB-Hz)
% 
% % A. 背景噪声
% rng(42); 
% white_noise = 0.01 * randn(1, N); 
% 
% % B. 真实手势波动 (Gesture - 4dB Deep Fade)
% gesture_center = 600;
% gesture_width = 30;
% % 主波谷
% main_lobe = -4.0 * exp(-((x_axis - gesture_center).^2)/(2*gesture_width^2));
% % 旁瓣
% side_lobes = 0.8 * exp(-((x_axis - (gesture_center-80)).^2)/(2*25^2)) ...
%            + 0.8 * exp(-((x_axis - (gesture_center+80)).^2)/(2*25^2));
%        
% gesture_sig = main_lobe + side_lobes;
% 
% % C. 脉冲噪声 (Spike)
% spike_sig = zeros(1, N);
% spike_loc = 250;
% spike_sig(spike_loc) = -2;     
% spike_sig(spike_loc-1) = -0.2;   
% spike_sig(spike_loc+1) = -0.2;
% 
% % --- 合成信号 ---
% raw_signal = baseline + gesture_sig + spike_sig + white_noise;
% processed_signal = baseline + gesture_sig + white_noise; 
% 
% % --- [修改] 定义阴影区域参数 ---
% shade_start = 450; 
% shade_end = 750;
% % [修改] 颜色加深: 0.92 -> 0.85
% shade_color = [0.85, 0.85, 0.85]; 
% y_lim_range = [43.5, 49.0];       
% 
% % --- 2. 绘图设置 ---
% fig = figure('Color', 'w', 'Position', [100, 100, 900, 600]);
% t = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
% 
% % ---------------- [子图 1: 处理前 Raw Data] ----------------
% nexttile;
% hold on; 
% 
% % 1. 先画背景阴影 (Darker)
% x_shade = [shade_start, shade_end, shade_end, shade_start];
% y_shade = [y_lim_range(1), y_lim_range(1), y_lim_range(2), y_lim_range(2)];
% fill(x_shade, y_shade, shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
% 
% % 2. [修改] 仅保留阴影标注，位置稍微居中
% text(gesture_center, 48.7, 'Gesture Active Region', ...
%     'FontName', 'Times New Roman', 'FontSize', 11, 'Color', 'k', ...
%     'HorizontalAlignment', 'center', 'FontWeight', 'bold');
% 
% % 3. 画基准线
% yline(baseline, 'r--', 'LineWidth', 1.0, 'Alpha', 0.6); 
% 
% % 4. 画原始信号
% plot(x_axis, raw_signal, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.8); 
% 
% % 标注: 噪声 (红色标注保留，绿色已移除)
% annotation('textarrow', [0.30, 0.26], [0.80, 0.72], 'String', {'Impulsive Noise', '(Small)'}, ...
%            'FontName', 'Times New Roman', 'FontSize', 11, 'Color', 'r', 'TextColor', 'r', 'LineWidth', 1.2);
% 
% hold off;
% 
% % 样式设置
% title('(a) Raw Observation Data (Signal Mixed with Noise)', 'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');
% ylabel('C/N_0 (dB-Hz)', 'FontName', 'Times New Roman', 'FontSize', 12);
% grid on; box on;
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2, 'Layer', 'top'); 
% xlim([1, N]); 
% ylim(y_lim_range); 
% 
% % ---------------- [子图 2: 处理后 Processed Data] ----------------
% nexttile;
% hold on;
% 
% % 1. 先画背景阴影
% fill(x_shade, y_shade, shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
% 
% % 2. 画处理后的信号
% plot(x_axis, processed_signal, 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth', 1.8); 
% 
% % 3. 绘制原来的 spike 位置的虚影
% plot(spike_loc-5:spike_loc+5, raw_signal(spike_loc-5:spike_loc+5), 'Color', [0.6 0.6 0.6], 'LineStyle', ':', 'LineWidth', 1.2);
% 
% hold off;
% 
% % 标注: 仅保留噪声去除说明
% text(spike_loc+20, 46.5, '\leftarrow Noise Removed', 'FontName', 'Times New Roman', 'FontSize', 11, 'Color', 'r', 'FontWeight', 'bold');
% % [已移除] Feature Preserved 的绿色标注
% 
% % 样式设置
% title('(b) Processed Data (Clean Signal)', 'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');
% xlabel('Sampling Points', 'FontName', 'Times New Roman', 'FontSize', 12);
% ylabel('C/N_0 (dB-Hz)', 'FontName', 'Times New Roman', 'FontSize', 12);
% grid on; box on;
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2, 'Layer', 'top'); 
% xlim([1, N]); 
% ylim(y_lim_range); 
% 
% % --- 导出 ---
% % exportgraphics(fig, 'Fig_Signal_Process_Final.png', 'Resolution', 300);