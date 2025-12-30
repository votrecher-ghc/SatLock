% =========================================================================
% waveform_reshaping.m (函数版)
% 功能: 中间层 - 信号二值化整形与可视化 (Matched with Step 1 v7.2)
% 描述:
%   该函数是 Step 1 (基准线清洗) 和 Step 2 (轨迹追踪) 的桥梁。
%   1. 读取 Step 1 算法清洗出的 cn0_clean_matrix。
%   2. 对 Step 1 计算的波动特征进行二值化和形态学整形 (Reshaping)。
%   3. 提供严格对齐的单星 4 子图对比。
%   4. [新增] 提供整体 GVI 趋势对比图，直观展示动作分段。
%
% [调用格式]:
%   [obs_waveform, step1_res_shaped] = waveform_reshaping(obs_raw, obs_clean, step1_res);
%
% [输入参数]:
%   1. obs_raw (struct): 
%      原始观测数据 (用于绘图 Subplot 1 对比)。
%   2. obs_clean (struct): 
%      (仅做占位或传递用，本函数主要读取 step1_res 中的清洗数据)。
%   3. step1_res (struct): 
%      Step 1 结果包 (含 cn0_clean_matrix, volatility_matrix, t_grid)。
%
% [返回值说明]:
%   1. obs_waveform (struct): 
%      返回 obs_clean 的副本。
%   2. step1_res_shaped (struct): 
%      更新后的结果包，其 .volatility_matrix 已被替换为整形后的方波数据。
% =========================================================================

function [obs_waveform, step1_res_shaped] = waveform_reshaping(obs_raw, obs_clean, step1_res)

fprintf('--> [Intermediate] 启动波形整形 (Waveform Reshaping v7.2)...\n');

%% 1. 初始化与解包
step1_res_shaped = step1_res; 
obs_waveform = obs_clean;    

% [关键] 从 step1_res 中提取核心矩阵
vol_mat   = step1_res.volatility_matrix;      
clean_mat = step1_res.cn0_clean_matrix; % 读取 Step 1 内部清洗结果      
t_grid    = step1_res.t_grid;                 
valid_sats = step1_res.valid_sats;
sampling_rate = step1_res.PARA.sampling_rate;
[num_samples, num_sats] = size(vol_mat);

% --- 整形参数 ---
SHAPE.binarize_th_db   = 2.0;   % [门限] 二值化阈值 (dB)
SHAPE.fixed_weight     = 10.0;  % [权重] 激活后的理想高度
SHAPE.gap_fill_sec     = 0.3;   % [闭运算] 填补裂缝 (秒)
SHAPE.gap_fill_pts     = round(SHAPE.gap_fill_sec * sampling_rate);
SHAPE.min_pulse_sec    = 0.12;  % [开运算] 剔除毛刺 (秒)
SHAPE.min_pulse_pts    = round(SHAPE.min_pulse_sec * sampling_rate);

fprintf('    整形参数: Th=%.1fdB, Fill=%.2fs, MinPulse=%.2fs\n', ...
    SHAPE.binarize_th_db, SHAPE.gap_fill_sec, SHAPE.min_pulse_sec);

%% 2. 核心整形逻辑 (Logic Only)
reshaped_mat = zeros(size(vol_mat));

for s = 1:num_sats
    raw_col_vol = vol_mat(:, s); 
    
    % --- 形态学处理 ---
    binary_col = raw_col_vol > SHAPE.binarize_th_db;
    
    % 1. 闭运算 (Fill Gaps)
    if SHAPE.gap_fill_pts > 0
        binary_col = func_fill_gaps(binary_col, SHAPE.gap_fill_pts);
    end
    
    % 2. 开运算 (Remove Spikes)
    if SHAPE.min_pulse_pts > 0
        binary_col = func_remove_spikes(binary_col, SHAPE.min_pulse_pts);
    end
    
    % 3. 赋权
    reshaped_col = double(binary_col) * SHAPE.fixed_weight;
    reshaped_mat(:, s) = reshaped_col;
end

% 更新结果包
step1_res_shaped.volatility_matrix = reshaped_mat;

%% 3. 单星绘图模块 (Individual Satellite Plots)
fprintf('--> [绘图] 生成单星对比图...\n');

% for s = 1:num_sats
%     sid = valid_sats{s};
%     
%     % 准备绘图数据
%     raw_col_vol  = vol_mat(:, s);
%     reshaped_col = reshaped_mat(:, s);
%     clean_curve  = clean_mat(:, s); 
%     
%     % 提取并插值 Raw SNR (用于 Subplot 1)
%     [raw_curve, has_raw] = extract_and_align_snr(obs_raw, sid, t_grid);
%     
%     % --- 绘图 (Strict 4 Subplots) ---
%     figure('Name', sprintf('Waveform Shaping: %s', sid), 'Position', [100+s*30, 50, 700, 900], 'Color', 'w');
%     
%     % [Subplot 1] 原始数据
%     ax1 = subplot(4, 1, 1); grid(ax1, 'on'); hold(ax1, 'on');
%     plot(ax1, t_grid, nan(size(t_grid)), 'HandleVisibility', 'off'); % 强制初始化
%     if has_raw
%         plot(ax1, t_grid, raw_curve, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
%     else
%         text(ax1, t_grid(1), 0, 'No Raw Data');
%     end
%     title(ax1, sprintf('[%s] 1. Original Raw SNR (Interpolated)', sid));
%     ylabel(ax1, 'SNR (dBHz)'); axis(ax1, 'tight');
%     
%     % [Subplot 2] 清洗数据
%     ax2 = subplot(4, 1, 2); grid(ax2, 'on'); hold(ax2, 'on');
%     plot(ax2, t_grid, nan(size(t_grid)), 'HandleVisibility', 'off');
%     plot(ax2, t_grid, clean_curve, 'b-', 'LineWidth', 1.5);
%     title(ax2, sprintf('[%s] 2. Baseline Cleaned SNR (Step 1 v7.2)', sid));
%     ylabel(ax2, 'SNR (dBHz)'); axis(ax2, 'tight');
%     
%     % [Subplot 3] 模拟波动
%     ax3 = subplot(4, 1, 3); grid(ax3, 'on'); hold(ax3, 'on');
%     plot(ax3, t_grid, raw_col_vol, 'k-', 'LineWidth', 1);
%     yline(ax3, SHAPE.binarize_th_db, 'r--', 'Threshold');
%     title(ax3, sprintf('[%s] 3. Analog Volatility (Feature)', sid));
%     ylabel(ax3, 'Vol (dB)'); axis(ax3, 'tight');
%     
%     % [Subplot 4] 整形波形
%     ax4 = subplot(4, 1, 4); grid(ax4, 'on'); hold(ax4, 'on');
%     area(ax4, t_grid, reshaped_col, 'FaceColor', [0.0 0.5 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'b');
%     title(ax4, sprintf('[%s] 4. Final Reshaped Waveform (Ideal)', sid));
%     ylabel(ax4, 'Weight'); ylim(ax4, [-1, SHAPE.fixed_weight + 2]); axis(ax4, 'tight');
%     
%     linkaxes([ax1, ax2, ax3, ax4], 'x');
%     drawnow;
% end

%% 4. [新增] 整体 GVI 绘图模块 (Overall GVI Plot)
fprintf('--> [绘图] 生成 Overall GVI Analysis 总图...\n');

% 计算累加 GVI
gvi_analog   = sum(vol_mat, 2, 'omitnan');       % 原始波动总和
gvi_reshaped = sum(reshaped_mat, 2, 'omitnan');  % 整形方波总和

figure('Name', 'Overall GVI Analysis', 'Position', [150, 150, 1000, 500], 'Color', 'w');
hold on; grid on;

% 画模拟 GVI (黑色曲线)
plot(t_grid, gvi_analog, 'k-', 'LineWidth', 1, 'DisplayName', 'Analog GVI (Raw Energy)');

% 画整形 GVI (蓝色填充)
area(t_grid, gvi_reshaped, 'FaceColor', [0.0 0.5 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'b', 'DisplayName', 'Reshaped GVI (Action Blocks)');

% 画阈值线 (红色虚线)
yline(step1_res.PARA.gvi_threshold, 'r--', 'Activation Threshold', 'LabelHorizontalAlignment', 'left', 'LineWidth', 1.5);

title('Overall Gesture Activity Analysis (Analog vs Reshaped)');
ylabel('GVI Index');
xlabel('Time');
legend('Location', 'best');
axis tight;
datetick('x', 'MM:SS', 'keepticks', 'keeplimits');

fprintf('✅ 整形完成。已生成单星图及总体 GVI 分析图。\n');

end

%% ================= Helper Functions =================

function [aligned_data, success] = extract_and_align_snr(obs_struct, sat_id, t_grid)
    aligned_data = NaN(size(t_grid)); success = false;
    if isempty(obs_struct), return; end
    
    target_code = '';
    for k = 1:min(50, length(obs_struct))
        if isfield(obs_struct(k).data, sat_id) && isfield(obs_struct(k).data.(sat_id), 'snr')
             fds = fieldnames(obs_struct(k).data.(sat_id).snr);
             if ~isempty(fds), target_code = fds{1}; break; end
        end
    end
    if isempty(target_code), return; end
    
    s_times = []; s_vals = [];
    for k = 1:length(obs_struct)
        if isfield(obs_struct(k).data, sat_id) && isfield(obs_struct(k).data.(sat_id).snr, target_code)
            val = obs_struct(k).data.(sat_id).snr.(target_code);
            if ~isnan(val) && val > 0
                s_times = [s_times; obs_struct(k).time];
                s_vals  = [s_vals; val];
            end
        end
    end
    
    if length(s_times) > 5
        [u_times, u_idx] = unique(s_times);
        aligned_data = interp1(u_times, s_vals(u_idx), t_grid, 'linear', NaN);
        success = true;
    end
end

function out_bw = func_fill_gaps(in_bw, max_gap_pts)
    out_bw = in_bw; padded = [1; in_bw; 1]; 
    gap_starts = find(diff(padded) == -1); gap_ends = find(diff(padded) == 1) - 1; 
    for k = 1:length(gap_starts)
        len = gap_ends(k) - gap_starts(k) + 1;
        if len <= max_gap_pts, out_bw(gap_starts(k)-1 : gap_ends(k)-1) = true; end
    end
end

function out_bw = func_remove_spikes(in_bw, min_pulse_pts)
    out_bw = in_bw; padded = [0; in_bw; 0];
    pulse_starts = find(diff(padded) == 1); pulse_ends = find(diff(padded) == -1) - 1; 
    for k = 1:length(pulse_starts)
        len = pulse_ends(k) - pulse_starts(k) + 1;
        if len < min_pulse_pts, out_bw(pulse_starts(k)-1 : pulse_ends(k)-1) = false; end
    end
end