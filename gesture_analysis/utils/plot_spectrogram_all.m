% =========================================================================
% plot_spectrogram_all.m (全星座频谱诊断工具 - 修复版)
% 功能: 绘制时频图 (Spectrogram)
% 修复: 解决了 imagesc 不支持 datetime 格式导致的 "XData must be vectors" 错误
% =========================================================================

function plot_spectrogram_all(obs_data, fs)

if nargin < 2, fs = 25; end % 默认采样率

% 获取变量名用于标题
input_var_name = inputname(1);
if isempty(input_var_name), input_var_name = 'OBS Data'; end

fprintf('--> [Spectrogram] 正在为 "%s" 生成全星座频谱图 (FS=%dHz)...\n', input_var_name, fs);

%% 1. 扫描卫星列表
if isempty(obs_data), warning('数据为空'); return; end

sat_map = containers.Map();
all_times = [obs_data.time];
t_start = min(all_times);
t_end   = max(all_times);

% 建立统一时间轴
t_grid = (t_start : seconds(1/fs) : t_end)';
num_samples = length(t_grid);

for k = 1:length(obs_data)
    if isfield(obs_data(k), 'data')
        sats = fieldnames(obs_data(k).data);
        for i = 1:length(sats)
            sat_map(sats{i}) = 1;
        end
    end
end
sat_list = sort(keys(sat_map));
num_sats = length(sat_list);

if num_sats == 0, warning('未找到卫星数据。'); return; end

fprintf('    共发现 %d 颗卫星，准备生成频谱矩阵...\n', num_sats);

%% 2. 频谱参数设置
window_len = round(1.0 * fs); 
noverlap   = round(0.9 * window_len); 
nfft       = 512; % 频率分辨率

%% 3. 计算布局并绘图
grid_cols = ceil(sqrt(num_sats));
grid_rows = ceil(num_sats / grid_cols);

figure('Name', sprintf('Spectrogram View: %s', input_var_name), ...
       'Color', 'w', 'Position', [50, 50, 1400, 900]);

ax_handles = []; % 用于联动

for s = 1:num_sats
    sid = sat_list{s};
    
    % --- A. 数据提取 ---
    raw_t = []; raw_snr = [];
    
    for k = 1:length(obs_data)
        if isfield(obs_data(k).data, sid) && isfield(obs_data(k).data.(sid), 'snr')
            snr_struct = obs_data(k).data.(sid).snr;
            fds = fieldnames(snr_struct);
            if ~isempty(fds)
                val = snr_struct.(fds{1});
                if ~isnan(val) && val > 0
                    raw_t = [raw_t; obs_data(k).time]; %#ok<AGROW>
                    raw_snr = [raw_snr; val];       %#ok<AGROW>
                end
            end
        end
    end
    
    % --- 创建子图 ---
    ax = subplot(grid_rows, grid_cols, s);
    
    if length(raw_t) > 50 % 数据太少不画
        % 1. 均匀插值
        [u_t, u_idx] = unique(raw_t);
        snr_interp = interp1(u_t, raw_snr(u_idx), t_grid, 'linear', NaN);
        
        % 2. 填补 NaN
        valid_mask = ~isnan(snr_interp);
        if sum(valid_mask) > 50
             filled_snr = interp1(find(valid_mask), snr_interp(valid_mask), (1:num_samples)', 'nearest', 'extrap');
             
             % 3. 去直流 (Detrend)
             baseline = movmean(filled_snr, 2 * fs);
             signal_ac = filled_snr - baseline;
             
             % 4. 计算 Spectrogram
             [~, F, T, P] = spectrogram(signal_ac, window_len, noverlap, nfft, fs);
             
             % 5. 转换时间轴为数值 (datenum) [关键修复]
             T_plot = t_grid(1) + seconds(T);
             T_num = datenum(T_plot); % 将 datetime 转为 double
             
             % 6. 绘图 (dB谱)
             freq_idx = F <= 5.0; 
             
             % 使用数值类型的 T_num 进行绘图
             imagesc(ax, T_num, F(freq_idx), 10*log10(abs(P(freq_idx, :))));
             
             axis(ax, 'xy'); % 频率轴向上
             colormap(ax, 'jet');
             caxis(ax, [-30, 5]); % 固定色标
             
             % 7. 格式化 X 轴为时间格式
             datetick(ax, 'x', 'MM:SS', 'keepticks', 'keeplimits');
        else
            text(ax, 0.5, 0.5, 'Sparse Data', 'HorizontalAlignment', 'center');
        end
        ax_handles = [ax_handles, ax]; %#ok<AGROW>
    else
        text(ax, 0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center');
    end
    
    % --- 修饰 ---
    title(ax, sid, 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'none');
    
    if s > (grid_rows-1)*grid_cols
         xlabel(ax, 'Time');
    else
         set(ax, 'XTickLabel', []);
    end
    
    if mod(s-1, grid_cols) == 0
        ylabel(ax, 'Hz');
    else
        set(ax, 'YTickLabel', []);
    end
end

% --- 交互联动 ---
if ~isempty(ax_handles)
    linkaxes(ax_handles, 'x');
end

sgtitle(sprintf('Spectrogram View: %s (0-5Hz Energy Flow)', input_var_name), ...
        'Interpreter', 'none', 'FontSize', 14);

fprintf('✅ 频谱图生成完成。\n');

end