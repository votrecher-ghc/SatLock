% =========================================================================
% plot_sn.m (通用可视化工具 - 分图版 - v2 修复版)
% 功能: 循环遍历 obs 结构体，将每颗卫星的 SNR 绘制在独立的子图中
% 修复: 解决了 linkaxes 因空子图导致的数据类型冲突报错
%
% [调用格式]:
%   plot_sn(obs_data);
% =========================================================================

function plot_sn(obs_data)

% 获取变量名用于标题
input_var_name = inputname(1);
if isempty(input_var_name), input_var_name = 'OBS Data'; end

fprintf('--> [plot_sn] 正在为 "%s" 生成分星独立波形图...\n', input_var_name);

%% 1. 扫描卫星列表
if isempty(obs_data), warning('数据为空'); return; end

sat_map = containers.Map();
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

if num_sats == 0
    warning('未找到卫星数据。');
    return;
end

fprintf('    共发现 %d 颗卫星，准备生成子图矩阵...\n', num_sats);

%% 2. 计算布局并绘图
% 计算合适的网格行列数 (例如 30颗星 -> 5x6 布局)
grid_cols = ceil(sqrt(num_sats));
grid_rows = ceil(num_sats / grid_cols);

figure('Name', sprintf('Split View: %s', input_var_name), ...
       'Color', 'w', 'Position', [50, 50, 1400, 900]); % 大窗口

ax_handles = []; % 存储坐标轴句柄，用于同步缩放

for s = 1:num_sats
    sid = sat_list{s};
    
    % --- 数据提取 ---
    t_seq = [];
    snr_seq = [];
    
    for k = 1:length(obs_data)
        if isfield(obs_data(k).data, sid) && isfield(obs_data(k).data.(sid), 'snr')
            snr_struct = obs_data(k).data.(sid).snr;
            fds = fieldnames(snr_struct);
            if ~isempty(fds)
                val = snr_struct.(fds{1});
                if ~isnan(val) % 不过滤0值
                    t_seq = [t_seq; obs_data(k).time]; %#ok<AGROW>
                    snr_seq = [snr_seq; val];       %#ok<AGROW>
                end
            end
        end
    end
    
    % --- 创建子图 ---
    ax = subplot(grid_rows, grid_cols, s);
    
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    
    if ~isempty(t_seq)
        % 绘图: 使用蓝色细线
        plot(ax, t_seq, snr_seq, 'b.-', 'LineWidth', 1, 'MarkerSize', 5);
        
        % 简单的统计显示
        avg_val = mean(snr_seq(snr_seq > 0));
        text(ax, t_seq(1), max(snr_seq)*0.9, sprintf('Avg: %.1f', avg_val), ...
             'FontSize', 8, 'Color', [0.3 0.3 0.3]);
         
        % [关键修复] 只有当子图里真的有数据(时间轴)时，才加入到联动列表
        ax_handles = [ax_handles, ax]; %#ok<AGROW>
    else
        text(ax, 0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center');
    end
    
    % --- 子图修饰 ---
    title(ax, sid, 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'none');
    
    % 智能时间轴标签：只处理有数据的轴，且尽量只在最后显示
    if ~isempty(t_seq)
        datetick(ax, 'x', 'MM:SS', 'keepticks', 'keeplimits');
    end
    
    if s <= (grid_rows-1)*grid_cols
        set(ax, 'XTickLabel', []); % 隐藏非最后一行的标签，减少杂乱
    else
         xlabel(ax, 'Time');
    end
    
    axis(ax, 'tight');
end

% --- 交互优化 ---
% 链接所有子图的 X 轴 (缩放任意一个，全部同步)
if ~isempty(ax_handles)
    linkaxes(ax_handles, 'x');
end

sgtitle(sprintf('Split Satellite View: %s (Linked Zoom)', input_var_name), ...
        'Interpreter', 'none', 'FontSize', 14);

fprintf('✅ 绘图完成。所有有效子图 X 轴已联动。\n');

end