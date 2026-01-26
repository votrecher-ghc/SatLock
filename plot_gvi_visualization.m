% =========================================================================
% plot_gvi_visualization.m
% 功能: 通用 GVI (GNSS Volatility Index) 可视化工具
% 描述:
%   接收包含波动矩阵的结果结构体，计算总能量 GVI 并绘制曲线。
%   适用于 step1_res (模拟GVI), step1_res_shaped (整形GVI) 等任意阶段数据。
%   
%   [特性]:
%   1. 自动计算 GVI: Sum(Volatility Matrix)。
%   2. 自动高亮分段: 如果输入包含 segments，用红色加粗显示动作区间。
%   3. 阈值显示: 如果输入包含 PARA.gvi_threshold，绘制蓝色虚线。
%
% [调用格式]:
%   plot_gvi_visualization(step1_res);
%   plot_gvi_visualization(step1_res_shaped, 'Shaped Waveform');
%
% [输入参数]:
%   1. res_struct (struct): 
%      必须包含:
%        - .volatility_matrix (矩阵): 用于计算 GVI。
%        - .t_grid (datetime/double): 时间轴。
%      可选包含:
%        - .segments (struct array): 用于绘制红色高亮区间。
%        - .PARA.gvi_threshold (double): 用于绘制阈值线。
%   2. custom_title (string/char, 可选): 
%      自定义图表标题后缀。如果不传，默认使用变量名。
% =========================================================================

function plot_gvi_visualization(res_struct, custom_title)

    %% 1. 参数解析与检查
    if nargin < 2
        % 尝试获取传入变量的名称
        custom_title = inputname(1);
        if isempty(custom_title)
            custom_title = 'GVI Data';
        end
    end

    if ~isstruct(res_struct)
        error('❌ 输入必须是结构体 (struct)。');
    end

    if ~isfield(res_struct, 'volatility_matrix') || ~isfield(res_struct, 't_grid')
        error('❌ 输入结构体必须包含 .volatility_matrix 和 .t_grid 字段。');
    end

    %% 2. 数据提取与计算
    % 提取核心矩阵
    vol_mat = res_struct.volatility_matrix;
    t_axis  = res_struct.t_grid;
    
    % [核心公式] GVI = 所有卫星波动之和
    gvi_curve = sum(vol_mat, 2, 'omitnan');
    
    % 为了美观，对模拟信号做轻微平滑 (如果是方波则不影响)
    % 简单的 5点滑动平均
    gvi_curve_smooth = movmean(gvi_curve, 5);
    
    %% 3. 绘图执行
    figure('Name', sprintf('GVI Plot: %s', custom_title), ...
           'Color', 'w', 'Position', [100, 100, 1000, 400]);
       
    hold on; grid on; box on;
    
    % --- 3.1 绘制背景基准线 (黑色) ---
    % 这是整体的 GVI 曲线
    plot(t_axis, gvi_curve_smooth, 'k-', 'LineWidth', 1, 'DisplayName', 'Overall GVI');
    
    % --- 3.2 绘制高亮分段 (红色) ---
    if isfield(res_struct, 'segments') && ~isempty(res_struct.segments)
        segs = res_struct.segments;
        has_label = false; % 防止图例重复
        
        for i = 1:length(segs)
            % 获取该段的时间索引
            s_idx = segs(i).start_idx;
            e_idx = segs(i).end_idx;
            
            % 边界保护
            s_idx = max(1, s_idx);
            e_idx = min(length(t_axis), e_idx);
            
            if s_idx >= e_idx, continue; end
            
            % 提取该段数据
            seg_t = t_axis(s_idx:e_idx);
            seg_v = gvi_curve_smooth(s_idx:e_idx);
            
            % 绘图 (红色加粗)
            if ~has_label
                plot(seg_t, seg_v, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Detected Action');
                has_label = true;
            else
                plot(seg_t, seg_v, 'r-', 'LineWidth', 2.5, 'HandleVisibility', 'off');
            end
            
            % 标注 ID
            [max_v, max_i] = max(seg_v);
            text(seg_t(max_i), max_v, sprintf('#%d', segs(i).id), ...
                 'Color', 'r', 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');
        end
    end
    
    % --- 3.3 绘制阈值线 (蓝色虚线) ---
    if isfield(res_struct, 'PARA') && isfield(res_struct.PARA, 'gvi_threshold')
        th_val = res_struct.PARA.gvi_threshold;
        yline(th_val, 'b--', 'LineWidth', 1.5, 'Label', 'Threshold', 'DisplayName', 'Trigger Level');
    end
    
    %% 4. 图表修饰
    title(sprintf('GVI Energy Curve: %s', custom_title), 'Interpreter', 'none', 'FontSize', 12);
    xlabel('Time (HH:MM:SS)');
    ylabel('GVI Amplitude (Sum of Volatility)');
    
    % 智能时间轴格式化
    if isdatetime(t_axis)
        datetick('x', 'MM:SS', 'keepticks', 'keeplimits');
    end
    
    legend('show', 'Location', 'best');
    
    % 自动调整 Y 轴范围，留一点头部空间
    max_y = max(gvi_curve_smooth);
    if max_y > 0
        ylim([0, max_y * 1.2]);
    end
    
    fprintf('✅ 已绘制 GVI 曲线: %s\n', custom_title);
end