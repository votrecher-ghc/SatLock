% =========================================================================
% plot_single_sat_snr.m
% 功能: 绘制指定卫星的指定（或自动检测）SNR 信号的时序波形图
%
% [调用格式]:
%   plot_single_sat_snr(obs_data, target_sat);            % 自动选择第一个SNR频点
%   plot_single_sat_snr(obs_data, target_sat, obs_code);  % 指定SNR频点 (如 'S1C')
%
% [输入参数]:
%   obs_data   - struct array, 解析后的观测数据结构体 (由 parse_rinex_obs 生成)
%   target_sat - string, 目标卫星编号 (例如 'G01', 'C02')
%   obs_code   - (可选) string, 指定的SNR观测码 (例如 'S1C', 'S2I')。
%                如果省略或为空，函数将自动探测该卫星包含的第一个SNR类型。
%
% [依赖]:
%   extract_satellite_observation.m
% =========================================================================

function plot_single_sat_snr(obs_data, target_sat, obs_code)

    % --- 1. 参数检查与自动推断观测码 ---
    % 如果用户没有传入 obs_code，或者传入为空，则尝试自动查找
    if nargin < 3 || isempty(obs_code)
        detected_code = '';
        
        % 遍历数据寻找目标卫星的 SNR 字段定义
        for k = 1:length(obs_data)
            if isfield(obs_data(k).data, target_sat) && ...
               isfield(obs_data(k).data.(target_sat), 'snr')
                
                snr_struct = obs_data(k).data.(target_sat).snr;
                field_list = fieldnames(snr_struct);
                
                if ~isempty(field_list)
                    detected_code = field_list{1}; % 默认取第一个
                    break; % 找到一个即可停止
                end
            end
        end
        
        if isempty(detected_code)
            warning('在数据中未找到卫星 %s 的任何 SNR 观测值，无法绘图。', target_sat);
            return;
        else
            obs_code = detected_code;
            fprintf('--> [提示] 未指定观测码，自动选择卫星 %s 的第一个 SNR 频点: [%s]\n', target_sat, obs_code);
        end
    end

    % --- 2. 调用工具函数提取数据 ---
    % 复用 extract_satellite_observation 进行数据抽取
    [time_vec, snr_vec] = extract_satellite_observation(obs_data, target_sat, obs_code);

    % --- 3. 检查数据有效性 ---
    if isempty(time_vec)
        fprintf('--> [警告] 卫星 %s 的观测码 %s 没有提取到有效数据。\n', target_sat, obs_code);
        return;
    end

    % --- 4. 绘图 ---
    figure_name = sprintf('Satellite SNR: %s - %s', target_sat, obs_code);
    figure('Name', figure_name, 'Color', 'w');
    
    plot(time_vec, snr_vec, 'b.-', 'LineWidth', 1, 'MarkerSize', 8);
    
    % --- 5. 图形修饰 ---
    grid on;
    box on;
    
    % 设置标题和轴标签
    title(sprintf('Satellite: %s | Signal: %s', target_sat, obs_code), ...
          'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none');
    xlabel('Time (UTC)');
    ylabel('SNR (dB-Hz)');
    
    % --- 6. 智能纵轴范围控制 (Smart Y-Limits) ---
    % 目的: 避免波形“顶满”纵轴，同时兼容普通SNR和整形后的方波数据
    
    max_val = max(snr_vec);
    min_val = min(snr_vec);
    
    if max_val < 25
        % [情况A] 判断为数字整形信号 (通常是 0 或 10)
        % 设置固定范围，保证 0 和 10 都有充足留白
        ylim([-2, 14]); 
    else
        % [情况B] 判断为普通模拟 SNR 信号 (通常 > 30 dBHz)
        % 动态计算留白 (Margin)
        data_span = max_val - min_val;
        
        % 至少留 2dB，或者数据波动幅度的 20%
        margin = max(2.0, data_span * 0.2); 
        
        % 应用范围
        ylim([min_val - margin, max_val + margin]);
    end
    
    % 时间轴格式化
    try
        datetick('x', 'HH:MM:SS', 'keepticks', 'keeplimits');
    catch
        % 忽略 datetick 报错
    end
    
    % 调整横坐标范围紧凑
    xlim([min(time_vec) max(time_vec)]);

    fprintf('✅ 绘图完成: %s (%d 个历元)\n', figure_name, length(time_vec));

end