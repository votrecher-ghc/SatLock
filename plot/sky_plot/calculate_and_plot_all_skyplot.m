% ============== calculate_and_plot_all_skyplot.m (全卫星显示+双列图例版) ==============
function [all_sats_data] = calculate_and_plot_all_skyplot(obs_data, nav_data)
    % 功能: 
    %   1. 进行单点定位以获取接收机概略位置。
    %   2. 绘制所有可见卫星 (G/C/R/E/J) 的天空图。
    %   3. 图例分为两列显示。
    
    fprintf('--> 开始执行整合分析流程...\n');
    
    % --- 1. 对所有历元进行定位解算 ---
    num_epochs = length(obs_data);
    positioning_results = NaN(num_epochs, 7); 
    
    h_wait = waitbar(0, '正在对所有历元进行定位解算，请稍候...');
    for epoch_idx = 1:length(obs_data)
        waitbar(epoch_idx / num_epochs, h_wait);
        try
            [ecef_pos, receiver_clk_err, ~, lat, lon, alt] = ...
                calculate_receiver_position(obs_data, nav_data, epoch_idx);
            if ~any(isnan(ecef_pos))
                positioning_results(epoch_idx, :) = [lat, lon, alt, ecef_pos', receiver_clk_err];
            end
        catch ME
            % 静默失败
        end
    end
    close(h_wait);
    
    positioning_results = positioning_results(~isnan(positioning_results(:,1)), :);
    success_count = size(positioning_results, 1);
    
    % --- 2. 打印关键的平均结果并绘图 ---
    if success_count > 1
        mean_ecef_pos = mean(positioning_results(:, 4:6), 1);
        [mean_lat, mean_lon, mean_alt] = ecef2geodetic(mean_ecef_pos(1), mean_ecef_pos(2), mean_ecef_pos(3));
        mean_clk_err = mean(positioning_results(:, 7));
        fprintf('\n================== [最终分析结果] ==================\n');
        fprintf('成功解算率: %d / %d 个历元\n', success_count, num_epochs);
        fprintf('------------------ 平均定位结果 ------------------\n');
        fprintf('  > 平均纬度: %.6f°\n', mean_lat);
        fprintf('  > 平均经度: %.6f°\n', mean_lon);
        fprintf('  > 平均高程: %.3f m\n', mean_alt);
        fprintf('  > 平均钟差: %.9f 秒\n', mean_clk_err);
        fprintf('====================================================\n\n');

        figure('Name', '接收机钟差变化图');
        plot(positioning_results(:,7), '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
        title(sprintf('接收机钟差变化 (%d 个点)', success_count));
        xlabel('成功定位的历元序号');
        ylabel('钟差 (秒)');
        grid on;
        
        receiver_pos_ecef = mean_ecef_pos';
        lat = mean_lat;
        lon = mean_lon;
    else
        error('‼️ 没有足够多的成功定位点，无法继续进行天空图计算。');
    end

    % --- 3. 计算并绘制天空图 ---
    fprintf('--> 正在计算并绘制所有卫星天空图...\n');
    
    all_sat_ids = [];
    for i = 1:length(obs_data)
        all_sat_ids = [all_sat_ids; fieldnames(obs_data(i).data)];
    end
    unique_sat_ids = unique(all_sat_ids);
    
    % [修改 1] 移除过滤，保留所有解析到的卫星 (G, C, R, E, J 等)
    % unique_sat_ids = unique_sat_ids(startsWith(unique_sat_ids, 'G') | startsWith(unique_sat_ids, 'C'));
    
    all_sats_data = struct();
    plotted_sats_count = 0;
    
    for sat_idx = 1:length(unique_sat_ids)
        sid = unique_sat_ids{sat_idx};
        az = []; el = [];
        
        for epoch_i = 1:length(obs_data)
            epoch = obs_data(epoch_i);
            if ~isfield(epoch.data, sid), continue; end
            sat_obs = epoch.data.(sid);
            
            % [修改 2] 动态寻找伪距信号，适配不同卫星系统
            pr_code = '';
            
            % 优先查找常用码
            if isfield(sat_obs.pseudorange, 'C2I'), pr_code = 'C2I';      % BDS B1I
            elseif isfield(sat_obs.pseudorange, 'C1C'), pr_code = 'C1C';  % GPS/GAL/QZSS L1
            elseif isfield(sat_obs.pseudorange, 'C1P'), pr_code = 'C1P';  % GLONASS
            end
            
            % 兜底策略：如果以上都没有，取第一个存在的字段
            if isempty(pr_code)
                available_codes = fieldnames(sat_obs.pseudorange);
                if ~isempty(available_codes)
                    pr_code = available_codes{1};
                end
            end
            
            % 检查数据有效性
            if isempty(pr_code) || ~isfield(sat_obs.pseudorange, pr_code) || isnan(sat_obs.pseudorange.(pr_code))
                continue; 
            end
            
            try
                [sat_pos, ~, ~] = calculate_satellite_state(epoch.time, sat_obs.pseudorange.(pr_code), sid, nav_data);
                vec = sat_pos - receiver_pos_ecef;
                [e, n, u] = ecef2enu(vec(1), vec(2), vec(3), lat, lon, 0);
                az_i = atan2d(e, n); if az_i < 0, az_i = az_i + 360; end
                el_i = asind(u / norm([e, n, u]));
                if el_i >= 0, az(end+1) = az_i; el(end+1) = el_i; end
            catch
                continue;
            end
        end
        
        if ~isempty(az)
            plotted_sats_count = plotted_sats_count + 1;
            all_sats_data(plotted_sats_count).ID = sid;
            all_sats_data(plotted_sats_count).azimuth = az;
            all_sats_data(plotted_sats_count).elevation = el;
        end
    end
    
    % --- 绘图部分 ---
    if plotted_sats_count > 0
        figure('Name', '所有卫星天空图 (带起止点标记)');
        pax = polaraxes;
        
        % 保持原有的极坐标设置
        set(pax, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise', 'RAxisLocation', 90, 'RDir', 'reverse');
        rlim(pax, [0 90]);
        rticks(pax, [0 30 60 90]);
        rticklabels(pax, {'90', '60', '30', '0'});
        pax.RAxis.Label.String = '仰角 (°)';
        title(sprintf('所有可见卫星天空轨迹图 (共 %d 颗)', plotted_sats_count));
        hold(pax, 'on');
        
        plot_handles = gobjects(plotted_sats_count, 1);
        legend_entries = cell(1, plotted_sats_count);
        
        % 保持原有的配色方案 lines()
        colors = lines(plotted_sats_count);
        
        for i = 1:plotted_sats_count
            theta = deg2rad(all_sats_data(i).azimuth);
            % 保持原有的半径映射 (直接使用仰角配合 reverse 轴)
            r = all_sats_data(i).elevation;
            
            % 绘制线条
            plot_handles(i) = polarplot(pax, theta, r, '-', 'LineWidth', 2.0, 'Color', colors(i,:));
            % 绘制起点 (空心圈)
            polarplot(pax, theta(1), r(1), 'o', 'MarkerEdgeColor', colors(i,:), 'MarkerSize', 8, 'LineWidth', 1.5);
            % 绘制终点 (实心圈)
            polarplot(pax, theta(end), r(end), 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:), 'MarkerSize', 8);
            % 绘制文字标签
            text(pax, theta(end), r(end), ['  ' all_sats_data(i).ID], 'FontSize', 10, 'FontWeight', 'bold', 'Color', colors(i,:));
            
            legend_entries{i} = all_sats_data(i).ID;
        end
        
        % [修改 3] 图例分成两列 (NumColumns = 2)
        legend(plot_handles, legend_entries, 'Location', 'eastoutside', 'NumColumns', 2);
        
        hold(pax, 'off');
    end
    fprintf('✅ 所有分析和绘图任务完成！\n\n');
end















% % ============== calculate_and_plot_all_skyplot.m (全系统通用版) ==============
% function [all_sats_data] = calculate_and_plot_all_skyplot(obs_data, nav_data)
%     % 功能:
%     %   - 计算并绘制所有已解析系统的卫星天空图 (G, R, C, E, J, S)。
%     %   - 自动适配伪距信号类型，不再局限于 C1C/C2I。
%     
%     fprintf('--> 开始执行整合分析流程 (全系统通用版)...\n');
%     
%     % --- 1. 对所有历元进行定位解算 ---
%     num_epochs = length(obs_data);
%     positioning_results = NaN(num_epochs, 7); 
%     
%     h_wait = waitbar(0, '正在对所有历元进行定位解算，请稍候...');
%     for epoch_idx = 1:length(obs_data)
%         waitbar(epoch_idx / num_epochs, h_wait);
%         try
%             % 调用定位函数
%             [ecef_pos, receiver_clk_err, ~, lat, lon, alt] = ...
%                 calculate_receiver_position(obs_data, nav_data, epoch_idx);
%             
%             if ~any(isnan(ecef_pos))
%                 positioning_results(epoch_idx, :) = [lat, lon, alt, ecef_pos', receiver_clk_err];
%             end
%         catch ME
%             % 静默失败，继续下一个历元
%         end
%     end
%     close(h_wait);
%     
%     % 统计成功率
%     positioning_results = positioning_results(~isnan(positioning_results(:,1)), :);
%     success_count = size(positioning_results, 1);
%     
%     % --- 2. 打印关键的平均结果并绘图 ---
%     if success_count > 1
%         mean_ecef_pos = mean(positioning_results(:, 4:6), 1);
%         [mean_lat, mean_lon, mean_alt] = ecef2geodetic(mean_ecef_pos(1), mean_ecef_pos(2), mean_ecef_pos(3));
%         mean_clk_err = mean(positioning_results(:, 7));
% 
%         fprintf('\n================== [最终分析结果] ==================\n');
%         fprintf('成功解算率: %d / %d 个历元\n', success_count, num_epochs);
%         fprintf('------------------ 平均定位结果 ------------------\n');
%         fprintf('  > 平均纬度: %.6f°\n', mean_lat);
%         fprintf('  > 平均经度: %.6f°\n', mean_lon);
%         fprintf('  > 平均高程: %.3f m\n', mean_alt);
%         fprintf('  > 平均钟差: %.9f 秒\n', mean_clk_err);
%         fprintf('====================================================\n\n');
% 
%         figure('Name', '接收机钟差变化图');
%         plot(positioning_results(:,7), '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
%         title(sprintf('接收机钟差变化 (%d 个点)', success_count));
%         xlabel('成功定位的历元序号');
%         ylabel('钟差 (秒)');
%         grid on;
%         
%         receiver_pos_ecef = mean_ecef_pos';
%         lat = mean_lat;
%         lon = mean_lon;
%     else
%         fprintf('‼️ 定位解算失败率过高 (成功历元 < 2)，无法计算天空图。\n');
%         error('‼️ 没有足够多的成功定位点，无法继续进行天空图计算。');
%     end
% 
%     % --- 3. 计算并绘制天空图 (全系统) ---
%     fprintf('--> 正在计算并绘制所有卫星天空图...\n');
%     
%     all_sat_ids = [];
%     for i = 1:length(obs_data)
%         all_sat_ids = [all_sat_ids; fieldnames(obs_data(i).data)];
%     end
%     unique_sat_ids = unique(all_sat_ids);
%     
%     % [修改点1] 放宽筛选：包含 G(GPS), R(GLONASS), C(BDS), E(Galileo), J(QZSS), S(SBAS)
%     % 只要是以这些字母开头的卫星ID都保留
%     mask = startsWith(unique_sat_ids, 'G') | ...
%            startsWith(unique_sat_ids, 'R') | ...
%            startsWith(unique_sat_ids, 'C') | ...
%            startsWith(unique_sat_ids, 'E') | ...
%            startsWith(unique_sat_ids, 'J') | ...
%            startsWith(unique_sat_ids, 'S');
%     unique_sat_ids = unique_sat_ids(mask);
%                                     
%     fprintf('    在数据中共找到 %d 颗不同的卫星 (涵盖所有解析系统)。\n', length(unique_sat_ids));
%     
%     all_sats_data = struct();
%     plotted_sats_count = 0;
%     
%     for sat_idx = 1:length(unique_sat_ids)
%         sid = unique_sat_ids{sat_idx};
%         az = []; el = [];
%         
%         for epoch_i = 1:length(obs_data)
%             epoch = obs_data(epoch_i);
%             if ~isfield(epoch.data, sid), continue; end
%             sat_obs = epoch.data.(sid);
%             
%             % [修改点2] 智能选择伪距信号
%             % 不再硬编码，而是查找该卫星有哪些伪距数据可用
%             available_codes = fieldnames(sat_obs.pseudorange);
%             if isempty(available_codes), continue; end
%             
%             pr_code = '';
%             sys = sid(1);
%             
%             % 优先级策略：如果有常用频点优先用，没有则用第一个可用的
%             if sys == 'C' && isfield(sat_obs.pseudorange, 'C2I')
%                 pr_code = 'C2I'; % 北斗优先 B1I (C2I)
%             elseif isfield(sat_obs.pseudorange, 'C1C')
%                 pr_code = 'C1C'; % GPS/GAL/QZSS/GLONASS 优先 L1 C/A
%             elseif isfield(sat_obs.pseudorange, 'C1P')
%                 pr_code = 'C1P'; % GLONASS 有时用 C1P
%             else
%                 % 兜底策略：如果以上首选都没有，直接取第一个存在的信号
%                 pr_code = available_codes{1};
%             end
%             
%             % 再次检查数值有效性
%             if isempty(pr_code) || isnan(sat_obs.pseudorange.(pr_code))
%                 continue; 
%             end
%             
%             try
%                 % 计算卫星位置
%                 [sat_pos, ~, ~] = calculate_satellite_state(epoch.time, sat_obs.pseudorange.(pr_code), sid, nav_data);
%                 
%                 % 转换为站心坐标系 (Az/El)
%                 vec = sat_pos - receiver_pos_ecef;
%                 [e, n, u] = ecef2enu(vec(1), vec(2), vec(3), lat, lon, 0);
%                 
%                 az_i = atan2d(e, n); if az_i < 0, az_i = az_i + 360; end
%                 el_i = asind(u / norm([e, n, u]));
%                 
%                 if el_i >= 0, az(end+1) = az_i; el(end+1) = el_i; end
%             catch 
%                 % 如果 calculate_satellite_state 报错（例如缺少该系统星历），跳过
%                 continue;
%             end
%         end
%         
%         if ~isempty(az)
%             plotted_sats_count = plotted_sats_count + 1;
%             all_sats_data(plotted_sats_count).ID = sid;
%             all_sats_data(plotted_sats_count).azimuth = az;
%             all_sats_data(plotted_sats_count).elevation = el;
%         end
%     end
%     
%     % --- 4. 绘制天空图 ---
%     if plotted_sats_count > 0
%         figure('Name', '全系统卫星天空图');
%         pax = polaraxes;
%         set(pax, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise', 'RAxisLocation', 90, 'RDir', 'reverse');
%         rlim(pax, [0 90]);
%         rticks(pax, [0 30 60 90]);
%         rticklabels(pax, {'90', '60', '30', '0'});
%         pax.RAxis.Label.String = '仰角 (°)';
%         title(sprintf('全系统卫星天空轨迹图 (共 %d 颗)', plotted_sats_count));
%         hold(pax, 'on');
%         
%         plot_handles = gobjects(plotted_sats_count, 1);
%         legend_entries = cell(1, plotted_sats_count);
%         
%         for i = 1:plotted_sats_count
%             % [修改点3] 扩展颜色定义，支持更多系统
%             sys_col = [0.5 0.5 0.5]; % 默认灰色
%             id_str = all_sats_data(i).ID;
%             
%             if startsWith(id_str, 'G'), sys_col = [0 0.4470 0.7410];     % GPS: 蓝
%             elseif startsWith(id_str, 'C'), sys_col = [0.8500 0.3250 0.0980]; % BDS: 橙
%             elseif startsWith(id_str, 'E'), sys_col = [0.4660 0.6740 0.1880]; % GAL: 绿
%             elseif startsWith(id_str, 'J'), sys_col = [0.4940 0.1840 0.5560]; % QZSS: 紫
%             elseif startsWith(id_str, 'R'), sys_col = [1 0 1];            % GLO: 品红
%             elseif startsWith(id_str, 'S'), sys_col = [0 1 1];            % SBAS: 青
%             end
%             
%             theta = deg2rad(all_sats_data(i).azimuth);
%             % r = 90 - all_sats_data(i).elevation; % 原版投影方式
%             r = all_sats_data(i).elevation;        % 标准极坐标投影 (0在中心, 90在边缘? 还是反过来?)
%             % 这里的极坐标轴设置了 'RDir', 'reverse'，所以 r=90(天顶)在圆心，r=0(地平)在边缘
%             % 直接传入仰角即可。
%             
%             % 绘制轨迹
%             plot_handles(i) = polarplot(pax, theta, r, '-', 'LineWidth', 2.0, 'Color', sys_col);
%             % 绘制起点 (空心圈)
%             polarplot(pax, theta(1), r(1), 'o', 'MarkerEdgeColor', sys_col, 'MarkerSize', 6, 'LineWidth', 1);
%             % 绘制终点 (实心圈)
%             polarplot(pax, theta(end), r(end), 'o', 'MarkerFaceColor', sys_col, 'MarkerEdgeColor', sys_col, 'MarkerSize', 6);
%             % 标注ID
%             text(pax, theta(end), r(end), ['  ' id_str], 'FontSize', 8, 'FontWeight', 'bold', 'Color', sys_col);
%             legend_entries{i} = id_str;
%         end
%         
%         % 图例处理：卫星太多时图例会遮挡，只显示前20个或者按系统显示图例更好
%         % 这里保留显示所有，但建议位置放旁边
%         legend(plot_handles, legend_entries, 'Location', 'eastoutside', 'NumColumns', 2);
%         hold(pax, 'off');
%     else
%         warning('没有卫星可用于绘制天空图。');
%     end
%     fprintf('✅ 所有分析和绘图任务完成！\n\n');
% end







% % ============== calculate_and_plot_all_skyplot.m (新增钟差绘图功能) ==============
% function [all_sats_data] = calculate_and_plot_all_skyplot(obs_data, nav_data)
% % 该版本在优化版基础上，新增了接收机钟差变化的绘图和统计。
% 
% fprintf('--> 开始执行整合分析流程 (优化版+钟差绘图)...\n');
% 
% % --- 1. [优化] 遍历所有历元，进行定位解算并存储所有结果 ---
% num_epochs = length(obs_data);
% % 结果矩阵列: [lat, lon, alt, x, y, z, clock_bias_sec]
% positioning_results = NaN(num_epochs, 7); 
% success_count = 0;
% 
% fprintf('\n=========================================================\n');
% fprintf('           步骤 1: 对所有 %d 个历元进行定位解算            \n', num_epochs);
% fprintf('=========================================================\n');
% h_wait = waitbar(0, '正在进行定位解算...');
% 
% for epoch_idx = 1:length(obs_data)
%     waitbar(epoch_idx / num_epochs, h_wait, sprintf('定位解算中... (%d / %d)', epoch_idx, num_epochs));
%     try
%         % [修改] 接收钟差输出 receiver_clk_err
%         [ecef_pos, receiver_clk_err, ~, lat, lon, alt] = ...
%             calculate_receiver_position(obs_data, nav_data, epoch_idx);
%         
%         if ~any(isnan(ecef_pos))
%             % [修改] 将钟差结果也存储起来
%             positioning_results(epoch_idx, :) = [lat, lon, alt, ecef_pos', receiver_clk_err];
%         end
%     catch ME
%         fprintf('历元 %d 定位失败: %s\n', epoch_idx, ME.message);
%     end
% end
% close(h_wait);
% 
% % 清理掉失败的历元 (NaN行)
% positioning_results = positioning_results(~isnan(positioning_results(:,1)), :);
% success_count = size(positioning_results, 1);
% 
% fprintf('\n>>> 定位解算完成！成功率: %d / %d 个历元 <<<\n\n', success_count, num_epochs);
% 
% % --- 2. [优化] 打印平均结果并绘制所有结果图 ---
% if success_count > 1
%     % 计算并打印平均值
%     mean_ecef_pos = mean(positioning_results(:, 4:6), 1);
%     [mean_lat, mean_lon, mean_alt] = ecef2geodetic(mean_ecef_pos(1), mean_ecef_pos(2), mean_ecef_pos(3));
%     mean_clk_err = mean(positioning_results(:, 7)); % [新增] 计算平均钟差
% 
%     fprintf('****************** 平均定位结果 ******************\n');
%     fprintf('  > 平均纬度: %.6f°\n', mean_lat);
%     fprintf('  > 平均经度: %.6f°\n', mean_lon);
%     fprintf('  > 平均高程: %.3f m\n', mean_alt);
%     fprintf('  > 平均钟差: %.9f 秒\n', mean_clk_err); % [新增] 打印平均钟差
%     fprintf('**************************************************\n\n');
% 
%     % 绘制轨迹图 (使用大地坐标)
%     figure('Name', '接收机定位轨迹图');
%     geoplot(positioning_results(:,1), positioning_results(:,2), 'b.-', 'MarkerSize', 10);
%     geobasemap satellite;
%     title(sprintf('接收机定位轨迹图 (%d 个点)', success_count));
%     grid on;
% 
%     % 绘制高程图 (使用大地坐标)
%     figure('Name', '接收机高程变化图');
%     plot(positioning_results(:,3), 'm-o', 'LineWidth', 1.5, 'MarkerSize', 4);
%     title(sprintf('接收机高程变化 (%d 个点)', success_count));
%     xlabel('成功定位的历元序号');
%     ylabel('高程 (米)');
%     grid on;
%     
%     % ====================【新增：绘制钟差变化图】====================
%     figure('Name', '接收机钟差变化图');
%     plot(positioning_results(:,7), 'g-o', 'LineWidth', 1.5, 'MarkerSize', 4);
%     title(sprintf('接收机钟差变化 (%d 个点)', success_count));
%     xlabel('成功定位的历元序号');
%     ylabel('钟差 (秒)');
%     grid on;
%     % =================================================================
%     
%     % [优化] 将平均ECEF和平均大地坐标作为天空图计算的参考
%     receiver_pos_ecef = mean_ecef_pos';
%     lat = mean_lat;
%     lon = mean_lon;
% 
% else
%     error('‼️ 没有足够多的成功定位点，无法继续进行天空图计算。');
% end
% 
% % --- 3. [优化] 使用平均位置，为所有卫星计算天空图数据 ---
% fprintf('\n=========================================================\n');
% fprintf('           步骤 2: 计算所有卫星的天空图轨迹            \n');
% fprintf('=========================================================\n');
% % (这部分代码保持不变，省略以保持简洁)
% % ...
% all_sat_ids = [];
% for i = 1:length(obs_data)
%     all_sat_ids = [all_sat_ids; fieldnames(obs_data(i).data)];
% end
% unique_sat_ids = unique(all_sat_ids);
% unique_sat_ids = unique_sat_ids(startsWith(unique_sat_ids, 'G') | startsWith(unique_sat_ids, 'C'));
% fprintf('    在数据中共找到 %d 颗不同的GPS/北斗卫星。\n', length(unique_sat_ids));
% fprintf('    使用平均参考位置 (ECEF): [%.2f, %.2f, %.2f]\n', receiver_pos_ecef);
% all_sats_data = struct();
% plotted_sats_count = 0;
% for sat_idx = 1:length(unique_sat_ids)
%     sid = unique_sat_ids{sat_idx};
%     az = []; el = [];
%     for epoch_i = 1:length(obs_data)
%         epoch = obs_data(epoch_i);
%         if ~isfield(epoch.data, sid), continue; end
%         
%         sat_obs = epoch.data.(sid);
%         pr_code = '';
%         if sid(1) == 'G' && isfield(sat_obs.pseudorange, 'C1C'), pr_code = 'C1C'; end
%         if sid(1) == 'C' && isfield(sat_obs.pseudorange, 'C2I'), pr_code = 'C2I'; end
%         if isempty(pr_code) || ~isfield(sat_obs.pseudorange, pr_code) || isnan(sat_obs.pseudorange.(pr_code)), continue; end
%         
%         try
%             [sat_pos, ~, ~] = calculate_satellite_state(epoch.time, sat_obs.pseudorange.(pr_code), sid, nav_data);
%             vec = sat_pos - receiver_pos_ecef;
%             [e, n, u] = ecef2enu(vec(1), vec(2), vec(3), lat, lon, 0);
%             
%             az_i = atan2d(e, n); if az_i < 0, az_i = az_i + 360; end
%             el_i = asind(u / norm([e, n, u]));
%             
%             if el_i >= 0
%                 az(end+1) = az_i;
%                 el(end+1) = el_i;
%             end
%         catch
%             continue;
%         end
%     end
%     if ~isempty(az)
%         plotted_sats_count = plotted_sats_count + 1;
%         all_sats_data(plotted_sats_count).ID = sid;
%         all_sats_data(plotted_sats_count).azimuth = az;
%         all_sats_data(plotted_sats_count).elevation = el;
%     end
% end
% 
% % --- 4. 绘制天空图 ---
% if plotted_sats_count == 0
%     error('‼️ 没有卫星可用于绘制天空图。');
% end
% fprintf('--> 数据计算完成，开始绘制最终天空图...\n');
% figure('Name', '所有卫星天空图 (带起止点标记)');
% pax = polaraxes;
% set(pax, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise', 'RAxisLocation', 90, 'RDir', 'reverse');
% rlim(pax, [0 90]);
% rticks(pax, [0 30 60 90]);
% rticklabels(pax, {'90', '60', '30', '0'});
% pax.RAxis.Label.String = '仰角 (°)';
% title('所有可见GPS/北斗卫星天空轨迹图');
% hold(pax, 'on');
% plot_handles = gobjects(plotted_sats_count, 1);
% legend_entries = cell(1, plotted_sats_count);
% colors = lines(plotted_sats_count);
% for i = 1:plotted_sats_count
%     theta = deg2rad(all_sats_data(i).azimuth);
%     r = 90 - all_sats_data(i).elevation;
%     plot_handles(i) = polarplot(pax, theta, r, '-', 'LineWidth', 2.0, 'Color', colors(i,:));
%     polarplot(pax, theta(1), r(1), 'o', 'MarkerEdgeColor', colors(i,:), 'MarkerSize', 8, 'LineWidth', 1.5);
%     polarplot(pax, theta(end), r(end), 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:), 'MarkerSize', 8);
%     text(pax, theta(end), r(end), ['  ' all_sats_data(i).ID], ...
%         'FontSize', 10, 'FontWeight', 'bold', 'Color', colors(i,:));
%     legend_entries{i} = all_sats_data(i).ID;
% end
% legend(plot_handles, legend_entries, 'Location', 'eastoutside');
% hold(pax, 'off');
% fprintf('✅ 所有分析和绘图任务完成！\n\n');
% 
% end




% 
% 
% % ============== calculate_and_plot_all_skyplot.m (简洁输出版) ==============
% function [all_sats_data] = calculate_and_plot_all_skyplot(obs_data, nav_data)
%     % 该版本只在最后输出关键的汇总信息和图形。
%     
%     fprintf('--> 开始执行整合分析流程...\n');
%     
%     % --- 1. 对所有历元进行定位解算 ---
%     num_epochs = length(obs_data);
%     positioning_results = NaN(num_epochs, 7); 
%     
%     h_wait = waitbar(0, '正在对所有历元进行定位解算，请稍候...');
%     for epoch_idx = 1:length(obs_data)
%         waitbar(epoch_idx / num_epochs, h_wait);
%         try
%             [ecef_pos, receiver_clk_err, ~, lat, lon, alt] = ...
%                 calculate_receiver_position(obs_data, nav_data, epoch_idx);
%             if ~any(isnan(ecef_pos))
%                 positioning_results(epoch_idx, :) = [lat, lon, alt, ecef_pos', receiver_clk_err];
%             end
%         catch ME
%             % 静默失败
%         end
%     end
%     close(h_wait);
%     
%     positioning_results = positioning_results(~isnan(positioning_results(:,1)), :);
%     success_count = size(positioning_results, 1);
%     
%     % --- 2. 打印关键的平均结果并绘图 ---
%     if success_count > 1
%         mean_ecef_pos = mean(positioning_results(:, 4:6), 1);
%         [mean_lat, mean_lon, mean_alt] = ecef2geodetic(mean_ecef_pos(1), mean_ecef_pos(2), mean_ecef_pos(3));
%         mean_clk_err = mean(positioning_results(:, 7));
% 
%         fprintf('\n================== [最终分析结果] ==================\n');
%         fprintf('成功解算率: %d / %d 个历元\n', success_count, num_epochs);
%         fprintf('------------------ 平均定位结果 ------------------\n');
%         fprintf('  > 平均纬度: %.6f°\n', mean_lat);
%         fprintf('  > 平均经度: %.6f°\n', mean_lon);
%         fprintf('  > 平均高程: %.3f m\n', mean_alt);
%         fprintf('  > 平均钟差: %.9f 秒\n', mean_clk_err);
%         fprintf('====================================================\n\n');
% 
% %         figure('Name', '接收机定位轨迹图');
% %         geoplot(positioning_results(:,1), positioning_results(:,2), 'b.-', 'MarkerSize', 10);
% %         geobasemap satellite;
% %         title(sprintf('接收机定位轨迹图 (%d 个点)', success_count));
% %         grid on;
% 
% %         figure('Name', '接收机高程变化图');
% %         plot(positioning_results(:,3), 'm-o', 'LineWidth', 1.5, 'MarkerSize', 4);
% %         title(sprintf('接收机高程变化 (%d 个点)', success_count));
% %         xlabel('成功定位的历元序号');
% %         ylabel('高程 (米)');
% %         grid on;
% 
%         figure('Name', '接收机钟差变化图');
%         plot(positioning_results(:,7), '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
%         title(sprintf('接收机钟差变化 (%d 个点)', success_count));
%         xlabel('成功定位的历元序号');
%         ylabel('钟差 (秒)');
%         grid on;
%         
%         receiver_pos_ecef = mean_ecef_pos';
%         lat = mean_lat;
%         lon = mean_lon;
%     else
%         error('‼️ 没有足够多的成功定位点，无法继续进行天空图计算。');
%     end
% 
%     % --- 3. 计算并绘制天空图 ---
%     fprintf('--> 正在计算并绘制卫星天空图...\n');
%     % (此部分代码与您提供的版本完全相同，省略以保持简洁)
%     all_sat_ids = [];
%     for i = 1:length(obs_data)
%         all_sat_ids = [all_sat_ids; fieldnames(obs_data(i).data)];
%     end
%     unique_sat_ids = unique(all_sat_ids);
%     unique_sat_ids = unique_sat_ids(startsWith(unique_sat_ids, 'G') | startsWith(unique_sat_ids, 'C'));
%     all_sats_data = struct();
%     plotted_sats_count = 0;
%     for sat_idx = 1:length(unique_sat_ids)
%         sid = unique_sat_ids{sat_idx};
%         az = []; el = [];
%         for epoch_i = 1:length(obs_data)
%             epoch = obs_data(epoch_i);
%             if ~isfield(epoch.data, sid), continue; end
%             sat_obs = epoch.data.(sid);
%             pr_code = '';
%             if sid(1) == 'G' && isfield(sat_obs.pseudorange, 'C1C'), pr_code = 'C1C'; end
%             if sid(1) == 'C' && isfield(sat_obs.pseudorange, 'C2I'), pr_code = 'C2I'; end
%             if isempty(pr_code) || ~isfield(sat_obs.pseudorange, pr_code) || isnan(sat_obs.pseudorange.(pr_code)), continue; end
%             try
%                 [sat_pos, ~, ~] = calculate_satellite_state(epoch.time, sat_obs.pseudorange.(pr_code), sid, nav_data);
%                 vec = sat_pos - receiver_pos_ecef;
%                 [e, n, u] = ecef2enu(vec(1), vec(2), vec(3), lat, lon, 0);
%                 az_i = atan2d(e, n); if az_i < 0, az_i = az_i + 360; end
%                 el_i = asind(u / norm([e, n, u]));
%                 if el_i >= 0, az(end+1) = az_i; el(end+1) = el_i; end
%             catch
%                 continue;
%             end
%         end
%         if ~isempty(az)
%             plotted_sats_count = plotted_sats_count + 1;
%             all_sats_data(plotted_sats_count).ID = sid;
%             all_sats_data(plotted_sats_count).azimuth = az;
%             all_sats_data(plotted_sats_count).elevation = el;
%         end
%     end
%     if plotted_sats_count > 0
%         figure('Name', '所有卫星天空图 (带起止点标记)');
%         pax = polaraxes;
%         set(pax, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise', 'RAxisLocation', 90, 'RDir', 'reverse');
%         rlim(pax, [0 90]);
%         rticks(pax, [0 30 60 90]);
%         rticklabels(pax, {'90', '60', '30', '0'});
%         pax.RAxis.Label.String = '仰角 (°)';
%         title('所有可见GPS/北斗卫星天空轨迹图');
%         hold(pax, 'on');
%         plot_handles = gobjects(plotted_sats_count, 1);
%         legend_entries = cell(1, plotted_sats_count);
%         colors = lines(plotted_sats_count);
%         for i = 1:plotted_sats_count
%             theta = deg2rad(all_sats_data(i).azimuth);
% 
% 
% 
% %             r = 90 - all_sats_data(i).elevation;
%             r = all_sats_data(i).elevation;
% 
%             plot_handles(i) = polarplot(pax, theta, r, '-', 'LineWidth', 2.0, 'Color', colors(i,:));
%             polarplot(pax, theta(1), r(1), 'o', 'MarkerEdgeColor', colors(i,:), 'MarkerSize', 8, 'LineWidth', 1.5);
%             polarplot(pax, theta(end), r(end), 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:), 'MarkerSize', 8);
%             text(pax, theta(end), r(end), ['  ' all_sats_data(i).ID], 'FontSize', 10, 'FontWeight', 'bold', 'Color', colors(i,:));
%             legend_entries{i} = all_sats_data(i).ID;
%         end
%         legend(plot_handles, legend_entries, 'Location', 'eastoutside');
%         hold(pax, 'off');
%     end
%     fprintf('✅ 所有分析和绘图任务完成！\n\n');
% end








