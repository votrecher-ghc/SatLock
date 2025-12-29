% ============== calculate_and_plot_skyplot.m (适配GPS+北斗修正版) ==============
function [az, el, sat_id_out] = calculate_and_plot_skyplot(obs_data, nav_data, target_satellite_id)
% 计算并绘制单颗卫星的天空图，适配GPS和北斗导航数据格式。

fprintf('--> 开始为卫星 %s 计算天空图数据...\n', target_satellite_id);

% --- 获取接收机位置 ---
try
    [receiver_pos_ecef, ~, ~] = calculate_receiver_position(obs_data, nav_data, 1);
catch ME
    error('计算接收机位置失败：%s', ME.message);
end

% --- 遍历历元 ---
num_epochs = length(obs_data);
az = NaN(1, num_epochs);
el = NaN(1, num_epochs);
valid_count = 0;

for i = 1:num_epochs
    current_epoch = obs_data(i);
    if ~isfield(current_epoch.data, target_satellite_id), continue; end

    sat_obs = current_epoch.data.(target_satellite_id);
    if ~(isfield(sat_obs, 'pseudorange') && isfield(sat_obs.pseudorange, 'C1C')), continue; end
    if isnan(sat_obs.pseudorange.C1C), continue; end

    try
        % 调用已更新的、支持双系统的 calculate_satellite_state 函数
        [sat_pos_ecef, ~, ~] = calculate_satellite_state(current_epoch.time, ...
            sat_obs.pseudorange.C1C, target_satellite_id, nav_data);

        [lat, lon, ~] = ecef2geodetic(receiver_pos_ecef(1), receiver_pos_ecef(2), receiver_pos_ecef(3));
        vec_ecef = sat_pos_ecef - receiver_pos_ecef;
        [e, n, u] = ecef2enu(vec_ecef(1), vec_ecef(2), vec_ecef(3), lat, lon, 0);

        current_az = atan2d(e, n); if current_az < 0, current_az = current_az + 360; end
        current_el = asind(u / norm([e, n, u]));

        if current_el >= 0
            valid_count = valid_count + 1;
            az(valid_count) = current_az;
            el(valid_count) = current_el;
        end
    catch ME
        % 如果calculate_satellite_state找不到星历而报错，这里会捕捉并静默跳过
        % 这样可以避免因单点数据问题导致整个程序中断
        % fprintf('警告：在处理卫星 %s 的历元 %d 时跳过，原因: %s\n', target_satellite_id, i, ME.message);
        continue;
    end
end

% --- 绘图 ---
if valid_count == 0
    warning('卫星 %s 没有足够的有效观测/导航数据，无法绘图。', target_satellite_id);
    az = []; el = []; sat_id_out = target_satellite_id;
    return;
end

az = az(1:valid_count);
el = el(1:valid_count);
sat_id_out = target_satellite_id;

figure('Name', sprintf('卫星 %s 天空图', target_satellite_id));
skyplot(az, el, 'LineWidth', 2);
title(sprintf('卫星 %s 天空轨迹图', target_satellite_id));
fprintf('✅ 卫星 %s 的天空图绘制完成！\n\n', target_satellite_id);
end