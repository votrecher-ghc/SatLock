% =========================================================================
% generate_ideal_multi_shape.m 
% 功能: 通用仿真信号生成器 (修复顶点信号连续性 & 完整绘图版)
% 描述: 
%   该函数用于在真实的观测数据基础上，注入模拟的手势遮挡信号。
%   通过修改 STAGES 逻辑，确保手部在转折点停顿时信号不会瞬间回升到基线。
%   保留完整的 Ground Truth 绘图功能。
%
% [调用格式]:
%   obs_sim = generate_ideal_multi_shape(obs_data, nav_data, TARGET_LETTER);
% =========================================================================
function obs_data = generate_ideal_multi_shape(obs_data, nav_data, TARGET_LETTER)

% --- [Part 0] 环境与数据检查 ---
addpath(genpath('calculate_clock_bias_and_positon'));
addpath(genpath('nav_parse'));
fprintf('--> 启动仿真 V7.1 (物理连续版): 目标 [%s]...\n', TARGET_LETTER);

%% ================= [Part 1] 仿真参数设置 =================
% 1. 信号参数
SIM.baseline_db    = 45;       % [基线] 平稳时的信号强度 (dB)
SIM.drop_depth_db  = 15;       % [波动] 遮挡时的下降深度 (dB)
SIM.noise_sigma    = 0.02;     % [噪声] 微小抖动

% 2. 几何参数
SIM.gesture_height = 0.30;     % [物理] 手势平面高度 (米)
SIM.arm_width      = 0.15;     % [物理] 手臂有效遮挡宽度 (米)
SIM.body_pos       = [0.0, -1.0]; % [物理] 身体位置 (接收机正南方1米)

% 3. 形状定义 (重点：同一字母内部的停留阶段设为 true，确保信号连续)
switch TARGET_LETTER
    case 'A'
        P1 = [-0.40, -0.50]; P2 = [ 0.00,  0.50]; P3 = [ 0.40, -0.50];
        P4 = [-0.20, -0.10]; P5 = [ 0.20, -0.10];
        STAGES = {
            P1, P2, 1.5, true;   % 左侧斜线上升
            P2, P2, 3.0, true;   % 顶点停留 (修正：保持信号下降)
            P2, P3, 1.5, true;   % 右侧斜线下降
            P3, P4, 3.0, false;  % 换笔移动 (抬手)
            P4, P5, 1.5, true    % 中间横线
        };
    case 'B'
        P1 = [-0.40, -1]; P2 = [-0.40,  1]; P3 = [ 1.5,  0.40]; 
        P4 = [-0.40,  0.00]; P5 = [ 1.5,  0.00]; P6 = [-0.40, -1]; 
        STAGES = {
            P1, P2, 1.5, true; P2, P2, 3.0, true;
            P2, P3, 1.5, true; P3, P3, 3.0, true;
            P3, P4, 3.0, true; P4, P4, 3.0, true;
            P4, P5, 1.5, true; P5, P5, 3.0, true;
            P5, P6, 1.5, true
        };
    case 'M'
        P1 = [-0.40, -0.40]; P2 = [-0.40,  0.40]; P3 = [ 0.00,  0.00];
        P4 = [ 0.40,  0.40]; P5 = [ 0.40, -0.40];
        STAGES = {
            P1, P2, 1.5, true;  P2, P2, 3.0, true;
            P2, P3, 1.5, true;  P3, P3, 3.0, true;
            P3, P4, 1.5, true;  P4, P4, 3.0, true;
            P4, P5, 1.5, true
        };
    case 'Star'
        P1 = [-0.30, -0.45]; P2 = [ 0.00,  0.55]; P3 = [ 0.30, -0.45];
        P4 = [-0.48,  0.15]; P5 = [ 0.48,  0.15];
        STAGES = {
            P1, P2, 1.5, true; P2, P2, 3.0, true;
            P2, P3, 1.5, true; P3, P3, 3.0, true;
            P3, P4, 1.5, true; P4, P4, 3.0, true;
            P4, P5, 1.5, true; P5, P5, 3.0, true;
            P5, P1, 1.5, true
        };
    case 'L'
        P_TL=[-0.3,0.4]; P_BL=[-0.3,-0.4]; P_BR=[0.3,-0.4];
        STAGES = {
            P_TL, P_BL, 1.5, true; 
            P_BL, P_BL, 3.0, true; 
            P_BL, P_BR, 1.5, true  
        };
    case 'X'
        P_TL=[-0.3,0.4]; P_TR=[0.3,0.4]; P_BL=[-0.3,-0.4]; P_BR=[0.3,-0.4];
        STAGES = {
            P_TL, P_BR, 1.5, true; 
            P_BR, P_TR, 3.0, false; 
            P_TR, P_BL, 1.5, true  
        };
    case 'Z'
        P_TL=[-0.3,0.4]; P_TR=[0.3,0.4]; P_BL=[-0.3,-0.4]; P_BR=[0.7,-0.4];
        STAGES = {
            P_TL, P_TR, 1.5, true; 
            P_TR, P_TR, 3.0, true; 
            P_TR, P_BL, 1.5, true; 
            P_BL, P_BL, 3.0, true; 
            P_BL, P_BR, 1.5, true  
        };
    case 'N'
        P_BL=[-0.3,-0.4]; P_TL=[-0.3,0.4]; P_BR=[0.3,-0.4]; P_TR=[0.3,0.4];
        STAGES = {
            P_BL, P_TL, 1.5, true; 
            P_TL, P_TL, 3.0, true; 
            P_TL, P_BR, 1.5, true; 
            P_BR, P_BR, 3.0, true; 
            P_BR, P_TR, 1.5, true  
        };
    otherwise
        error('未定义的字母: %s', TARGET_LETTER);
end

% 计算总仿真时长
total_sim_duration = 0;
for k=1:size(STAGES, 1), total_sim_duration = total_sim_duration + STAGES{k,3}; end

%% ================= [Part 2] 数据准备 =================
ideal_obs = obs_data;
raw_times = [ideal_obs.time];
num_samples = length(raw_times);
start_idx = round(num_samples * 0.3); 
sampling_rate = 25; 
end_idx = start_idx + round(total_sim_duration * sampling_rate);

% 提取支持的卫星系统 (G/C/R/E/J)
all_sat_ids = {}; for i=1:min(100, length(ideal_obs)), if ~isempty(ideal_obs(i).data), all_sat_ids = [all_sat_ids, fieldnames(ideal_obs(i).data)']; end; end
unique_sat_ids = unique(all_sat_ids); 
valid_sats = {}; for i=1:length(unique_sat_ids), sid=unique_sat_ids{i}; if ismember(sid(1),['G','C','R','E','J']), valid_sats{end+1}=sid; end; end

% 计算接收机参考坐标
rec_pos_acc = [0,0,0]; count = 0;
for t = start_idx : 10 : end_idx
    try [rp,~,~]=calculate_receiver_position(obs_data, nav_data, t); rec_pos_acc=rec_pos_acc+rp; count=count+1; catch, end
end
rec_pos_mean = rec_pos_acc / count;
[lat0, lon0, alt0] = ecef2geodetic(rec_pos_mean(1), rec_pos_mean(2), rec_pos_mean(3));

%% ================= [Part 3] 信号注入核心循环 =================
num_sim_pts = end_idx - start_idx + 1;
gt_trace_x = NaN(num_sim_pts, 1);
gt_trace_y = NaN(num_sim_pts, 1);
gt_pen_down = false(num_sim_pts, 1);

for t_idx = 1 : num_samples
    is_pen_down = false;
    current_hand_pos = [NaN, NaN];
    
    if t_idx >= start_idx && t_idx <= end_idx
        dt = (t_idx - start_idx) / sampling_rate;
        elapsed = 0;
        for k = 1:size(STAGES, 1)
            dur = STAGES{k,3};
            if dt <= (elapsed + dur)
                local_prog = max(0, min(1, (dt - elapsed) / dur));
                current_hand_pos = STAGES{k,1} + (STAGES{k,2} - STAGES{k,1}) * local_prog;
                is_pen_down = STAGES{k,4};
                
                trace_idx = t_idx - start_idx + 1;
                gt_trace_x(trace_idx) = current_hand_pos(1);
                gt_trace_y(trace_idx) = current_hand_pos(2);
                gt_pen_down(trace_idx) = is_pen_down;
                break;
            end
            elapsed = elapsed + dur;
        end
    end
    
    if isempty(ideal_obs(t_idx).data), continue; end
    try [~, ~, sat_states] = calculate_receiver_position(obs_data, nav_data, t_idx); catch, continue; end
    
    for s = 1:length(valid_sats)
        sid = valid_sats{s};
        if ~isfield(ideal_obs(t_idx).data, sid), continue; end
        sim_val = SIM.baseline_db + randn() * SIM.noise_sigma;
        
        if is_pen_down && ~isnan(current_hand_pos(1)) && isfield(sat_states, sid)
            sat_p = sat_states.(sid).position;
            [e, n, u] = ecef2enu(sat_p(1)-rec_pos_mean(1), sat_p(2)-rec_pos_mean(2), sat_p(3)-rec_pos_mean(3), lat0, lon0, alt0);
            if u > 0
                scale = SIM.gesture_height / u;
                P = [scale * e, scale * n];
                A = current_hand_pos; B = SIM.body_pos;
                vec_AB = B - A; vec_AP = P - A;
                len_sq = sum(vec_AB.^2);
                if len_sq > 0
                    t_proj = max(0, min(1, dot(vec_AP, vec_AB) / len_sq));
                    dist_to_arm = norm(P - (A + t_proj * vec_AB));
                    if dist_to_arm < SIM.arm_width
                        sim_val = SIM.baseline_db - SIM.drop_depth_db + randn() * SIM.noise_sigma;
                    end
                end
            end
        end
        snr_struct = ideal_obs(t_idx).data.(sid).snr;
        fields = fieldnames(snr_struct);
        for f = 1:length(fields), ideal_obs(t_idx).data.(sid).snr.(fields{f}) = sim_val; end
    end
end
obs_data = ideal_obs;

%% ================= [Part 4] Ground Truth 可视化 =================
fprintf('--> 生成 Ground Truth 轨迹图...\n');
figure('Name', sprintf('Ideal Ground Truth Trajectory [%s]', TARGET_LETTER), 'Position', [100, 100, 600, 600], 'Color', 'w');
ax = axes; hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
xlabel('East (m)'); ylabel('North (m)');

% 1. 画接收机位置和身体参考点
plot(ax, 0, 0, '^', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', 'Receiver');
plot(ax, SIM.body_pos(1), SIM.body_pos(2), 'bs', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'Body/Shoulder Ref');

% 2. 画轨迹 (区分落笔和抬手)
plot_x_up = gt_trace_x; plot_y_up = gt_trace_y;
plot_x_up(gt_pen_down) = NaN; plot_y_up(gt_pen_down) = NaN;
plot_x_down = gt_trace_x; plot_y_down = gt_trace_y;
plot_x_down(~gt_pen_down) = NaN; plot_y_down(~gt_pen_down) = NaN;

plot(ax, plot_x_up, plot_y_up, 'k--', 'LineWidth', 1, 'Color', [0.6 0.6 0.6], 'DisplayName', 'Pen Up (Interval Move)');
plot(ax, plot_x_down, plot_y_down, 'b-', 'LineWidth', 3, 'DisplayName', 'Pen Down (Target Stroke)');

% 3. 标记起点和终点
f_idx = find(~isnan(gt_trace_x), 1, 'first');
if ~isempty(f_idx), plot(ax, gt_trace_x(f_idx), gt_trace_y(f_idx), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10, 'DisplayName', 'Start'); end
l_idx = find(~isnan(gt_trace_x), 1, 'last');
if ~isempty(l_idx), plot(ax, gt_trace_x(l_idx), gt_trace_y(l_idx), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'DisplayName', 'End'); end

title({sprintf('Ideal Ground Truth Trajectory (Letter %s)', TARGET_LETTER), 'Blue: Signal Drop Regions, Gray: Intervals'});
legend('Location', 'best');

% 自动缩放
valid_x = gt_trace_x(~isnan(gt_trace_x));
valid_y = gt_trace_y(~isnan(gt_trace_y));
if ~isempty(valid_x), m_range = max([max(abs(valid_x)), max(abs(valid_y)), 0.6]); else, m_range = 0.6; end
xlim([-m_range*1.2, m_range*1.2]); ylim([-m_range*1.2, m_range*1.2] + SIM.body_pos(2)/2);

fprintf('✅ 仿真与绘图全部完成。\n');
end