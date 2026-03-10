function [traj_x, traj_y, traj_t, traj_conf, debug_info] = ...
    run_gesture_analysis_data_driven(obs_waveform, nav_data, step1_res_shaped, user_cfg)
% RUN_GESTURE_ANALYSIS_DATA_DRIVEN
% 第二层：轨迹推演（纯数据驱动）
% 说明：
%   1) 不使用目标字母；
%   2) 不使用 groundtruth；
%   3) 基于 inverse-beam 做物理一致反演。

cfg = default_cfg();
if nargin >= 4 && isstruct(user_cfg)
    cfg = merge_cfg(cfg, user_cfg);
end

% 强制关闭任何模板吸附，避免“看答案”式推演。
cfg.track.template_snap_enable = false;

[traj_x, traj_y, traj_t, traj_conf, debug_info] = ...
    run_gesture_analysis_inverse_beam(obs_waveform, nav_data, step1_res_shaped, cfg);

% 记录可审计标记，便于确认推演过程未使用真值信息。
if ~isstruct(debug_info)
    debug_info = struct();
end
debug_info.inference_mode = 'data_driven_inverse_beam';
debug_info.uses_groundtruth = false;
debug_info.uses_target_letter = false;
end

% -------------------------------------------------------------------------
function cfg = default_cfg()
cfg = struct();

% 调试与可视化
cfg.debug = struct();
cfg.debug.verbose = true;
cfg.debug.plot = true;    % true 时直接弹出图

% 物理与搜索范围：按 50cm x 50cm 感知范围收敛
cfg.model = struct();
cfg.model.max_hand_radius = 0.40;
cfg.model.center_prior_pos = [0.0, 0.0];
cfg.model.center_prior_weight = 0.02;

cfg.grid = struct();
cfg.grid.x_min = -0.35;
cfg.grid.x_max = 0.35;
cfg.grid.y_min = -0.35;
cfg.grid.y_max = 0.35;
cfg.grid.step = 0.015;

% 轨迹约束参数
cfg.track = struct();
cfg.track.lambda_smooth = 12.0;
cfg.track.final_smooth_pts = 3;
cfg.track.max_jump_m = 0.18;
cfg.track.max_speed_mps = 2.0;
cfg.track.use_active_interpolation = true;
cfg.track.use_process_window_output = true;
cfg.track.output_pad_frames = 10;
cfg.track.use_draw_mask_output = true;
cfg.track.use_draw_energy_gate = true;
cfg.track.output_energy_quantile = 0.35;
cfg.track.output_conf_quantile = 0.25;
cfg.track.drawing_conf_quantile = 0.30;
cfg.track.drawing_energy_quantile = 0.45;
cfg.track.enforce_piecewise_linear = false;
cfg.track.template_snap_enable = false;
cfg.track.endpoint_lock_enable = false;
end

% -------------------------------------------------------------------------
function dst = merge_cfg(dst, src)
keys = fieldnames(src);
for i = 1:numel(keys)
    k = keys{i};
    if isstruct(src.(k))
        if ~isfield(dst, k) || ~isstruct(dst.(k))
            dst.(k) = src.(k);
        else
            dst.(k) = merge_cfg(dst.(k), src.(k));
        end
    else
        dst.(k) = src.(k);
    end
end
end

