function [obs_data, nav_data, meta] = prepare_input_data(obs_filepath, nav_filepath, sim_cfg)
% PREPARE_INPUT_DATA
% 第一层：数据处理
% 功能：解析 obs/nav，并可选注入理想仿真手势数据。

if nargin < 3 || isempty(sim_cfg)
    sim_cfg = struct();
end

if ~isfield(sim_cfg, 'enable')
    sim_cfg.enable = true;
end
if ~isfield(sim_cfg, 'target_letter')
    sim_cfg.target_letter = 'N';
end
if ~isfield(sim_cfg, 'max_span_x')
    sim_cfg.max_span_x = 0.50; % 左右总跨度 50cm
end
if ~isfield(sim_cfg, 'max_span_y')
    sim_cfg.max_span_y = 0.50; % 上下总跨度 50cm
end

fprintf('--> 正在解析观测文件: %s\n', obs_filepath);
obs_data = parse_rinex_obs(obs_filepath);
fprintf('--> 正在解析导航文件: %s\n', nav_filepath);
nav_data = parse_rinex_nav_multi_gnss(nav_filepath);

if sim_cfg.enable
    letter = char(sim_cfg.target_letter);
    fprintf('--> 注入理想仿真手势: %s (跨度约束 %.2fm x %.2fm)\n', ...
        upper(letter), sim_cfg.max_span_x, sim_cfg.max_span_y);
    obs_data = generate_ideal_multi_shape(obs_data, nav_data, letter, sim_cfg);
end

meta = struct();
meta.obs_filepath = obs_filepath;
meta.nav_filepath = nav_filepath;
meta.sim_enabled = logical(sim_cfg.enable);
meta.target_letter = sim_cfg.target_letter;
meta.max_span_x = sim_cfg.max_span_x;
meta.max_span_y = sim_cfg.max_span_y;
end

