function [pair_tbl, case_tbl, out_dir] = benchmark_sim_oracle_multi_dataset(data_root, letters, max_pairs, user_cfg)
% BENCHMARK_SIM_ORACLE_MULTI_DATASET
% Discover matched obs/nav pairs, then run benchmark_sim_oracle_vs_baselines.
%
% Example:
%   [pair_tbl, case_tbl, out_dir] = benchmark_sim_oracle_multi_dataset( ...
%       'D:/Matproject/SatLock/data', {'A','L','Z','N','Star'}, 4);

if nargin < 1 || isempty(data_root)
    data_root = fullfile('data');
end
if nargin < 2 || isempty(letters)
    letters = {'A', 'B', 'M', 'Star', 'L', 'X', 'Z', 'N'};
end
if nargin < 3 || isempty(max_pairs)
    max_pairs = 4;
end
if nargin < 4
    user_cfg = struct();
end
if ischar(letters) || isstring(letters)
    letters = cellstr(letters);
end
max_pairs = max(1, round(max_pairs));

this_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(this_dir);
addpath(genpath(repo_dir));

if ~isfolder(data_root)
    data_root = fullfile(repo_dir, data_root);
end
if ~isfolder(data_root)
    error('data_root not found: %s', data_root);
end

results_root = fullfile(this_dir, 'results');
if ~exist(results_root, 'dir')
    mkdir(results_root);
end
stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
out_dir = fullfile(results_root, ['sim_oracle_multi_', stamp]);
mkdir(out_dir);

fprintf('\n== Multi-dataset Sim-Oracle Benchmark ==\n');
fprintf('data_root : %s\n', data_root);
fprintf('letters   : %s\n', strjoin(upper(letters), ', '));
fprintf('max_pairs : %d\n\n', max_pairs);

pairs = discover_obs_nav_pairs(data_root, max_pairs);
if isempty(pairs)
    error('No reliable obs/nav pairs found under: %s', data_root);
end

fprintf('[Pair Discovery] %d pairs selected:\n', numel(pairs));
for i = 1:numel(pairs)
    fprintf('  [%d] %s\n', i, pairs(i).folder);
    fprintf('      obs: %s\n', pairs(i).obs_path);
    fprintf('      nav: %s (%s)\n', pairs(i).nav_path, pairs(i).pair_rule);
end
fprintf('\n');

case_tbl = table();
fail_logs = {};

for i = 1:numel(pairs)
    p = pairs(i);
    fprintf('[Pair %d/%d] Running benchmark...\n', i, numel(pairs));
    fprintf('  folder: %s\n', p.folder);
    try
        t = benchmark_sim_oracle_vs_baselines(p.obs_path, p.nav_path, letters, user_cfg);
        if isempty(t)
            fail_logs{end+1, 1} = sprintf('Pair %d returned empty table.', i); %#ok<AGROW>
            continue;
        end

        n = height(t);
        t.pair_id = repmat(i, n, 1);
        t.folder = repmat({p.folder}, n, 1);
        t.obs_path = repmat({p.obs_path}, n, 1);
        t.nav_path = repmat({p.nav_path}, n, 1);
        t.pair_rule = repmat({p.pair_rule}, n, 1);
        t.known_rmse_gain_m = t.baseline_rmse_m - t.known_rmse_m;
        t.auto_rmse_gain_m = t.baseline_rmse_m - t.auto_rmse_m;
        run_dir = fileparts(t.plot_path{1});
        t.run_dir = repmat({run_dir}, n, 1);

        case_tbl = [case_tbl; t]; %#ok<AGROW>
    catch ME
        warn_msg = sprintf('Pair %d failed: %s', i, ME.message);
        warning('%s', warn_msg);
        fail_logs{end+1, 1} = warn_msg; %#ok<AGROW>
    end
end

if isempty(case_tbl)
    error('All pair runs failed. See MATLAB warnings for details.');
end

pair_tbl = aggregate_pair_metrics(case_tbl);

case_csv = fullfile(out_dir, 'case_summary.csv');
pair_csv = fullfile(out_dir, 'pair_summary.csv');
writetable(case_tbl, case_csv);
writetable(pair_tbl, pair_csv);

if ~isempty(fail_logs)
    fail_txt = fullfile(out_dir, 'fail_logs.txt');
    fid = fopen(fail_txt, 'w');
    if fid > 0
        for i = 1:numel(fail_logs)
            fprintf(fid, '%s\n', fail_logs{i});
        end
        fclose(fid);
    end
end

save(fullfile(out_dir, 'summary.mat'), 'case_tbl', 'pair_tbl', 'pairs', 'fail_logs', 'out_dir');

fprintf('\n== Multi-dataset benchmark complete ==\n');
fprintf('Output dir    : %s\n', out_dir);
fprintf('Case summary  : %s\n', case_csv);
fprintf('Pair summary  : %s\n\n', pair_csv);
disp(pair_tbl);
end

% -------------------------------------------------------------------------
function pair_tbl = aggregate_pair_metrics(case_tbl)
pair_ids = unique(case_tbl.pair_id(:));
n_pair = numel(pair_ids);

pair_tbl = table('Size', [n_pair, 14], ...
    'VariableTypes', {'double','string','string','string','double','double','double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'pair_id','folder','obs_path','nav_path', ...
    'n_cases','mean_baseline_rmse_m','mean_known_rmse_m','mean_auto_rmse_m', ...
    'mean_known_err_pct','mean_auto_err_pct','known_pass20_rate','auto_pass20_rate', ...
    'mean_known_gain_m','mean_auto_gain_m'});

for k = 1:n_pair
    pid = pair_ids(k);
    idx = (case_tbl.pair_id == pid);
    sub = case_tbl(idx, :);

    pair_tbl.pair_id(k) = pid;
    pair_tbl.folder(k) = string(sub.folder{1});
    pair_tbl.obs_path(k) = string(sub.obs_path{1});
    pair_tbl.nav_path(k) = string(sub.nav_path{1});
    pair_tbl.n_cases(k) = height(sub);

    pair_tbl.mean_baseline_rmse_m(k) = mean(sub.baseline_rmse_m, 'omitnan');
    pair_tbl.mean_known_rmse_m(k) = mean(sub.known_rmse_m, 'omitnan');
    pair_tbl.mean_auto_rmse_m(k) = mean(sub.auto_rmse_m, 'omitnan');
    pair_tbl.mean_known_err_pct(k) = mean(sub.known_err_pct, 'omitnan');
    pair_tbl.mean_auto_err_pct(k) = mean(sub.auto_err_pct, 'omitnan');
    pair_tbl.known_pass20_rate(k) = mean(double(sub.known_pass20), 'omitnan');
    pair_tbl.auto_pass20_rate(k) = mean(double(sub.auto_pass20), 'omitnan');
    pair_tbl.mean_known_gain_m(k) = mean(sub.known_rmse_gain_m, 'omitnan');
    pair_tbl.mean_auto_gain_m(k) = mean(sub.auto_rmse_gain_m, 'omitnan');
end

pair_tbl = sortrows(pair_tbl, 'mean_known_gain_m', 'descend');
end

% -------------------------------------------------------------------------
function pairs = discover_obs_nav_pairs(data_root, max_pairs)
dirs = strsplit(genpath(data_root), pathsep);
pairs = struct('folder', {}, 'obs_path', {}, 'nav_path', {}, 'pair_rule', {});

for i = 1:numel(dirs)
    d = strtrim(dirs{i});
    if isempty(d) || ~isfolder(d)
        continue;
    end

    obs_list = dir(fullfile(d, '*.obs'));
    nav_list = dir(fullfile(d, '*.nav'));
    if isempty(obs_list) || isempty(nav_list)
        continue;
    end

    [obs_path, nav_path, pair_rule] = pick_pair_for_folder(d, obs_list, nav_list);
    if isempty(obs_path) || isempty(nav_path)
        continue;
    end

    p = struct();
    p.folder = d;
    p.obs_path = obs_path;
    p.nav_path = nav_path;
    p.pair_rule = pair_rule;
    pairs(end+1) = p; %#ok<AGROW>

    if numel(pairs) >= max_pairs
        break;
    end
end
end

% -------------------------------------------------------------------------
function [obs_path, nav_path, pair_rule] = pick_pair_for_folder(folder_path, obs_list, nav_list)
obs_path = '';
nav_path = '';
pair_rule = '';

obs_names = {obs_list.name};
nav_names = {nav_list.name};

for i = 1:numel(obs_names)
    [~, stem, ~] = fileparts(obs_names{i});
    nav_name = [stem, '.nav'];
    j = find(strcmpi(nav_names, nav_name), 1, 'first');
    if ~isempty(j)
        obs_path = fullfile(folder_path, obs_names{i});
        nav_path = fullfile(folder_path, nav_names{j});
        pair_rule = 'same_stem';
        return;
    end
end

if numel(nav_list) == 1
    pick_idx = choose_preferred_obs(obs_names);
    obs_path = fullfile(folder_path, obs_names{pick_idx});
    nav_path = fullfile(folder_path, nav_names{1});
    pair_rule = 'single_nav_in_folder';
end
end

% -------------------------------------------------------------------------
function pick_idx = choose_preferred_obs(obs_names)
pick_idx = 1;
if isempty(obs_names)
    return;
end

names_low = lower(obs_names);
score = zeros(numel(names_low), 1);

for i = 1:numel(names_low)
    nm = names_low{i};
    if contains(nm, 'fingure') || contains(nm, 'mixfingure')
        score(i) = score(i) + 3;
    end
    if contains(nm, 'a_') || contains(nm, '_a_') || strcmp(nm, 'a.obs')
        score(i) = score(i) + 2;
    end
    if contains(nm, 'l_') || contains(nm, 'z_') || contains(nm, 'n_') || contains(nm, 'star')
        score(i) = score(i) + 1;
    end
end

[~, pick_idx] = max(score);
end
