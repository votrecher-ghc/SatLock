function varargout = gesture_template_library(action, varargin)
% GESTURE_TEMPLATE_LIBRARY
% Shared gesture template utilities for synthetic generation, GT building,
% and shape matching.

action = lower(string(action));

switch action
    case "all"
        varargout{1} = {'A', 'C', 'M', 'Star', 'L', 'X', 'Z', 'N', 'V', 'Rectangle', 'LeftSwipe', 'RightSwipe'};

    case "label"
        varargout{1} = canonical_label(varargin{1});

    case "stages"
        label = canonical_label(varargin{1});
        stages = raw_stages(label);
        if numel(varargin) >= 2 && ~isempty(varargin{2})
            span_cfg = normalize_span_cfg(varargin{2});
            [stages, meta] = normalize_stages_to_span(stages, span_cfg.max_span_x, span_cfg.max_span_y);
        else
            meta = stage_meta(stages, 1.0);
        end
        varargout{1} = stages;
        if nargout >= 2
            varargout{2} = meta;
        end

    case "groundtruth"
        label = canonical_label(varargin{1});
        num_samples = varargin{2};
        span_cfg = [];
        if numel(varargin) >= 3
            span_cfg = varargin{3};
        end
        [gt_x, gt_y, gt_pen_down, meta] = build_ground_truth(label, num_samples, span_cfg);
        varargout{1} = gt_x;
        varargout{2} = gt_y;
        varargout{3} = gt_pen_down;
        if nargout >= 4
            varargout{4} = meta;
        end

    case "trace"
        label = canonical_label(varargin{1});
        n_points = 120;
        if numel(varargin) >= 2 && ~isempty(varargin{2})
            n_points = varargin{2};
        end
        span_cfg = [];
        if numel(varargin) >= 3
            span_cfg = varargin{3};
        end
        [trace_x, trace_y] = build_pen_trace(label, n_points, span_cfg);
        varargout{1} = trace_x;
        varargout{2} = trace_y;

    case "span"
        label = canonical_label(varargin{1});
        span_cfg = [];
        if numel(varargin) >= 2
            span_cfg = varargin{2};
        end
        [dx, dy] = trace_span(label, span_cfg);
        varargout{1} = dx;
        varargout{2} = dy;

    otherwise
        error('gesture_template_library:UnsupportedAction', 'Unsupported action: %s', action);
end

end

function label = canonical_label(label_in)
label = upper(strtrim(char(string(label_in))));
switch label
    case {'STAR', 'PENTAGRAM'}
        label = 'Star';
    case {'RECTANGLE', 'RECT', 'BOX'}
        label = 'Rectangle';
    case {'LEFTSWIPE', 'LEFT_SWIPE', 'LSWIPE', 'LEFT', '左滑'}
        label = 'LeftSwipe';
    case {'RIGHTSWIPE', 'RIGHT_SWIPE', 'RSWIPE', 'RIGHT', '右滑'}
        label = 'RightSwipe';
    otherwise
        if strcmp(label, 'V')
            label = 'V';
        elseif any(strcmp(label, {'A', 'B', 'C', 'L', 'M', 'N', 'X', 'Z'}))
            label = label;
        else
            error('gesture_template_library:UnsupportedLabel', ...
                'Unsupported gesture template: %s', char(label_in));
        end
end
end

function span_cfg = normalize_span_cfg(span_cfg)
if isempty(span_cfg) || ~isstruct(span_cfg)
    span_cfg = struct();
end
if ~isfield(span_cfg, 'max_span_x') || isempty(span_cfg.max_span_x)
    span_cfg.max_span_x = 0.50;
end
if ~isfield(span_cfg, 'max_span_y') || isempty(span_cfg.max_span_y)
    span_cfg.max_span_y = 0.50;
end
end

function stages = raw_stages(label)
switch label
    case 'A'
        P1 = [-0.40, -0.50];
        P2 = [0.00, 0.50];
        P3 = [0.40, -0.50];
        P4 = [-0.20, -0.10];
        P5 = [0.20, -0.10];
        stages = {
            P1, P2, 1.5, true;
            P2, P2, 3.0, true;
            P2, P3, 1.5, true;
            P3, P4, 3.0, false;
            P4, P5, 1.5, true
        };

    case 'B'
        P1 = [-0.40, -1.00];
        P2 = [-0.40, 1.00];
        P3 = [1.50, 0.40];
        P4 = [-0.40, 0.00];
        P5 = [1.50, 0.00];
        P6 = [-0.40, -1.00];
        stages = {
            P1, P2, 1.5, true;
            P2, P2, 3.0, true;
            P2, P3, 1.5, true;
            P3, P3, 3.0, true;
            P3, P4, 3.0, true;
            P4, P4, 3.0, true;
            P4, P5, 1.5, true;
            P5, P5, 3.0, true;
            P5, P6, 1.5, true
        };

    case 'C'
        P1 = [0.35, 0.40];
        P2 = [-0.35, 0.40];
        P3 = [-0.35, -0.40];
        P4 = [0.35, -0.40];
        stages = {
            P1, P2, 1.5, true;
            P2, P2, 3.0, true;
            P2, P3, 1.8, true;
            P3, P3, 3.0, true;
            P3, P4, 1.5, true
        };

    case 'L'
        P1 = [-0.30, 0.40];
        P2 = [-0.30, -0.40];
        P3 = [0.30, -0.40];
        stages = {
            P1, P2, 1.5, true;
            P2, P2, 3.0, true;
            P2, P3, 1.5, true
        };

    case 'M'
        P1 = [-0.40, -0.40];
        P2 = [-0.40, 0.40];
        P3 = [0.00, 0.00];
        P4 = [0.40, 0.40];
        P5 = [0.40, -0.40];
        stages = {
            P1, P2, 1.5, true;
            P2, P2, 3.0, true;
            P2, P3, 1.5, true;
            P3, P3, 3.0, true;
            P3, P4, 1.5, true;
            P4, P4, 3.0, true;
            P4, P5, 1.5, true
        };

    case 'N'
        P1 = [-0.30, -0.40];
        P2 = [-0.30, 0.40];
        P3 = [0.30, -0.40];
        P4 = [0.30, 0.40];
        stages = {
            P1, P2, 1.5, true;
            P2, P2, 3.0, true;
            P2, P3, 1.5, true;
            P3, P3, 3.0, true;
            P3, P4, 1.5, true
        };

    case 'Rectangle'
        P1 = [-0.35, 0.30];
        P2 = [0.35, 0.30];
        P3 = [0.35, -0.30];
        P4 = [-0.35, -0.30];
        stages = {
            P1, P2, 1.4, true;
            P2, P2, 2.0, true;
            P2, P3, 1.2, true;
            P3, P3, 2.0, true;
            P3, P4, 1.4, true;
            P4, P4, 2.0, true;
            P4, P1, 1.2, true
        };

    case 'LeftSwipe'
        P1 = [0.38, 0.00];
        P2 = [-0.38, 0.00];
        stages = {
            P1, P2, 1.6, true
        };

    case 'RightSwipe'
        P1 = [-0.38, 0.00];
        P2 = [0.38, 0.00];
        stages = {
            P1, P2, 1.6, true
        };

    case 'Star'
        P1 = [-0.30, -0.45];
        P2 = [0.00, 0.55];
        P3 = [0.30, -0.45];
        P4 = [-0.48, 0.15];
        P5 = [0.48, 0.15];
        stages = {
            P1, P2, 1.5, true;
            P2, P2, 3.0, true;
            P2, P3, 1.5, true;
            P3, P3, 3.0, true;
            P3, P4, 1.5, true;
            P4, P4, 3.0, true;
            P4, P5, 1.5, true;
            P5, P5, 3.0, true;
            P5, P1, 1.5, true
        };

    case 'V'
        P1 = [-0.35, 0.40];
        P2 = [0.00, -0.40];
        P3 = [0.35, 0.40];
        stages = {
            P1, P2, 1.5, true;
            P2, P2, 3.0, true;
            P2, P3, 1.5, true
        };

    case 'X'
        P1 = [-0.30, 0.40];
        P2 = [0.30, -0.40];
        P3 = [0.30, 0.40];
        P4 = [-0.30, -0.40];
        stages = {
            P1, P2, 1.5, true;
            P2, P3, 3.0, false;
            P3, P4, 1.5, true
        };

    case 'Z'
        P1 = [-0.30, 0.40];
        P2 = [0.30, 0.40];
        P3 = [-0.30, -0.40];
        P4 = [0.70, -0.40];
        stages = {
            P1, P2, 1.5, true;
            P2, P2, 3.0, true;
            P2, P3, 1.5, true;
            P3, P3, 3.0, true;
            P3, P4, 1.5, true
        };

    otherwise
        error('gesture_template_library:UnsupportedLabel', ...
            'Unsupported gesture template: %s', label);
end
end

function [stages_out, meta] = normalize_stages_to_span(stages_in, max_span_x, max_span_y)
stages_out = stages_in;
meta = struct('scale', 1.0, 'span_x', NaN, 'span_y', NaN);
if isempty(stages_in)
    return;
end

pts = zeros(0, 2);
for i = 1:size(stages_in, 1)
    pts(end+1, :) = stages_in{i, 1}; %#ok<AGROW>
    pts(end+1, :) = stages_in{i, 2}; %#ok<AGROW>
end

min_xy = min(pts, [], 1);
max_xy = max(pts, [], 1);
span = max_xy - min_xy;
center = 0.5 * (min_xy + max_xy);

sx = max_span_x / max(span(1), eps);
sy = max_span_y / max(span(2), eps);
s = min([sx, sy, 1.0]);

for i = 1:size(stages_in, 1)
    p1 = stages_in{i, 1};
    p2 = stages_in{i, 2};
    stages_out{i, 1} = (p1 - center) * s;
    stages_out{i, 2} = (p2 - center) * s;
end

meta = stage_meta(stages_out, s);
end

function meta = stage_meta(stages, scale)
pts = zeros(0, 2);
for i = 1:size(stages, 1)
    pts(end+1, :) = stages{i, 1}; %#ok<AGROW>
    pts(end+1, :) = stages{i, 2}; %#ok<AGROW>
end

min_xy = min(pts, [], 1);
max_xy = max(pts, [], 1);
span = max_xy - min_xy;

meta = struct();
meta.scale = scale;
meta.span_x = span(1);
meta.span_y = span(2);
end

function [gt_x, gt_y, gt_pen_down, meta] = build_ground_truth(label, num_samples, span_cfg)
span_cfg = normalize_span_cfg(span_cfg);
sampling_rate = 25;
start_ratio = 0.30;

[stages, stage_info] = gesture_template_library('stages', label, span_cfg);

gt_x = nan(num_samples, 1);
gt_y = nan(num_samples, 1);
gt_pen_down = false(num_samples, 1);

total_dur = 0;
for k = 1:size(stages, 1)
    total_dur = total_dur + stages{k, 3};
end

start_idx = max(1, round(num_samples * start_ratio));
end_idx = min(num_samples, start_idx + round(total_dur * sampling_rate));

for t_idx = start_idx:end_idx
    dt = (t_idx - start_idx) / sampling_rate;
    elapsed = 0;
    for k = 1:size(stages, 1)
        p1 = stages{k, 1};
        p2 = stages{k, 2};
        dur = stages{k, 3};
        pen = stages{k, 4};
        if dt <= (elapsed + dur)
            alpha = max(0, min(1, (dt - elapsed) / max(dur, eps)));
            pos = p1 + (p2 - p1) * alpha;
            gt_x(t_idx) = pos(1);
            gt_y(t_idx) = pos(2);
            gt_pen_down(t_idx) = logical(pen);
            break;
        end
        elapsed = elapsed + dur;
    end
end

meta = struct();
meta.sampling_rate = sampling_rate;
meta.start_ratio = start_ratio;
meta.start_idx = start_idx;
meta.end_idx = end_idx;
meta.stage_info = stage_info;
meta.label = label;
end

function [trace_x, trace_y] = build_pen_trace(label, n_points, span_cfg)
span_cfg = normalize_span_cfg(span_cfg);
[stages, ~] = gesture_template_library('stages', label, span_cfg);

pts = zeros(0, 2);
for i = 1:size(stages, 1)
    if ~logical(stages{i, 4})
        continue;
    end
    p1 = stages{i, 1};
    p2 = stages{i, 2};
    seg_len = norm(p2 - p1);
    if seg_len < 1e-9
        if isempty(pts)
            pts = p1;
        else
            pts(end+1, :) = p1; %#ok<AGROW>
        end
        continue;
    end
    n_seg = max(4, round(24 * seg_len / max(0.25, seg_len)));
    s = linspace(0, 1, n_seg).';
    seg = p1 + (p2 - p1) .* s;
    if ~isempty(pts)
        seg = seg(2:end, :);
    end
    pts = [pts; seg]; %#ok<AGROW>
end

if isempty(pts)
    trace_x = nan(n_points, 1);
    trace_y = nan(n_points, 1);
    return;
end

[trace_x, trace_y] = resample_polyline(pts(:, 1), pts(:, 2), n_points);
end

function [dx, dy] = trace_span(label, span_cfg)
[trace_x, trace_y] = gesture_template_library('trace', label, 200, span_cfg);
valid = isfinite(trace_x) & isfinite(trace_y);
if ~any(valid)
    dx = NaN;
    dy = NaN;
    return;
end
dx = max(trace_x(valid)) - min(trace_x(valid));
dy = max(trace_y(valid)) - min(trace_y(valid));
end

function [xr, yr] = resample_polyline(x, y, n_points)
x = x(:);
y = y(:);
valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);

if isempty(x)
    xr = nan(n_points, 1);
    yr = nan(n_points, 1);
    return;
end
if numel(x) == 1
    xr = repmat(x, n_points, 1);
    yr = repmat(y, n_points, 1);
    return;
end

d = hypot(diff(x), diff(y));
keep = [true; d > eps];
x = x(keep);
y = y(keep);
if numel(x) == 1
    xr = repmat(x, n_points, 1);
    yr = repmat(y, n_points, 1);
    return;
end
d = hypot(diff(x), diff(y));
s = [0; cumsum(d)];
if s(end) <= eps
    xr = repmat(x(1), n_points, 1);
    yr = repmat(y(1), n_points, 1);
    return;
end

sq = linspace(0, s(end), n_points).';
xr = interp1(s, x, sq, 'linear');
yr = interp1(s, y, sq, 'linear');
end
