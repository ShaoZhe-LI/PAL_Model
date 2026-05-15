%% ============================================================
% Manually draw separated x/y axes only
%
% Output:
%   AxisOnly_z1_xy_manual/X_axis_only.pdf
%   AxisOnly_z1_xy_manual/Y_axis_only.pdf
%
% Key idea:
%   Do NOT use MATLAB native axes ticks/labels.
%   Draw axis line, ticks, tick labels, and axis names manually.
%   This avoids cropping and gives exact whitespace control.
%% ============================================================

clear; clc; close all;

%% ===================== output folder =====================
out_dir = 'AxisOnly_z1_xy_manual';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% ===================== common config =====================
cfg = struct();

cfg.xy_min = -0.5;
cfg.xy_max =  0.5;
cfg.ticks  = [-0.5, -0.25, 0, 0.25, 0.5];

cfg.tick_labels = {'-0.5', '-0.25', '0', '0.25', '0.5'};

cfg.font_size_tick  = 34;
cfg.font_size_label = 40;
cfg.line_width      = 3.0;
cfg.tick_len        = 0.045;    % normalized length

% figure size
cfg.x_fig_pos = [100 100 1200 210];
cfg.y_fig_pos = [100 100 300 1200];

%% ===================== x-axis layout =====================
% normalized coordinates in figure canvas
% x-axis upper whitespace controlled by x_axis_y.
% Larger x_axis_y -> less top whitespace.
cfg.x_left   = 0.07;
cfg.x_right  = 0.97;
cfg.x_axis_y = 0.82;

cfg.x_tick_label_y = 0.50;
cfg.x_label_y      = 0.16;

%% ===================== y-axis layout =====================
% y-axis right whitespace controlled by y_axis_x.
% Larger y_axis_x -> less right whitespace.
% ylabel clipping controlled by y_label_x.
cfg.y_bottom = 0.04;
cfg.y_top    = 0.97;
cfg.y_axis_x = 0.82;

cfg.y_tick_label_x = 0.69;
cfg.y_label_x      = 0.16;

%% ===================== draw x-axis =====================
fig_x = draw_x_axis_manual(cfg);
fp_x_pdf = fullfile(out_dir, 'X_axis_only.pdf');
fp_x_png = fullfile(out_dir, 'X_axis_only.png');
export_safe(fig_x, fp_x_pdf, fp_x_png);
close(fig_x);

fprintf('Saved:\n  %s\n  %s\n', fp_x_pdf, fp_x_png);

%% ===================== draw y-axis =====================
fig_y = draw_y_axis_manual(cfg);
fp_y_pdf = fullfile(out_dir, 'Y_axis_only.pdf');
fp_y_png = fullfile(out_dir, 'Y_axis_only.png');
export_safe(fig_y, fp_y_pdf, fp_y_png);
close(fig_y);

fprintf('Saved:\n  %s\n  %s\n', fp_y_pdf, fp_y_png);

fprintf('\nDone.\nOutput folder:\n%s\n', out_dir);

%% ============================================================
% Local functions
%% ============================================================

function fig = draw_x_axis_manual(cfg)

fig = figure( ...
    'Name', 'X_axis_only_manual', ...
    'Color', 'w', ...
    'Position', cfg.x_fig_pos);

ax = axes('Parent', fig, ...
    'Position', [0 0 1 1], ...
    'XLim', [0 1], ...
    'YLim', [0 1], ...
    'Visible', 'off');

hold(ax, 'on');

% axis line
line(ax, [cfg.x_left cfg.x_right], [cfg.x_axis_y cfg.x_axis_y], ...
    'Color', 'k', ...
    'LineWidth', cfg.line_width, ...
    'Clipping', 'off');

% ticks and labels
for ii = 1:numel(cfg.ticks)
    t = cfg.ticks(ii);

    x = map_value(t, cfg.xy_min, cfg.xy_max, cfg.x_left, cfg.x_right);

    % tick line downward
    line(ax, [x x], [cfg.x_axis_y cfg.x_axis_y - cfg.tick_len], ...
        'Color', 'k', ...
        'LineWidth', cfg.line_width, ...
        'Clipping', 'off');

    text(ax, x, cfg.x_tick_label_y, cfg.tick_labels{ii}, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', cfg.font_size_tick, ...
        'Interpreter', 'latex', ...
        'Clipping', 'off');
end

% x label
text(ax, 0.5, cfg.x_label_y, '$x$ (m)', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', cfg.font_size_label, ...
    'Interpreter', 'latex', ...
    'Clipping', 'off');

end

function fig = draw_y_axis_manual(cfg)

fig = figure( ...
    'Name', 'Y_axis_only_manual', ...
    'Color', 'w', ...
    'Position', cfg.y_fig_pos);

ax = axes('Parent', fig, ...
    'Position', [0 0 1 1], ...
    'XLim', [0 1], ...
    'YLim', [0 1], ...
    'Visible', 'off');

hold(ax, 'on');

% axis line
line(ax, [cfg.y_axis_x cfg.y_axis_x], [cfg.y_bottom cfg.y_top], ...
    'Color', 'k', ...
    'LineWidth', cfg.line_width, ...
    'Clipping', 'off');

% ticks and labels
for ii = 1:numel(cfg.ticks)
    t = cfg.ticks(ii);

    y = map_value(t, cfg.xy_min, cfg.xy_max, cfg.y_bottom, cfg.y_top);

    % tick line leftward
    line(ax, [cfg.y_axis_x cfg.y_axis_x - cfg.tick_len], [y y], ...
        'Color', 'k', ...
        'LineWidth', cfg.line_width, ...
        'Clipping', 'off');

    text(ax, cfg.y_tick_label_x, y, cfg.tick_labels{ii}, ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', cfg.font_size_tick, ...
        'Interpreter', 'latex', ...
        'Clipping', 'off');
end

% y label
text(ax, cfg.y_label_x, 0.5, '$y$ (m)', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'Rotation', 90, ...
    'FontSize', cfg.font_size_label, ...
    'Interpreter', 'latex', ...
    'Clipping', 'off');

end

function x = map_value(v, vmin, vmax, omin, omax)

x = omin + (v - vmin) ./ (vmax - vmin) .* (omax - omin);

end

function export_safe(fig_handle, pdf_path, png_path)

set(fig_handle, 'Color', 'w');
set(fig_handle, 'InvertHardcopy', 'off');

pos = get(fig_handle, 'Position');
fig_w_in = pos(3) / 100;
fig_h_in = pos(4) / 100;

set(fig_handle, 'PaperUnits', 'inches');
set(fig_handle, 'PaperSize', [fig_w_in fig_h_in]);
set(fig_handle, 'PaperPosition', [0 0 fig_w_in fig_h_in]);

print(fig_handle, pdf_path, '-dpdf', '-painters');
print(fig_handle, png_path, '-dpng', '-r300');

end