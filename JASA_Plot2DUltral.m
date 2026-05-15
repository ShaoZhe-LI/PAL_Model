%% ============================================================
% MAIN: King only + two xOz 2D plots
%
% 仅计算 King/FHT 超声场，并绘制：
% (1) xOz: |p|
% (2) xOz: SPL
%
% 说明：
% - 不计算 DIM
% - 不画 xOy 图
% - 画图时只对“显示网格”做抽样降采样，不做插值
% - 所有图都直接使用原始计算结果中的采样点
% - 在最终的 SPL 2D 图上标出 3 个预设绿点
% - 图像按真实坐标等比例显示
% - 同时保存 PNG 和 PDF
%% ============================================================

clear; clc; close all;

%% -------------------- save figures (GLOBAL control) --------------------
global SAVE_FIG SAVE_DIR
SAVE_FIG = true;   % true 保存；false 不保存

%% -------------------- display downsampling --------------------
plot_cfg = struct();
plot_cfg.max_nx_xoz = 400;   % xOz 图最多保留的 x 方向点数
plot_cfg.max_nz_xoz = 400;   % xOz 图最多保留的 z 方向点数

%% -------------------- preset marker points --------------------
% 每行: [x, z]
mark_style.size = 60;
mark_style.face_color = [0 1 0];
mark_style.edge_color = [0 0.5 0];
mark_style.line_width = 0.8;

%% -------------------- medium --------------------
medium.c0 = 343;
medium.rho0 = 1.21;
medium.beta = 1.2;
medium.pref = 2e-5;
medium.use_absorp = true;
medium.atten_handle = @(f) AbsorpAttenCoef(f);

%% -------------------- source --------------------
source.profile = 'Vortex-m'; % 'Uniform' | 'Focus' | 'Vortex-m' | 'Poly' | 'Custom'
source.a = 0.05;
source.v0 = 0.108;
source.v_ratio = 1;
source.m = 3;
source.F = 0.2;

source.f1 = 40e3;
source.fa = 0.5e3;

if source.m == 0
    mark_pts_xz = [ ...
        0.0, 0.3;
        0.0, 1.0;
        0.3, 0.3];
end
if source.m == 1
    mark_pts_xz = [ ...
        0.02, 0.3;
        0.06, 1.0;
        0.3, 0.3];
end
if source.m == 3
    mark_pts_xz = [ ...
        0.05, 0.3;
        0.15, 1.0;
        0.3, 0.3];
end

%% -------------------- calc parameters --------------------
calc = struct();

% --- FHT / King ---
calc.fht.N_FHT = 16384;
calc.fht.rho_max = 1;
calc.fht.Nh_scale = 1.2;
if source.m == 0
    calc.fht.NH_scale = 1.2;
elseif source.m == 1
    calc.fht.NH_scale = 4;
elseif source.m == 3
    calc.fht.NH_scale = 4;
else
    calc.fht.NH_scale = 4;
end
calc.fht.Nh_v_scale = 1.1;
calc.fht.delta = medium.c0 / source.f1 / 4;
calc.fht.zu_max = 2 + calc.fht.delta;
calc.fht.za_max = 0.5;

% --- DIM fields kept only for compatibility, not used here ---
calc.dim.use_freq = 'f2';
calc.dim.dis_coe = 16;
calc.dim.margin = 1;
calc.dim.src_discretization = 'polar';

% --- King analytic spectrum stability ---
calc.king.gspec_method = 'transform';
calc.king.eps_kzz = 1e-3;
calc.king.eps_phase = calc.king.eps_kzz;
calc.king.kz_min = 1e-12;
calc.king.band_refine.enable = false;

% --- ASM settings (unused here) ---
calc.asm.pad_factor = 16;
calc.asm.kzz_eps = 1e-12;

%% -------------------- save setup --------------------
if SAVE_FIG
    save_root = 'JASA_PlotUltra2D';
    m_str = sprintf('m%d', source.m);
    SAVE_DIR = fullfile(save_root, m_str);

    if ~exist(SAVE_DIR, 'dir')
        mkdir(SAVE_DIR);
    end

    local_write_runinfo_txt_kingonly(SAVE_DIR, medium, source, calc, plot_cfg, mark_pts_xz);
else
    SAVE_DIR = '';
end

%% ==================== PROFILER: helpers ====================
get_mem_mb = @() local_get_mem_mb();
fmt_mem = @(x) local_fmt_mem(x);

%% ============================================================
% KING ONLY
%% ============================================================
fprintf('\n==================== KING ONLY ====================\n');

mem0 = get_mem_mb();
t0 = tic;

obs_grid = struct();
res = calc_ultrasound_field(source, medium, calc, obs_grid, 'king');

tAll = toc(t0);
mem1 = get_mem_mb();

fprintf('Total time (king): %.3f s\n', tAll);
fprintf('Total memory: start %s, end %s, delta %s\n', ...
    fmt_mem(mem0), fmt_mem(mem1), fmt_mem(mem1-mem0));

%% ============================================================
% Extract King field
%% ============================================================
p_rz = res.king.p_f1;       % Nr x Nz
z_ultra = res.king.z(:).';  % 1 x Nz
r_ultra = res.king.rho(:);  % Nr x 1

fprintf('King field size: [%d, %d]\n', size(p_rz,1), size(p_rz,2));

%% ============================================================
% xOz plots only
% - plotting uses original computed samples only, with display downsampling
%% ============================================================
rho_max = calc.fht.rho_max;
idx_r = find(r_ultra <= rho_max, 1, 'last');
if isempty(idx_r)
    idx_r = numel(r_ultra);
end

r_use = r_ultra(1:idx_r);
p_use = p_rz(1:idx_r, :);

% mirror to xOz using original samples
x_axis_full = [flipud(-r_use); r_use];
p_xz_full   = [flipud(p_use);  p_use];

% display downsampling only
idx_x_ds = local_make_ds_index(numel(x_axis_full), plot_cfg.max_nx_xoz);
idx_z_ds = local_make_ds_index(numel(z_ultra),    plot_cfg.max_nz_xoz);

x_axis = x_axis_full(idx_x_ds);
z_plot = z_ultra(idx_z_ds);
p_xz   = p_xz_full(idx_x_ds, idx_z_ds);

AMP_xz = abs(p_xz);
SPL_xz = 20*log10(AMP_xz / medium.pref / sqrt(2));

%% -------------------- Figure 1: xOz AMP --------------------
fig1 = figure('Name','King: xOz AMP (|p|)', 'Position', [80 80 900 700], 'Color', 'w');
pcolor(z_plot, x_axis, AMP_xz);
shading flat;
colormap(MyColor('vik'));
cb = colorbar;
cb.Title.Interpreter = 'latex';
cb.Title.String = '$|p|$ (Pa)';
set(cb, 'FontSize', 20);

xlim([0 z_ultra(end)]);
ylim([-rho_max rho_max]);
xlabel('$z$ (m)', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$x$ (m)', 'Interpreter', 'latex', 'FontSize', 22);
daspect([1 1 1]);
set(gca, 'LineWidth', 2, 'FontSize', 20, 'TickLabelInterpreter', 'latex');

local_save_fig_png_pdf(fig1, 'King_xOz_AMP');

%% -------------------- Figure 2: xOz SPL --------------------
fig2 = figure('Name','King: xOz SPL', 'Position', [120 100 900 700], 'Color', 'w');
pcolor(z_plot, x_axis, SPL_xz);
shading flat;
colormap(MyColor('vik'));
cb = colorbar;
cb.Title.Interpreter = 'none';
cb.Title.String = 'SPL (dB)';
set(cb, 'FontSize', 20);
clim([80 130]);

xlim([0 z_ultra(end)]);
ylim([-rho_max rho_max]);
xlabel('$z$ (m)', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$x$ (m)', 'Interpreter', 'latex', 'FontSize', 22);
daspect([1 1 1]);
set(gca, 'LineWidth', 2, 'FontSize', 20, 'TickLabelInterpreter', 'latex');

hold on;
local_plot_green_markers(mark_pts_xz, [0 z_ultra(end)], [-rho_max rho_max], mark_style);
hold off;

local_save_fig_png_pdf(fig2, 'King_xOz_SPL');

%% ============================================================
% Summary
%% ============================================================
fprintf('\n==================== SUMMARY ====================\n');
fprintf('Only King/FHT computed.\n');
fprintf('Total time: %.3f s\n', tAll);
fprintf('xOz display size = %d x %d\n', numel(x_axis), numel(z_plot));
fprintf('Marked %d preset points on SPL figure.\n', size(mark_pts_xz,1));
fprintf('NOTE: memory readings may be NaN if OS API is unavailable.\n');

%% ==================== local functions ====================
function idx = local_make_ds_index(n, n_keep)
if n <= n_keep
    idx = 1:n;
else
    idx = unique(round(linspace(1, n, n_keep)));
end
idx = idx(:);
end

function local_plot_green_markers(mark_pts_xz, xlim_use, ylim_use, mark_style)
% mark_pts_xz: [x, z]
for k = 1:size(mark_pts_xz,1)
    x_now = mark_pts_xz(k,1);
    z_now = mark_pts_xz(k,2);

    if z_now < xlim_use(1) || z_now > xlim_use(2)
        continue;
    end
    if x_now < ylim_use(1) || x_now > ylim_use(2)
        continue;
    end

    scatter(z_now, x_now, ...
        mark_style.size, ...
        'o', ...
        'MarkerFaceColor', mark_style.face_color, ...
        'MarkerEdgeColor', mark_style.edge_color, ...
        'LineWidth', mark_style.line_width);
end
end

function mem_mb = local_get_mem_mb()
mem_mb = NaN;

if ispc
    try
        m = memory();
        mem_mb = double(m.MemUsedMATLAB) / (1024^2);
        return;
    catch
    end
end

try
    rt = java.lang.Runtime.getRuntime();
    used = double(rt.totalMemory() - rt.freeMemory());
    mem_mb = used / (1024^2);
catch
    mem_mb = NaN;
end
end

function s = local_fmt_mem(x)
if isnan(x)
    s = 'NaN';
else
    s = sprintf('%.1f MB', x);
end
end

function local_save_fig_png_pdf(fig_handle, base_name)
global SAVE_FIG SAVE_DIR
if isempty(SAVE_FIG) || ~SAVE_FIG
    return;
end
if isempty(SAVE_DIR) || ~exist(SAVE_DIR, 'dir')
    return;
end

fname = local_sanitize_filename(base_name);
fp_png = fullfile(SAVE_DIR, [fname, '.png']);
fp_pdf = fullfile(SAVE_DIR, [fname, '.pdf']);

set(fig_handle, 'Color', 'w');
set(fig_handle, 'PaperPositionMode', 'auto');

try
    exportgraphics(fig_handle, fp_png, 'Resolution', 300);
catch
    print(fig_handle, fp_png, '-dpng', '-r300');
end

try
    exportgraphics(fig_handle, fp_pdf, 'ContentType', 'vector');
catch
    print(fig_handle, fp_pdf, '-dpdf', '-painters');
end
end

function s = local_sanitize_filename(s)
s = char(string(s));
s = strrep(s, ' ', '_');
bad = '<>:"/\|?*';
for k = 1:numel(bad)
    s = strrep(s, bad(k), '_');
end
end

function local_write_runinfo_txt_kingonly(save_dir, medium, source, calc, plot_cfg, mark_pts_xz)
fp = fullfile(save_dir, 'run_info.txt');
fid = fopen(fp, 'w');
if fid < 0
    warning('Cannot create run_info.txt in %s', save_dir);
    return;
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, '===== RUN INFO (KING ONLY) =====\n');
fprintf(fid, 'Time: %s\n\n', datestr(datetime('now'), 'yyyy-mm-dd HH:MM:SS'));

fprintf(fid, 'medium.c0 = %.15g\n', medium.c0);
fprintf(fid, 'medium.rho0 = %.15g\n', medium.rho0);
fprintf(fid, 'medium.beta = %.15g\n', medium.beta);
fprintf(fid, 'medium.pref = %.15g\n', medium.pref);
fprintf(fid, 'medium.use_absorp = %d\n\n', medium.use_absorp);

fprintf(fid, 'source.profile = %s\n', char(string(source.profile)));
fprintf(fid, 'source.a = %.15g\n', source.a);
fprintf(fid, 'source.v0 = %.15g\n', source.v0);
fprintf(fid, 'source.v_ratio = %.15g\n', source.v_ratio);
fprintf(fid, 'source.m = %d\n', source.m);
fprintf(fid, 'source.F = %.15g\n', source.F);
fprintf(fid, 'source.f1 = %.15g\n', source.f1);
fprintf(fid, 'source.fa = %.15g\n\n', source.fa);

fprintf(fid, 'calc.fht.N_FHT = %.15g\n', calc.fht.N_FHT);
fprintf(fid, 'calc.fht.rho_max = %.15g\n', calc.fht.rho_max);
fprintf(fid, 'calc.fht.Nh_scale = %.15g\n', calc.fht.Nh_scale);
fprintf(fid, 'calc.fht.NH_scale = %.15g\n', calc.fht.NH_scale);
fprintf(fid, 'calc.fht.Nh_v_scale = %.15g\n', calc.fht.Nh_v_scale);
fprintf(fid, 'calc.fht.zu_max = %.15g\n', calc.fht.zu_max);
fprintf(fid, 'calc.fht.za_max = %.15g\n', calc.fht.za_max);
fprintf(fid, 'calc.fht.delta = %.15g\n\n', calc.fht.delta);

fprintf(fid, 'calc.king.gspec_method = %s\n', char(string(calc.king.gspec_method)));
fprintf(fid, 'calc.king.eps_kzz = %.15g\n', calc.king.eps_kzz);
fprintf(fid, 'calc.king.kz_min = %.15g\n', calc.king.kz_min);
fprintf(fid, 'calc.king.band_refine.enable = %d\n\n', calc.king.band_refine.enable);

fprintf(fid, 'plot_cfg.max_nx_xoz = %d\n', plot_cfg.max_nx_xoz);
fprintf(fid, 'plot_cfg.max_nz_xoz = %d\n\n', plot_cfg.max_nz_xoz);

fprintf(fid, 'mark_pts_xz = [\n');
for k = 1:size(mark_pts_xz,1)
    fprintf(fid, '    %.15g, %.15g;\n', mark_pts_xz(k,1), mark_pts_xz(k,2));
end
fprintf(fid, '];\n');

fprintf(fid, '\n===== END =====\n');
end