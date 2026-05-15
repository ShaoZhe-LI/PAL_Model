%% ============================================================
% Transform-only calculation of 2D audio field
%
% 功能：
%   1) 仅使用变换法计算 PAL 音频声场
%   2) 计算二维 rho-z 平面上的音频结果
%   3) 绘制 2D SPL 图像
%   4) 绘图用数据自动降采样，每个维度最多数百点
%
% 外部依赖：
%   - AbsorpAttenCoef.m
%   - solve_kappa0.m
%   - m_FHT.m
%
% 不再包含：
%   - direct integration
%   - direct ultrasound cache
%   - local-effect correction
%   - pressure / velocity 点值重算
%% ============================================================

clear; clc; close all;

%% ===================== case 设置 =====================
% i = 1: m1 = 0, m2 = 0
% i = 2: m1 = 3, m2 = 3
% i = 3: m1 = 0, m2 = 3
% i = 4: m1 = 2, m2 = 5

for i = 2:2

clearvars -except i
close all;

%% ===================== 保存与显示 =====================
save_figures      = true;
save_calc_results = true;
save_params_txt   = true;
show_figures      = true;

%% ===================== 基本物理参数 =====================
a     = 0.05;
v0    = 0.108;
c     = 343;
rho0  = 1.21;
beta  = 1.2;
pref  = 2e-5;

fu = 40e3;
fa = 0.5e3;
f1 = fu;
f2 = fu + fa;

if i == 1
    m1 = 0;
    m2 = 0;
elseif i == 2
    m1 = 3;
    m2 = 3;
elseif i == 3
    m1 = 0;
    m2 = 3;
elseif i == 4
    m1 = 2;
    m2 = 5;
else
    error('Unknown case index i.');
end

ma = m2 - m1;

fprintf('\n============================================================\n');
fprintf('Transform-only audio field calculation\n');
fprintf('case %d: m1 = %d, m2 = %d, ma = %d\n', i, m1, m2, ma);
fprintf('============================================================\n');

%% ===================== FHT / 计算参数 =====================
N_FHT = 8192 * 1;
delta = 0.002 / 1;

rho_max = 0.5;
zu_max  = 15.0;
za_max  = 3.0 + delta;

green_R_min = 1e-12;

%% ===================== 绘图降采样参数 =====================
max_plot_rho_pts = 500;
max_plot_z_pts   = 600;

plot_z_max   = za_max - delta;
plot_rho_max = rho_max;

use_normalized_spl = false;
% false: 绝对 SPL
% true : 归一化 SPL，最大值为 0 dB

plot_dynamic_range_dB = 60;
% 仅用于 normalized SPL 时，例如显示 [-60, 0] dB

%% ===================== 保存路径 =====================
time_tag = datestr(now, 'mmdd_HHMMSS');

save_root_parent = 'result_transform_only_2d';
save_root = fullfile(save_root_parent, ...
    sprintf('AudioTransform2D_%s_m1_%d_m2_%d', time_tag, m1, m2));

if (save_figures || save_calc_results || save_params_txt) && ~exist(save_root, 'dir')
    mkdir(save_root);
end

%% ===================== 构造参数结构 =====================
source_cfg = build_source_cfg(a, v0, m1, m2, f1, fa, f2);
medium_cfg = build_medium_cfg(c, rho0, beta, pref);

%% ===================== 保存参数 =====================
if save_params_txt
    meta = struct();
    meta.a = a;
    meta.v0 = v0;
    meta.c = c;
    meta.rho0 = rho0;
    meta.beta = beta;
    meta.pref = pref;

    meta.fu = fu;
    meta.fa = fa;
    meta.f1 = f1;
    meta.f2 = f2;

    meta.m1 = m1;
    meta.m2 = m2;
    meta.ma = ma;

    meta.N_FHT = N_FHT;
    meta.delta = delta;
    meta.rho_max = rho_max;
    meta.zu_max = zu_max;
    meta.za_max = za_max;
    meta.green_R_min = green_R_min;

    meta.max_plot_rho_pts = max_plot_rho_pts;
    meta.max_plot_z_pts = max_plot_z_pts;
    meta.plot_z_max = plot_z_max;
    meta.plot_rho_max = plot_rho_max;
    meta.use_normalized_spl = use_normalized_spl;
    meta.plot_dynamic_range_dB = plot_dynamic_range_dB;

    meta.save_root = save_root;
    meta.time_tag = time_tag;

    local_write_all_params_txt(fullfile(save_root, 'all_parameters.txt'), meta);
end

%% ============================================================
% 变换法计算 2D audio field
%% ============================================================
t_all = tic;

out = compute_audio_field_transform_only( ...
    N_FHT, delta, ...
    rho_max, zu_max, za_max, ...
    source_cfg, medium_cfg, ...
    green_R_min);

time_transform = toc(t_all);

fprintf('\nTransform-only calculation finished. elapsed = %.2f s\n', time_transform);

rho_audio = out.rho_audio_grid(:);
z_audio   = out.z_audio_grid(:);
pa_audio  = out.pa_transform_full;

spl_audio = local_pressure_to_spl(abs(pa_audio), pref);

if use_normalized_spl
    spl_plot_full = 20*log10(abs(pa_audio) ./ max(abs(pa_audio(:)), [], 'omitnan') + eps);
else
    spl_plot_full = spl_audio;
end

%% ============================================================
% 绘图降采样
%% ============================================================
idx_rho_show = find(rho_audio <= plot_rho_max);
idx_z_show   = find(z_audio   <= plot_z_max);

rho_show_full = rho_audio(idx_rho_show);
z_show_full   = z_audio(idx_z_show);

spl_show_full = spl_plot_full(idx_rho_show, idx_z_show);

ds_rho = max(1, ceil(numel(rho_show_full) / max_plot_rho_pts));
ds_z   = max(1, ceil(numel(z_show_full)   / max_plot_z_pts));

idx_rho_plot = 1:ds_rho:numel(rho_show_full);
idx_z_plot   = 1:ds_z:numel(z_show_full);

rho_plot = rho_show_full(idx_rho_plot);
z_plot   = z_show_full(idx_z_plot);

spl_plot = spl_show_full(idx_rho_plot, idx_z_plot);

fprintf('Plot downsampling:\n');
fprintf('  rho: %d -> %d points, ds = %d\n', ...
    numel(rho_show_full), numel(rho_plot), ds_rho);
fprintf('  z  : %d -> %d points, ds = %d\n', ...
    numel(z_show_full), numel(z_plot), ds_z);

%% ============================================================
% 绘制 2D 音频 SPL 图像
%% ============================================================
fig1 = figure('Color','w', 'Position', [80 80 900 520]);

imagesc(z_plot, rho_plot, spl_plot);
axis xy;
xlabel('z (m)');
ylabel('\rho (m)');

if use_normalized_spl
    title(sprintf('Transform-only normalized audio field, m_1 = %d, m_2 = %d', m1, m2));
    clim([-plot_dynamic_range_dB, 0]);
    cb = colorbar;
    ylabel(cb, 'Normalized SPL (dB)');
else
    title(sprintf('Transform-only audio SPL, m_1 = %d, m_2 = %d', m1, m2));
    cb = colorbar;
    ylabel(cb, 'SPL (dB)');
end

set(gca, 'FontSize', 12, 'LineWidth', 0.9, 'Box', 'on');
colormap turbo;

%% ============================================================
% 另存一张归一化图，方便看指向性结构
%% ============================================================
spl_norm_full = 20*log10(abs(pa_audio) ./ max(abs(pa_audio(:)), [], 'omitnan') + eps);
spl_norm_show = spl_norm_full(idx_rho_show, idx_z_show);
spl_norm_plot = spl_norm_show(idx_rho_plot, idx_z_plot);

fig2 = figure('Color','w', 'Position', [120 120 900 520]);

imagesc(z_plot, rho_plot, spl_norm_plot);
axis xy;
xlabel('z (m)');
ylabel('\rho (m)');
title(sprintf('Transform-only normalized audio field, m_1 = %d, m_2 = %d', m1, m2));
clim([-plot_dynamic_range_dB, 0]);
cb = colorbar;
ylabel(cb, 'Normalized SPL (dB)');

set(gca, 'FontSize', 12, 'LineWidth', 0.9, 'Box', 'on');
colormap turbo;

%% ============================================================
% 保存结果
%% ============================================================
if save_figures
    saveas(fig1, fullfile(save_root, 'Audio2D_SPL_absolute.png'));
    savefig(fig1, fullfile(save_root, 'Audio2D_SPL_absolute.fig'));

    saveas(fig2, fullfile(save_root, 'Audio2D_SPL_normalized.png'));
    savefig(fig2, fullfile(save_root, 'Audio2D_SPL_normalized.fig'));
end

if ~show_figures
    close(fig1);
    close(fig2);
end

if save_calc_results
    save(fullfile(save_root, 'audio_transform_only_2d_results.mat'), ...
        'out', ...
        'rho_audio', 'z_audio', 'pa_audio', 'spl_audio', ...
        'rho_plot', 'z_plot', 'spl_plot', 'spl_norm_plot', ...
        'm1', 'm2', 'ma', ...
        'N_FHT', 'delta', 'rho_max', 'zu_max', 'za_max', ...
        'time_transform', ...
        '-v7.3');
end

fprintf('\nAll done.\n');
fprintf('Results folder: %s\n', save_root);

end

%% ============================================================
% 变换法计算完整 rho-z 音频场
%% ============================================================
function out = compute_audio_field_transform_only( ...
    N_FHT, delta, ...
    rho_max, zu_max, za_max, ...
    source_cfg, medium_cfg, ...
    green_R_min)

c     = medium_cfg.c0;
rho0  = medium_cfg.rho0;
beta  = medium_cfg.beta;

f1 = source_cfg.f1;
f2 = source_cfg.f2;
fa = source_cfg.fa;

m1 = source_cfg.m1;
m2 = source_cfg.m2;
ma = m2 - m1;

w1 = 2*pi*f1;
w2 = 2*pi*f2;
wa = 2*pi*fa;

k1 = w1/c + 1j*AbsorpAttenCoef(f1);
k2 = w2/c + 1j*AbsorpAttenCoef(f2);
ka = wa/c + 1j*AbsorpAttenCoef(fa);

Nh = 2 * rho_max;
NH = 4 * w2 / c;

n_FHT = 0:N_FHT-1;
[a_solve, k0, x1, x0] = solve_kappa0(N_FHT, n_FHT);

rho_audio = (x1 * Nh).';
z_ultra   = 0:delta:zu_max;
z_audio   = 0:delta:za_max;

Nz = numel(z_ultra);

fprintf('\n[FHT grid]\n');
fprintf('  N_FHT = %d\n', N_FHT);
fprintf('  delta = %.6g m\n', delta);
fprintf('  rho points = %d\n', numel(rho_audio));
fprintf('  z_ultra points = %d\n', numel(z_ultra));
fprintf('  z_audio points = %d\n', numel(z_audio));

%% ===================== 源面速度谱 =====================
fprintf('\nBuilding source velocity spectra...\n');

a  = source_cfg.a;
v0 = source_cfg.v0;

Nh_v = 1.1 * a;
NH_v = NH;
rho_v_grid = (x1 * Nh_v).';

switch source_cfg.profile
    case 'Uniform'
        vs1 = v0 * double(rho_v_grid <= a);
        vs2 = vs1;

    case 'Vortex-m'
        vs1 = v0 * double(rho_v_grid <= a);
        vs2 = vs1;

    case 'Focus'
        F = source_cfg.F;
        vs1 = v0 * exp(-1j*real(k1)*sqrt(rho_v_grid.^2 + F^2)) ...
            .* double(rho_v_grid <= a);
        vs2 = v0 * exp(-1j*real(k2)*sqrt(rho_v_grid.^2 + F^2)) ...
            .* double(rho_v_grid <= a);

    otherwise
        error('Unknown source profile: %s', source_cfg.profile);
end

Vs1 = m_FHT(vs1, N_FHT, 1, Nh_v, NH_v, a_solve, x0, x1, k0, m1);
Vs2 = m_FHT(vs2, N_FHT, 1, Nh_v, NH_v, a_solve, x0, x1, k0, m2);

%% ===================== 超声 Green 变换 =====================
fprintf('\nComputing ultrasonic Green transforms...\n');

G1_raw = build_green_space_g_transform_raw( ...
    rho_audio, z_ultra, k1, green_R_min, ...
    N_FHT, Nh, NH, a_solve, x0, x1, k0);

G2_raw = build_green_space_g_transform_raw( ...
    rho_audio, z_ultra, k2, green_R_min, ...
    N_FHT, Nh, NH, a_solve, x0, x1, k0);

G1 = (-4*pi*1j) * G1_raw;
G2 = (-4*pi*1j) * G2_raw;

%% ===================== 超声压力场 =====================
fprintf('\nComputing ultrasonic pressure fields by transform method...\n');

p1_transform = compute_ultrasound_pressure_from_G( ...
    G1, Vs1, N_FHT, Nz, NH, Nh, a_solve, x0, x1, k0, m1, rho0, c, k1);

p2_transform = compute_ultrasound_pressure_from_G( ...
    G2, Vs2, N_FHT, Nz, NH, Nh, a_solve, x0, x1, k0, m2, rho0, c, k2);

%% ===================== 非线性源项 q =====================
fprintf('\nComputing nonlinear source term q...\n');

q_transform = conj(p1_transform) .* p2_transform ...
    * beta*wa/(1j*rho0^2*c^4);

mask_rho = double(rho_audio(:) <= rho_max);
q_transform = q_transform .* repmat(mask_rho, 1, size(q_transform, 2));

%% ===================== 音频 Green 变换 =====================
fprintf('\nComputing audio Green transform...\n');

Ga_raw = build_green_space_g_transform_raw( ...
    rho_audio, z_ultra, ka, green_R_min, ...
    N_FHT, Nh, NH, a_solve, x0, x1, k0);

Ga = (-4*pi*1j) * Ga_raw;

%% ===================== 音频场 =====================
fprintf('\nComputing audio field by transform method...\n');

[pa_transform_full, phia_transform_full] = compute_paW_from_Gr00( ...
    q_transform, z_ultra, z_audio, Ga, ...
    N_FHT, Nh, NH, a_solve, x0, x1, k0, ma, ...
    delta, rho0, wa);

%% ===================== 输出 =====================
out = struct();

out.rho_audio_grid = rho_audio(:);
out.z_ultra_grid   = z_ultra(:);
out.z_audio_grid   = z_audio(:);

out.p1_transform = p1_transform;
out.p2_transform = p2_transform;
out.q_transform  = q_transform;

out.pa_transform_full   = pa_transform_full;
out.phia_transform_full = phia_transform_full;

out.k1 = k1;
out.k2 = k2;
out.ka = ka;

out.w1 = w1;
out.w2 = w2;
out.wa = wa;

out.m1 = m1;
out.m2 = m2;
out.ma = ma;
end

%% ============================================================
% 用给定谱域 G 计算超声压力
%% ============================================================
function p_out = compute_ultrasound_pressure_from_G( ...
    G_spec, Vs, N_FHT, Nz, NH, Nh, a_solve, x0, x1, k0, m_use, rho0, c0, k_use)

F = G_spec .* Vs;

phi = -1j * m_FHT( ...
    F, N_FHT, Nz, ...
    NH, Nh, ...
    a_solve, x0, x1, k0, ...
    m_use);

p_out = 1j * rho0 * c0 * real(k_use) .* phi;
end

%% ============================================================
% 空间域 Green function 做 0 阶 Hankel 变换
%% ============================================================
function G_raw = build_green_space_g_transform_raw( ...
    rho_vec, z_vec, k_use, green_R_min, ...
    N_FHT, Nh, NH, a_solve, x0, x1, k0)

[RHO, Z] = ndgrid(rho_vec, z_vec);

RR = sqrt(RHO.^2 + Z.^2);
RR_use = max(RR, green_R_min);

g_space = exp(1j * k_use * RR_use) ./ (4*pi * RR_use);

G_raw = m_FHT( ...
    g_space, ...
    N_FHT, numel(z_vec), ...
    Nh, NH, ...
    a_solve, x0, x1, k0, ...
    0);
end

%% ============================================================
% 从给定 Gr00 计算音频 pa_W
%% ============================================================
function [pa_W, phia_W] = compute_paW_from_Gr00( ...
    q_full, z, z_audio, Gr00, ...
    N_FHT, Nh, NH, a_solve, x0, x1, k0, ma, ...
    delta, rho0, wa)

Nz  = numel(z);
Nza = numel(z_audio);

absz = [-fliplr(z(2:end)) z];
Nz1  = numel(absz);
N_conv = Nz1 + Nza - 1;

Qr00 = m_FHT( ...
    q_full, ...
    N_FHT, Nz, ...
    Nh, NH, ...
    a_solve, x0, x1, k0, ...
    ma);

Qr0 = [fliplr(Qr00(:,2:end)) Qr00];

[M_0, N_0] = size(Qr0);
Qr = [Qr0, zeros(M_0, N_conv-N_0)];
Q  = (fft(Qr.')).';

Gr0 = [fliplr(Gr00(:,2:end)) Gr00];
Gr  = [Gr0, zeros(M_0, N_conv-N_0)];
G   = (fft(Gr.')).';

Pa = Q .* G;

par0 = (ifft(Pa.')).';
par  = par0(:, N_conv-Nza+1:N_conv);

phia0 = m_FHT( ...
    par, ...
    N_FHT, Nza, ...
    NH, Nh, ...
    a_solve, x0, x1, k0, ...
    ma);

phia_W = -phia0 * delta * 1j / 2;
pa_W   = 1j * rho0 * wa * phia_W;
end

%% ============================================================
% 构造 source cfg
%% ============================================================
function source_cfg = build_source_cfg(a, v0, m1, m2, f1, fa, f2)

source_cfg = struct();

source_cfg.profile = 'Vortex-m';

source_cfg.a = a;
source_cfg.v0 = v0;
source_cfg.v_ratio = 1;

source_cfg.m1 = m1;
source_cfg.m2 = m2;
source_cfg.m  = m1;

source_cfg.F = 0.2;

source_cfg.f1 = f1;
source_cfg.fa = fa;
source_cfg.f2 = f2;

source_cfg.internal = struct();
end

%% ============================================================
% 构造 medium cfg
%% ============================================================
function medium_cfg = build_medium_cfg(c, rho0, beta, pref)

medium_cfg = struct();

medium_cfg.c0 = c;
medium_cfg.rho0 = rho0;
medium_cfg.beta = beta;
medium_cfg.pref = pref;

medium_cfg.use_absorp = true;
medium_cfg.atten_handle = @(f) AbsorpAttenCoef(f);

medium_cfg.internal = struct();
end

%% ============================================================
% 压力幅值转 SPL
%% ============================================================
function spl = local_pressure_to_spl(p_amp, pref)

spl = 20 * log10(p_amp ./ max(pref, eps) ./ sqrt(2) + eps);
end

%% ============================================================
% 保存全部参数到 txt
%% ============================================================
function local_write_all_params_txt(txt_path, S)

fid = fopen(txt_path, 'w');

if fid < 0
    error('Cannot open txt file for writing: %s', txt_path);
end

cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, '============================================================\n');
fprintf(fid, 'ALL PARAMETERS\n');
fprintf(fid, 'Generated time: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, '============================================================\n\n');

local_dump_any(fid, 'params', S, 0);
end

%% ============================================================
% 递归写任意变量
%% ============================================================
function local_dump_any(fid, name, val, indent)

sp = repmat(' ', 1, indent);

if isstruct(val)
    fprintf(fid, '%s%s = struct\n', sp, name);
    fn = fieldnames(val);
    for ii = 1:numel(fn)
        local_dump_any(fid, sprintf('%s.%s', name, fn{ii}), val.(fn{ii}), indent + 2);
    end
    return;
end

if islogical(val) && isscalar(val)
    fprintf(fid, '%s%s = %s\n', sp, name, mat2str(val));
    return;
end

if isnumeric(val)
    if isscalar(val)
        if isreal(val)
            fprintf(fid, '%s%s = %.16g\n', sp, name, val);
        else
            fprintf(fid, '%s%s = %.16g%+.16gi\n', sp, name, real(val), imag(val));
        end
    else
        sz = size(val);
        fprintf(fid, '%s%s : numeric array, size = [%s]\n', sp, name, num2str(sz));

        vec = val(:);
        nshow = min(numel(vec), 20);

        if isreal(vec)
            fprintf(fid, '%s  first %d values = [', sp, nshow);
            fprintf(fid, ' %.8g', vec(1:nshow));
            fprintf(fid, ' ]\n');
        else
            fprintf(fid, '%s  first %d values = [', sp, nshow);
            for kk = 1:nshow
                fprintf(fid, ' %.8g%+.8gi', real(vec(kk)), imag(vec(kk)));
            end
            fprintf(fid, ' ]\n');
        end
    end
    return;
end

if ischar(val) || (isstring(val) && isscalar(val))
    fprintf(fid, '%s%s = %s\n', sp, name, char(string(val)));
    return;
end

if isa(val, 'function_handle')
    fprintf(fid, '%s%s = %s\n', sp, name, func2str(val));
    return;
end

if iscell(val)
    sz = size(val);
    fprintf(fid, '%s%s : cell, size = [%s]\n', sp, name, num2str(sz));

    nshow = min(numel(val), 10);
    for ii = 1:nshow
        local_dump_any(fid, sprintf('%s{%d}', name, ii), val{ii}, indent + 2);
    end
    return;
end

try
    fprintf(fid, '%s%s = %s\n', sp, name, evalc('disp(val)'));
catch
    fprintf(fid, '%s%s = <unprintable type: %s>\n', sp, name, class(val));
end
end