%% ============================================================
% MAIN: King vs DIM (single call: method='both') + profiler (time & memory)
% + EXTRA FIGS (FHT only):
% (1) xOz (from King/FHT): |p| and SPL on (x,z)
% (2) xOy @ z≈1 m (from King/FHT): SPL and phase on (x,y)
% via planar phase extension: p(x,y)=p(rho)*exp(j*m*phi)
% ============================================================

clear; clc; close all;


%% -------------------- save figures (GLOBAL control) --------------------
global SAVE_PNG SAVE_DIR
SAVE_PNG = false;   % <<< 全局开关：true 保存；false 不保存

%% -------------------- medium --------------------
medium.c0 = 343;
medium.rho0 = 1.21;
medium.beta = 1.2;
medium.pref = 2e-5;
medium.use_absorp = true;
medium.atten_handle = @(f) AbsorpAttenCoef(f);

%% -------------------- source --------------------
source.profile = 'Vortex-m'; % 'Uniform' | 'Focus' | 'Vortex-m' | 'Poly' | 'Custom'
source.a = 0.1;
source.v0 = 0.172;
source.v_ratio = 1;
source.m = 10;
source.F = 0.2;
% source.poly_n = 0;

source.f1 = 42e3;
source.fa = 2e3;

%% -------------------- calc parameters --------------------
calc = struct();

% --- FHT / King ---
calc.fht.N_FHT = 32768 * 1;
calc.fht.rho_max = 1.2;
calc.fht.Nh_scale = 1.2;
calc.fht.NH_scale = 4;
calc.fht.Nh_v_scale = 1.1;
calc.fht.zu_max = 1.1;
calc.fht.za_max = 1;

% --- DIM source discretization ---
calc.dim.use_freq = 'f2';
calc.dim.dis_coe = 16 * 1;
calc.dim.margin = 1;
calc.dim.src_discretization = 'polar';   % 'cart' (default) | 'polar'


% --- King analytic spectrum stability ---
calc.king.gspec_method = 'analytic'; % 'analytic' or 'transform'，优选前者
calc.king.eps_kzz = 1e-3; % keep your naming
% map to function-expected fields:
calc.king.eps_phase = calc.king.eps_kzz;
calc.king.kz_min = 1e-12;

calc.king.band_refine.enable = true;
% 本行决定是否局部加细以改善近轴结果
% 计算发现，局部加细比等倍率的全局加细更加耗时，且效果不如后者，后续应改善
% 相应的，局部加细占用内存会少于全局加细

% --- ASM settings (used only if calc.dim.method = 'asm') ---
calc.asm.pad_factor = 16;
calc.asm.kzz_eps = 1e-12;

%% -------------------- fig setup --------------------
fig.unwrap = false;

%% -------------------- save setup --------------------
if SAVE_PNG
    tstr = datestr(datetime('now'), 'mmdd_HHMM');
    a_str = sprintf('%.2fm', source.a);
    m_str = sprintf('m=%d', source.m);
    SAVE_DIR = sprintf('%s_%s__%s', a_str, m_str, tstr);

    if ~exist(SAVE_DIR, 'dir')
        mkdir(SAVE_DIR);
    end

    local_write_runinfo_txt(SAVE_DIR, medium, source, calc, fig);
else
    SAVE_DIR = '';
end

%% -------------------- target plane --------------------
z_target = 1;
ds = 16 * 2;

%% ==================== PROFILER: helpers ====================
get_mem_mb = @() local_get_mem_mb(); % returns NaN if not available
fmt_mem = @(x) local_fmt_mem(x);

%% ============================================================
% PREPARE OBS GRID (need rho/z grids before calling 'both')
% - Use make_source_velocity() to get FHT grids (NO field computation here)
% ============================================================
[~, fht_tmp] = make_source_velocity(source, medium, calc); %#ok<ASGLU>

rho_full = fht_tmp.xh(:);       % Nr x 1
zK_full  = fht_tmp.z_ultra(:);  % Nz x 1

[~, iz] = min(abs(zK_full - z_target));
z_use = zK_full(iz);

rho_ds = rho_full(1:ds:end);

% Build observation grid for DIM
obs_grid.dim.x = rho_ds(:).';
obs_grid.dim.y = 0;
obs_grid.dim.z = z_use;
obs_grid.dim.block_size = 200000;   % used by rayleigh only; ignored by asm

%% ============================================================
% SINGLE CALL: method = 'both' (one res stores king+dim)
%% ============================================================
fprintf('\n==================== BOTH (King + DIM) ====================\n');
mem0 = get_mem_mb(); t0 = tic;

res = calc_ultrasound_field( ...
    source, medium, calc, obs_grid, ...
    'both');

tAll = toc(t0); mem1 = get_mem_mb();

fprintf('Total time (both): %.3f s\n', tAll);
fprintf('Total memory: start %s, end %s, delta %s\n', fmt_mem(mem0), fmt_mem(mem1), fmt_mem(mem1-mem0));

%% -------------------- extract fields (use f1 here) --------------------
% King: take same z index (nearest) and downsample to rho_ds
zK = res.king.z(:);
rhoK = res.king.rho(:);

[~, izK] = min(abs(zK - z_use));
pK_full = res.king.p_f1(:, izK);     % Nr x 1
pK = pK_full(1:ds:end);              % match rho_ds

% DIM:
if strcmpi(res.calc.dim.method,'rayleigh')
    % DIM-Rayleigh output is Ny x Nx x Nz with Ny=1, Nx=numel(rho_ds), Nz=1
    pD = squeeze(res.dim.p_f1(1,:,1)).';  % Nx x 1
else
    % ASM output: p_f1 is Ny x Nx x Nz on native uniform grid
    xA = res.dim.x(:);
    yA = res.dim.y(:);
    zA = res.dim.z(:);

    [~, iy0] = min(abs(yA - 0));
    [~, iz0] = min(abs(zA - z_use));

    pA_line = squeeze(res.dim.p_f1(iy0,:,iz0)).'; % numel(xA) x 1
    % interpolate to rho_ds ONLY in MAIN
    pD = interp1(xA, pA_line, rho_ds, 'linear', 0);
end

%% -------------------- metrics --------------------
magK = abs(pK);
magD = abs(pD);

if fig.unwrap
    phK = unwrap(angle(pK));
    phD = unwrap(angle(pD));
else
    phK = (angle(pK));
    phD = (angle(pD));
end

eps0 = 1e-12;
rel_err_log = log10( abs(pD - pK) ./ (abs(pK) + eps0) );

%% -------------------- plots: 1D compare (rho line) --------------------
figure('Name',sprintf('King vs DIM-%s (z≈1 m)', upper(string(res.calc.dim.method))), ...
    'position',[100 100 1200 1200]);

subplot(4,1,1);
plot(rho_ds, 20*log10(magK / medium.pref / sqrt(2)), 'LineWidth',1.5); hold on;
plot(rho_ds, 20*log10(magD / medium.pref / sqrt(2)), '--', 'LineWidth',1.5);
grid on;
xlabel('\rho (m)'); ylabel('SPL (dB)');
title(sprintf('Magnitude @ z = %.2f m', z_use));
legend('King','DIM');

subplot(4,1,2);
plot(rho_ds, magK, 'LineWidth',1.5); hold on;
plot(rho_ds, magD, '--', 'LineWidth',1.5);
grid on;
xlabel('\rho (m)'); ylabel('|p| (Pa)');
title(sprintf('Magnitude @ z = %.2f m', z_use));
legend('King','DIM');

subplot(4,1,3);
plot(rho_ds, phK, 'LineWidth',1.5); hold on;
plot(rho_ds, phD, '--', 'LineWidth',1.5);
grid on;
xlabel('\rho (m)'); ylabel('Phase (rad)');
title('Phase');
legend('King','DIM');

subplot(4,1,4);
plot(rho_ds, rel_err_log, 'LineWidth',1.5);
grid on;
xlabel('\rho (m)');
ylabel('log_{10} relative error');
title('log_{10}(|p_{DIM}-p_{King}| / |p_{King}|)');

local_save_fig_png(gcf, sprintf('King_vs_DIM_1D_z%.2fm', z_use));

fprintf('\n==================== SUMMARY ====================\n');
fprintf('dim_method = %s\n', string(res.calc.dim.method));
fprintf('z_target = %.1f m, z_use = %.2f m (nearest FHT grid)\n', z_target, z_use);
fprintf('Total (both) time: %.3f s\n', tAll);
fprintf('NOTE: memory readings may be NaN if OS API is unavailable.\n');

%% ============================================================
% 2D FIG #1: xOz (FHT/King only)
% - show |p| (AMP) and SPL on xOz
% ============================================================

% choose which frequency field to plot
p_rz = res.king.p_f1;     % Nr x Nz
z_ultra = res.king.z(:).';% 1 x Nz
r_ultra = res.king.rho(:);% Nr x 1

rho_max = calc.fht.rho_max;
idx_r = find(r_ultra <= rho_max, 1, 'last');
if isempty(idx_r); idx_r = numel(r_ultra); end

r_use = r_ultra(1:idx_r);      % Nr_use x 1
p_use = p_rz(1:idx_r, :);      % Nr_use x Nz

% mirror to xOz (no interpolation)
x_axis = [flipud(-r_use); r_use];   % (2*Nr_use) x 1
p_xz = [flipud(p_use); p_use];      % (2*Nr_use) x Nz

AMP_xz = abs(p_xz);
SPL_xz = 20*log10(AMP_xz / medium.pref / sqrt(2)); % RMS SPL

% ---- AMP (|p|) map ----
figure('Name','FHT: xOz AMP (|p|)');
pcolor(z_ultra, x_axis, AMP_xz);
colormap(MyColor('vik'));
shading interp;
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$|p|$ (Pa)';
set(clb,'Fontsize',20);

xlim([0 z_ultra(end)]); ylim([-rho_max rho_max]);
fontsize(gca,24,'points');
xlabel('$z$ (m)', 'Interpreter','latex','Fontsize',21);
ylabel('$x$ (m)', 'Interpreter','latex','Fontsize',21);
pbaspect([z_ultra(end), 2*rho_max, 1]);
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');

local_save_fig_png(gcf, 'FHT_xOz_AMP');

% ---- SPL map ----
figure('Name','FHT: xOz SPL');
pcolor(z_ultra, x_axis, SPL_xz);
colormap(MyColor('vik'));
shading interp;
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$\mathrm{SPL}$ (dB re $p_{\mathrm{ref}}$, RMS)';
set(clb,'Fontsize',20);

xlim([0 z_ultra(end)]); ylim([-rho_max rho_max]);
fontsize(gca,24,'points');
xlabel('$z$ (m)', 'Interpreter','latex','Fontsize',21);
ylabel('$x$ (m)', 'Interpreter','latex','Fontsize',21);
pbaspect([z_ultra(end), 2*rho_max, 1]);
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');

local_save_fig_png(gcf, 'FHT_xOz_SPL');

%% ============================================================
% 2D FIG #2: xOy @ z≈1 m (FHT vs DIM) — SAME colorbar per figure
% - Construct polar grid (rho,theta) and do phase extension:
%   p2D(r,theta) = p_line(r) * exp(1j*m*theta)
% ============================================================

% ---- view settings ----
r_boundary = 0.30; % radius to show (m)
theta_fig = 0:0.01:2*pi;

% ---- m ----
if isfield(res,'source') && isfield(res.source,'m_used')
    m_use = res.source.m_used;
else
    m_use = source.m;
end

% ---- radial samples (NO interp) ----
rho_fig = rho_ds(:).';          % 1 x Nr_fig
pK_line = pK(:).';              % 1 x Nr_fig
pD_line = pD(:).';              % 1 x Nr_fig

% crop to r_boundary
idx_rb = find(rho_fig <= r_boundary, 1, 'last');
if isempty(idx_rb); idx_rb = 1; end
rho_fig = rho_fig(1:idx_rb);
pK_line = pK_line(1:idx_rb);
pD_line = pD_line(1:idx_rb);

% polar grid -> Cartesian
[TH, R] = meshgrid(theta_fig, rho_fig);
[X, Y] = pol2cart(TH, R);

% phase extension
pK_2D = (pK_line(:) * ones(1, numel(theta_fig))) .* exp(1i*m_use*TH);
pD_2D = (pD_line(:) * ones(1, numel(theta_fig))) .* exp(1i*m_use*TH);

% fields
AMP_K = abs(pK_2D);
AMP_D = abs(pD_2D);

SPL_K = 20*log10(AMP_K / medium.pref / sqrt(2)); % RMS SPL
SPL_D = 20*log10(AMP_D / medium.pref / sqrt(2));

PH_K = angle(pK_2D)/pi;
PH_D = angle(pD_2D)/pi;

% ---- unified color limits ----
spl_lim = [min([SPL_K(:); SPL_D(:)]), max([SPL_K(:); SPL_D(:)])];
amp_lim = [min([AMP_K(:); AMP_D(:)]), max([AMP_K(:); AMP_D(:)])];
ph_lim  = [-1 1];

%% ---- FIG: xOy SPL (two subplots, SAME colorbar range) ----
figure('Name',sprintf('xOy SPL @ z=%.2f m, m=%d (FHT vs DIM)', z_use, m_use), ...
    'position',[100 100 1400 650]);

subplot(1,2,1);
pcolor(X, Y, SPL_K); shading flat
xlim([-1.2*r_boundary 1.2*r_boundary]);
ylim([-1.2*r_boundary 1.2*r_boundary]);
pbaspect([1 1 1])
colormap(MyColor('vik'))
clim(spl_lim);
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$\mathrm{SPL}$ (dB re $p_{\mathrm{ref}}$, RMS)';
set(clb,'Fontsize',18);
set(gca,'position',[0.07 0.12 0.38 0.78]);
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');
xlabel('$x$ (m)','Interpreter','latex','Fontsize',18);
ylabel('$y$ (m)','Interpreter','latex','Fontsize',18);
title('FHT (from King)','Interpreter','latex','Fontsize',20);

subplot(1,2,2);
pcolor(X, Y, SPL_D); shading flat
xlim([-1.2*r_boundary 1.2*r_boundary]);
ylim([-1.2*r_boundary 1.2*r_boundary]);
pbaspect([1 1 1])
colormap(MyColor('vik'))
clim(spl_lim);
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$\mathrm{SPL}$ (dB re $p_{\mathrm{ref}}$, RMS)';
set(clb,'Fontsize',18);
set(gca,'position',[0.55 0.12 0.38 0.78]);
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');
xlabel('$x$ (m)','Interpreter','latex','Fontsize',18);
ylabel('$y$ (m)','Interpreter','latex','Fontsize',18);
title(sprintf('DIM-%s', upper(string(res.calc.dim.method))), 'Interpreter','latex','Fontsize',20);

sgtitle(sprintf('$xOy$ SPL @ $z=%.2f$ m, $m=%d$, $r\\le%.2f$ m', z_use, m_use, r_boundary), ...
    'Interpreter','latex','Fontsize',20);

local_save_fig_png(gcf, sprintf('xOy_SPL_z%.2fm_m%d', z_use, m_use));

%% ---- FIG: xOy AMP (|p|) (two subplots, SAME colorbar range) ----
figure('Name',sprintf('xOy AMP @ z=%.2f m, m=%d (FHT vs DIM)', z_use, m_use), ...
    'position',[100 100 1400 650]);

subplot(1,2,1);
pcolor(X, Y, AMP_K); shading flat
xlim([-1.2*r_boundary 1.2*r_boundary]);
ylim([-1.2*r_boundary 1.2*r_boundary]);
pbaspect([1 1 1])
colormap(MyColor('vik'))
clim(amp_lim);
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$|p|$ (Pa)';
set(clb,'Fontsize',18);
set(gca,'position',[0.07 0.12 0.38 0.78]);
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');
xlabel('$x$ (m)','Interpreter','latex','Fontsize',18);
ylabel('$y$ (m)','Interpreter','latex','Fontsize',18);
title('FHT (from King)','Interpreter','latex','Fontsize',20);

subplot(1,2,2);
pcolor(X, Y, AMP_D); shading flat
xlim([-1.2*r_boundary 1.2*r_boundary]);
ylim([-1.2*r_boundary 1.2*r_boundary]);
pbaspect([1 1 1])
colormap(MyColor('vik'))
clim(amp_lim);
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$|p|$ (Pa)';
set(clb,'Fontsize',18);
set(gca,'position',[0.55 0.12 0.38 0.78]);
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');
xlabel('$x$ (m)','Interpreter','latex','Fontsize',18);
ylabel('$y$ (m)','Interpreter','latex','Fontsize',18);
title(sprintf('DIM-%s', upper(string(res.calc.dim.method))), 'Interpreter','latex','Fontsize',20);

sgtitle(sprintf('$xOy$ AMP ($|p|$) @ $z=%.2f$ m, $m=%d$, $r\\le%.2f$ m', z_use, m_use, r_boundary), ...
    'Interpreter','latex','Fontsize',20);

local_save_fig_png(gcf, sprintf('xOy_AMP_z%.2fm_m%d', z_use, m_use));

%% ---- FIG: xOy Phase/pi (two subplots, SAME colorbar range) ----
figure('Name',sprintf('xOy Phase/pi @ z=%.2f m, m=%d (FHT vs DIM)', z_use, m_use), ...
    'position',[100 100 1400 650]);

subplot(1,2,1);
pcolor(X, Y, PH_K); shading flat
xlim([-1.2*r_boundary 1.2*r_boundary]);
ylim([-1.2*r_boundary 1.2*r_boundary]);
pbaspect([1 1 1])
colormap('hsv')
clim(ph_lim);
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$\angle p/\pi$';
set(clb,'Fontsize',18);
set(gca,'position',[0.07 0.12 0.38 0.78]);
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');
xlabel('$x$ (m)','Interpreter','latex','Fontsize',18);
ylabel('$y$ (m)','Interpreter','latex','Fontsize',18);
title('FHT (from King)','Interpreter','latex','Fontsize',20);

subplot(1,2,2);
pcolor(X, Y, PH_D); shading flat
xlim([-1.2*r_boundary 1.2*r_boundary]);
ylim([-1.2*r_boundary 1.2*r_boundary]);
pbaspect([1 1 1])
colormap('hsv')
clim(ph_lim);
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$\angle p/\pi$';
set(clb,'Fontsize',18);
set(gca,'position',[0.55 0.12 0.38 0.78]);
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');
xlabel('$x$ (m)','Interpreter','latex','Fontsize',18);
ylabel('$y$ (m)','Interpreter','latex','Fontsize',18);
title(sprintf('DIM-%s', upper(string(res.calc.dim.method))), 'Interpreter','latex','Fontsize',20);

sgtitle(sprintf('$xOy$ Phase$/\\pi$ @ $z=%.2f$ m, $m=%d$, $r\\le%.2f$ m', z_use, m_use, r_boundary), ...
    'Interpreter','latex','Fontsize',20);

local_save_fig_png(gcf, sprintf('xOy_Phase_z%.2fm_m%d', z_use, m_use));


%%
xH = res.fht.xH;
N_FHT = res.calc.fht.N_FHT;
r_idx = 1:ds:N_FHT;
G1_ana_plot = res.king.G1(r_idx,iz);
Vs1 = res.king.Vs1;


figure('Name','Radial');
plot(xH(r_idx), abs(G1_ana_plot), 'LineWidth', 1.2);
hold on
yyaxis right
plot(xH(r_idx), abs(Vs1(r_idx)));
hold off
legend('ana', 'V');

local_save_fig_png(gcf, 'Radial_G&V');

%% ===================== EXTRA FIG (normalized + mark <1e-3) =====================
% - Same layout as the first figure (4 subplots)
% - Subplot(2): normalized magnitude (each by its own max)
% - Mark regions where |p|/max < 1e-3
%   * Mark on subplot(2), (3), and (4)
% - Legend includes explicit label for marked region
% - Save with suffix "__norm"

% ---------- normalization ----------
magK = abs(pK);
magD = abs(pD);

magK_n = magK ./ max(magK + eps0);
magD_n = magD ./ max(magD + eps0);

thr = 1e-3;

% ---------- union mask (King OR DIM below threshold) ----------
maskU = (magK_n < thr) | (magD_n < thr);
segU  = local_mask_to_segments(rho_ds, maskU);

% ---------- plot ----------
fig2 = figure('Name',sprintf('King vs DIM-%s (z≈1 m) [NORM]', ...
    upper(string(res.calc.dim.method))), ...
    'position',[150 120 1200 1200]);

% ==================== (1) SPL ====================
subplot(4,1,1);
plot(rho_ds, 20*log10(magK / medium.pref / sqrt(2)), 'LineWidth',1.5); hold on;
plot(rho_ds, 20*log10(magD / medium.pref / sqrt(2)), '--', 'LineWidth',1.5);
grid on;
xlabel('\rho (m)'); ylabel('SPL (dB)');
title(sprintf('Magnitude @ z = %.2f m', z_use));
legend('King','DIM','Location','best');

% ==================== (2) normalized magnitude ====================
subplot(4,1,2);
plot(rho_ds, magK_n, 'LineWidth',1.5); hold on;
plot(rho_ds, magD_n, '--', 'LineWidth',1.5);
grid on;
xlabel('\rho (m)'); ylabel('|p| / max(|p|)');
title(sprintf('Normalized magnitude (thr = %.0e)', thr));

yl = ylim;
hPatch = local_mark_segments_strong(segU, yl);  % <--- mark
ylim(yl);

legend([ ...
    findobj(gca,'Type','Line','-and','LineStyle','-'), ...
    findobj(gca,'Type','Line','-and','LineStyle','--'), ...
    hPatch(1)], ...
    {'King','DIM','|p|/max < 1e-3'}, ...
    'Location','best');

% ==================== (3) phase ====================
subplot(4,1,3);
plot(rho_ds, phK, 'LineWidth',1.5); hold on;
plot(rho_ds, phD, '--', 'LineWidth',1.5);
grid on;
xlabel('\rho (m)'); ylabel('Phase (rad)');
title('Phase');

yl = ylim;
local_mark_segments_strong(segU, yl);
ylim(yl);

legend('King','DIM','Location','best');

% ==================== (4) relative error ====================
subplot(4,1,4);
rel_err_log = log10( abs(pD - pK) ./ (abs(pK) + eps0) );
plot(rho_ds, rel_err_log, 'LineWidth',1.5);
grid on;
xlabel('\rho (m)');
ylabel('log_{10} relative error');
title('log_{10}(|p_{DIM}-p_{King}| / |p_{King}|)');

yl = ylim;
local_mark_segments_strong(segU, yl);
ylim(yl);

legend('Error','|p|/max < 1e-3','Location','best');

% ---------- save with suffix ----------
local_save_fig_png(fig2, sprintf('King_vs_DIM_1D_z%.2fm__norm', z_use));


%% ============================================================
% EXTRA FIGS (PHASE-FIX for FHT only; DIM unchanged)
% 1) 1D compare (4 subplots): ONLY subplot(3) uses FHT phase-extrapolated
% 2) 2D phase compare (2 subplots): FHT uses phase-extrapolated radial line,
%    DIM uses RAW radial line
%
% NOTE: visualization only (m large -> near-axis low-|p| phase unreliable)
% ============================================================

thr_phase = thr;   % normalized amplitude threshold for phase fixing
Lfit      = 200;    % linear-fit points used for inward extrapolation
thr_u     = 5 * thr_phase;  % 倍数越大越平滑，越小越接近原值 3 5 10

%% -------------------- (A) 1D: 4 subplots (ONLY phase of FHT is replaced) --------------------
pK_raw = pK(:);           % King line (Nx1)
pD_raw = pD(:);           % DIM line  (Nx1)

% ONLY fix FHT phase
pK_fix = local_phase_extrapolate_low_amp(pK_raw, rho_ds(:), thr_phase, thr_u); % King only

magK = abs(pK_raw);
magD = abs(pD_raw);

if fig.unwrap
    phK_fix = unwrap(angle(pK_fix));   % fixed King phase
    phD_raw = unwrap(angle(pD_raw));   % raw DIM phase (unchanged)
else
    phK_fix = (angle(pK_fix));   % fixed King phase
    phD_raw = (angle(pD_raw));   % raw DIM phase (unchanged)
end

rel_err_log = log10( abs(pD_raw - pK_raw) ./ (abs(pK_raw) + eps0) );

fig1_fix = figure('Name',sprintf('King vs DIM-%s (z≈1 m) [FHT-PHASE-EXTRAP]', ...
    upper(string(res.calc.dim.method))), 'position',[120 120 1200 1200]);

subplot(4,1,1);
plot(rho_ds, 20*log10(magK / medium.pref / sqrt(2)), 'LineWidth',1.5); hold on;
plot(rho_ds, 20*log10(magD / medium.pref / sqrt(2)), '--', 'LineWidth',1.5);
grid on; xlabel('\rho (m)'); ylabel('SPL (dB)');
title(sprintf('Magnitude @ z = %.2f m', z_use));
legend('King','DIM','Location','best');

subplot(4,1,2);
plot(rho_ds, magK, 'LineWidth',1.5); hold on;
plot(rho_ds, magD, '--', 'LineWidth',1.5);
grid on; xlabel('\rho (m)'); ylabel('|p| (Pa)');
title(sprintf('Magnitude @ z = %.2f m', z_use));
legend('King','DIM','Location','best');

subplot(4,1,3);
plot(rho_ds, phK_fix, 'LineWidth',1.5); hold on;
plot(rho_ds, phD_raw, '--', 'LineWidth',1.5);
grid on; xlabel('\rho (m)'); ylabel('Phase (rad)');
title(sprintf('Phase: ONLY King extrapolated for |p|/max<%.0e (viz only)', thr_phase));
legend('King (phase-fixed)','DIM','Location','best');

subplot(4,1,4);
plot(rho_ds, rel_err_log, 'LineWidth',1.5);
grid on; xlabel('\rho (m)'); ylabel('log_{10} relative error');
title('log_{10}(|p_{DIM}-p_{King}| / |p_{King}|)');

local_save_fig_png(fig1_fix, sprintf('King_vs_DIM_1D_z%.2fm__FHTphaseExtrap', z_use));


%% -------------------- (B) 2D: phase/pi compare (FHT fixed vs DIM raw) --------------------
% Build rho_fig from rho_ds (NO interpolation), then crop to r_boundary
rho_fig = rho_ds(:).';
pK_line_raw = pK_raw(:).';
pD_line_raw = pD_raw(:).';

idx_rb = find(rho_fig <= r_boundary, 1, 'last');
if isempty(idx_rb); idx_rb = numel(rho_fig); end
rho_fig     = rho_fig(1:idx_rb);
pK_line_raw = pK_line_raw(1:idx_rb);
pD_line_raw = pD_line_raw(1:idx_rb);

% ONLY fix FHT radial line
pK_line_fix = local_phase_extrapolate_low_amp(pK_line_raw, rho_fig, thr_phase, thr_u); % King only

% Polar grid -> Cartesian
[TH, R] = meshgrid(theta_fig, rho_fig);
[X, Y]  = pol2cart(TH, R);

% Phase extension for vortex-m
pK_2D_fix = (pK_line_fix(:) * ones(1, numel(theta_fig))) .* exp(1i*m_use*TH); % FHT fixed
pD_2D_raw = (pD_line_raw(:) * ones(1, numel(theta_fig))) .* exp(1i*m_use*TH); % DIM raw

PH_K = angle(pK_2D_fix)/pi;
PH_D = angle(pD_2D_raw)/pi;

figure('Name',sprintf('xOy Phase/pi @ z=%.2f m, m=%d (FHT fixed vs DIM)', z_use, m_use), ...
    'position',[100 100 1400 650]);

ph_lim = [-1 1];

subplot(1,2,1);
pcolor(X, Y, PH_K); shading flat
xlim([-1.2*r_boundary 1.2*r_boundary]);
ylim([-1.2*r_boundary 1.2*r_boundary]);
pbaspect([1 1 1])
colormap('hsv')
clim(ph_lim);
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$\angle p/\pi$';
set(clb,'Fontsize',18);
set(gca,'position',[0.07 0.12 0.38 0.78]);
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');
xlabel('$x$ (m)','Interpreter','latex','Fontsize',18);
ylabel('$y$ (m)','Interpreter','latex','Fontsize',18);
title(sprintf('FHT (phase-fixed, thr=%.0e)', thr_phase), 'Interpreter','latex','Fontsize',20);

subplot(1,2,2);
pcolor(X, Y, PH_D); shading flat
xlim([-1.2*r_boundary 1.2*r_boundary]);
ylim([-1.2*r_boundary 1.2*r_boundary]);
pbaspect([1 1 1])
colormap('hsv')
clim(ph_lim);
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$\angle p/\pi$';
set(clb,'Fontsize',18);
set(gca,'position',[0.55 0.12 0.38 0.78]);
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');
xlabel('$x$ (m)','Interpreter','latex','Fontsize',18);
ylabel('$y$ (m)','Interpreter','latex','Fontsize',18);
title(sprintf('DIM-%s', upper(string(res.calc.dim.method))), 'Interpreter','latex','Fontsize',20);

sgtitle(sprintf('$xOy$ Phase$/\\pi$ (FHT phase extrapolated only) @ $z=%.2f$ m, $m=%d$, $r\\le%.2f$ m', ...
    z_use, m_use, r_boundary), 'Interpreter','latex','Fontsize',20);

local_save_fig_png(gcf, sprintf('xOy_Phase_z%.2fm_m%d__FHTphaseExtrap', z_use, m_use));














%% ==================== local functions ====================
function mem_mb = local_get_mem_mb()
% Returns MATLAB process memory usage in MB if possible; otherwise NaN.
mem_mb = NaN;

% Prefer Windows: use memory() for MATLAB process peak/usage
if ispc
    try
        m = memory();
        mem_mb = double(m.MemUsedMATLAB) / (1024^2);
        return;
    catch
    end
end

% macOS/Linux: Java heap usage (trend only, not full RSS)
try
    rt = java.lang.Runtime.getRuntime();
    used = double(rt.totalMemory() - rt.freeMemory()); % bytes
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

function local_save_fig_png(fig_handle, base_name)
% Save current figure as PNG into SAVE_DIR when SAVE_PNG is true.
global SAVE_PNG SAVE_DIR
if isempty(SAVE_PNG) || ~SAVE_PNG
    return;
end
if isempty(SAVE_DIR) || ~exist(SAVE_DIR,'dir')
    return;
end

% sanitize filename (avoid illegal chars)
fname = local_sanitize_filename(base_name);
fp = fullfile(SAVE_DIR, [fname, '.png']);

% ensure white background (optional)
set(fig_handle, 'Color', 'w');

% export (prefer exportgraphics if available)
try
    exportgraphics(fig_handle, fp, 'Resolution', 300);
catch
    print(fig_handle, fp, '-dpng', '-r300');
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

function local_write_runinfo_txt(save_dir, medium, source, calc, fig)
% Write run configuration into a text file under save_dir.
% Primary section: explicitly listed user parameters (fixed order)
% Secondary section: auto-dump of remaining fields

fp = fullfile(save_dir, 'run_info.txt');
fid = fopen(fp, 'w');
if fid < 0
    warning('Cannot create run_info.txt in %s', save_dir);
    return;
end
cleanupObj = onCleanup(@() fclose(fid));

written = struct();   % record written fields to avoid duplication

fprintf(fid, '===== RUN INFO =====\n');
fprintf(fid, 'Time: %s\n\n', datestr(datetime('now'), 'yyyy-mm-dd HH:MM:SS'));

%% =========================================================
fprintf(fid, '===== PRIMARY PARAMETERS (USER SPECIFIED) =====\n\n');

% -------------------- medium --------------------
fprintf(fid, '%% -------------------- medium --------------------\n');
fprintf(fid, 'medium.c0 = %.15g;\n', medium.c0);           written.medium.c0 = true;
fprintf(fid, 'medium.rho0 = %.15g;\n', medium.rho0);       written.medium.rho0 = true;
fprintf(fid, 'medium.beta = %.15g;\n', medium.beta);       written.medium.beta = true;
fprintf(fid, 'medium.pref = %.15g;\n', medium.pref);       written.medium.pref = true;
fprintf(fid, 'medium.use_absorp = %s;\n', local_bool_str(medium.use_absorp));
written.medium.use_absorp = true;

if isfield(medium,'atten_handle')
    fprintf(fid, 'medium.atten_handle = @(f) AbsorpAttenCoef(f);\n');
    written.medium.atten_handle = true;
end
fprintf(fid, '\n');

% -------------------- source --------------------
fprintf(fid, '%% -------------------- source --------------------\n');
fprintf(fid, 'source.profile = ''%s'';\n', char(string(source.profile)));
written.source.profile = true;

fprintf(fid, 'source.a = %.15g;\n', source.a);             written.source.a = true;
fprintf(fid, 'source.v0 = %.15g;\n', source.v0);           written.source.v0 = true;
fprintf(fid, 'source.v_ratio = %.15g;\n', source.v_ratio); written.source.v_ratio = true;
fprintf(fid, 'source.m = %d;\n', source.m);                written.source.m = true;
fprintf(fid, 'source.F = %.15g;\n', source.F);             written.source.F = true;

fprintf(fid, 'source.f1 = %.15g;\n', source.f1);           written.source.f1 = true;
fprintf(fid, 'source.fa = %.15g;\n', source.fa);           written.source.fa = true;
fprintf(fid, '\n');

% -------------------- calc.fht --------------------
fprintf(fid, '%% -------------------- calc.fht --------------------\n');
cf = calc.fht;
fprintf(fid, 'calc.fht.N_FHT = %.15g;\n', cf.N_FHT);        written.calc.fht.N_FHT = true;
fprintf(fid, 'calc.fht.rho_max = %.15g;\n', cf.rho_max);    written.calc.fht.rho_max = true;
fprintf(fid, 'calc.fht.Nh_scale = %.15g;\n', cf.Nh_scale);  written.calc.fht.Nh_scale = true;
fprintf(fid, 'calc.fht.NH_scale = %.15g;\n', cf.NH_scale);  written.calc.fht.NH_scale = true;
fprintf(fid, 'calc.fht.Nh_v_scale = %.15g;\n', cf.Nh_v_scale);
written.calc.fht.Nh_v_scale = true;
fprintf(fid, 'calc.fht.zu_max = %.15g;\n', cf.zu_max);      written.calc.fht.zu_max = true;
fprintf(fid, 'calc.fht.za_max = %.15g;\n', cf.za_max);      written.calc.fht.za_max = true;
fprintf(fid, '\n');

% -------------------- calc.dim --------------------
fprintf(fid, '%% -------------------- calc.dim --------------------\n');
fprintf(fid, 'calc.dim.dis_coe = %.15g;\n', calc.dim.dis_coe);
written.calc.dim.dis_coe = true;

fprintf(fid, 'calc.dim.src_discretization = ''%s'';\n', ...
    char(string(calc.dim.src_discretization)));
written.calc.dim.src_discretization = true;
fprintf(fid, '\n');

% -------------------- calc.king --------------------
fprintf(fid, '%% -------------------- calc.king --------------------\n');
fprintf(fid, 'calc.king.gspec_method = ''%s'';\n', ...
    char(string(calc.king.gspec_method)));
written.calc.king.gspec_method = true;

fprintf(fid, 'calc.king.band_refine.enable = %s;\n', ...
    local_bool_str(calc.king.band_refine.enable));
written.calc.king.band_refine.enable = true;
fprintf(fid, '\n');

% -------------------- fig --------------------
fprintf(fid, '%% -------------------- fig --------------------\n');
fprintf(fid, 'fig.unwrap = %s;\n', local_bool_str(fig.unwrap));
fprintf(fid, '\n');

fprintf(fid, '\n');

%% =========================================================
fprintf(fid, '===== SECONDARY PARAMETERS (AUTO DUMP) =====\n\n');

local_dump_struct(fid, medium, 'medium', written);
local_dump_struct(fid, source, 'source', written);
local_dump_struct(fid, calc,   'calc',   written);

fprintf(fid, '\n===== END =====\n');
end

function local_dump_struct(fid, S, prefix, written)
fn = fieldnames(S);
for i = 1:numel(fn)
    f = fn{i};

    % skip already written fields
    if isfield(written, prefix) && isfield(written.(prefix), f)
        continue;
    end

    v = S.(f);
    name = sprintf('%s.%s', prefix, f);

    if isnumeric(v) && isscalar(v)
        fprintf(fid, '%s = %.15g;\n', name, v);
    elseif islogical(v) && isscalar(v)
        fprintf(fid, '%s = %s;\n', name, local_bool_str(v));
    elseif ischar(v) || isstring(v)
        fprintf(fid, '%s = ''%s'';\n', name, char(string(v)));
    elseif isstruct(v)
        local_dump_struct(fid, v, name, struct());
    elseif isa(v,'function_handle')
        fprintf(fid, '%s = %s;\n', name, func2str(v));
    end
end
end

function s = local_bool_str(x)
if logical(x)
    s = 'true';
else
    s = 'false';
end
end

function seg = local_mask_to_segments(x, mask)
% seg: Kx2 array, each row [xL xR] for contiguous true regions
mask = mask(:);
x = x(:);
if numel(mask) ~= numel(x)
    error('local_mask_to_segments:DimMismatch', 'mask and x must have same length.');
end

d = diff([false; mask; false]);
i1 = find(d==1);
i2 = find(d==-1)-1;

seg = zeros(numel(i1), 2);
for k = 1:numel(i1)
    seg(k,:) = [x(i1(k)), x(i2(k))];
end
end

function h = local_mark_segments_strong(seg, yl)
% Dark, visible shading for low-amplitude regions
% Returns patch handles for legend use

h = gobjects(0);
if isempty(seg), return; end

for k = 1:size(seg,1)
    xL = seg(k,1);
    xR = seg(k,2);

    h(end+1) = patch( ...
        [xL xR xR xL], ...
        [yl(1) yl(1) yl(2) yl(2)], ...
        [0.65 0.65 0.65], ...   % darker gray
        'EdgeColor','none', ...
        'FaceAlpha',0.45);      %#ok<AGROW>
end
end

function p_out = local_phase_extrapolate_low_amp(p_in, rho, thr, thr_u)
% Scheme A (recommended):
% - Build extrapolated phase using PCHIP on reliable region (|p|/max >= thr)
% - Replace low-amplitude phase by extrapolated phase
% - Use smooth blending in [thr, thr_u] to avoid a kink at the boundary
% - Magnitude is unchanged (visualization only)
%
% Inputs:
%   p_in  : complex vector (Nx1 or 1xN)
%   rho   : real vector (same length as p_in)
%   thr   : normalized threshold (e.g., 1e-3)
%   thr_u : upper threshold for blending (e.g., 5e-3). If omitted, thr_u=5*thr
%
% Output:
%   p_out : same size as p_in, magnitude unchanged, phase modified near axis

if nargin < 4 || isempty(thr_u)
    thr_u = 5*thr;
end
thr_u = max(thr_u, thr*(1+1e-12));  % ensure thr_u > thr

p = p_in(:);
r = rho(:);

if numel(p) ~= numel(r)
    error('local_phase_extrapolate_low_amp:DimMismatch', 'p_in and rho must have same length.');
end

A = abs(p);
Amax = max(A);
if ~(isfinite(Amax)) || Amax <= 0
    p_out = p_in;
    return;
end
an = A / Amax;

% unwrap phase from raw complex data
phi_raw = unwrap(angle(p));

% reliable indices
I = find(an >= thr);
if numel(I) < 10
    % not enough support for a stable pchip
    p_out = p_in;
    return;
end

% ensure r is strictly increasing for pchip (it should be in your FHT grid)
% if not, sort (robust)
if any(diff(r) <= 0)
    [rS, idxS] = sort(r);
    pS = p(idxS);
    AS = abs(pS);
    anS = AS / max(AS);
    phiS = unwrap(angle(pS));
    IS = find(anS >= thr);
    if numel(IS) < 10
        p_out = p_in; return;
    end
    phi_ext_S = pchip(rS(IS), phiS(IS), rS);
    % blending weights on sorted grid
    wS = (anS - thr) / (thr_u - thr);
    wS = min(max(wS, 0), 1);
    phi_new_S = (1-wS).*phi_ext_S + wS.*phiS;
    p2S = AS .* exp(1j*phi_new_S);
    % unsort back
    p2 = zeros(size(p2S));
    p2(idxS) = p2S;
else
    % extrapolated phase from reliable points (shape-preserving)
    phi_ext = pchip(r(I), phi_raw(I), r);

    % blending weights
    w = (an - thr) / (thr_u - thr);
    w = min(max(w, 0), 1);

    % phase mix: low amp -> extrapolated, high amp -> raw
    phi_new = (1-w).*phi_ext + w.*phi_raw;

    % rebuild with same magnitude
    p2 = A .* exp(1j*phi_new);
end

% restore original shape
p_out = reshape(p2, size(p_in));
end
