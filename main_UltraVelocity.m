%% ============================================================
% MAIN: King vs DIM velocity @ z≈1 m
% - 1D compare on phi=0 line (y=0, x=rho): v_rho, v_phi, v_z
%   Each component: 3 subplots (magnitude / phase / log relative error)
%
% - 2D xOy @ z≈1 m (King/FHT only, via phase extension exp(j*m*theta)):
%   (1) |v| + quiver(Re{v_x,v_y})
%   (2) phase of v_z ( /pi )
%   (3) s = v·v (NO conjugate): magnitude + phase(/pi)
%
% Dependencies:
%   - make_source_velocity.m
%   - calc_ultrasound_velocity_field.m   (your latest version incl. DIM block)
%   - m_FHT.m, solve_kappa0.m
%   - MyColor.m, fontsize.m
% ============================================================

clear; clc; close all;

%% -------------------- save figures (GLOBAL control) --------------------
global SAVE_PNG SAVE_DIR
SAVE_PNG = false;

%% -------------------- medium --------------------
medium.c0 = 343;
medium.rho0 = 1.21;
medium.beta = 1.2;
medium.pref = 2e-5;
medium.use_absorp = true;
medium.atten_handle = @(f) AbsorpAttenCoef(f);

%% -------------------- source --------------------
source.profile = 'Vortex-m';
source.a = 0.1;
source.v0 = 0.172;
source.v_ratio = 1;
source.m = 5;
source.F = 0.2;

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

% --- DIM source discretization (must be consistent with make_source_velocity) ---
calc.dim.use_freq = 'f2';
calc.dim.dis_coe = 16 * 1;
calc.dim.margin = 1;
calc.dim.src_discretization = 'polar'; % 'cart' or 'polar'

% --- King analytic spectrum stability ---
calc.king.gspec_method = 'analytic';
calc.king.eps_kzz = 1e-3;
calc.king.eps_phase = calc.king.eps_kzz;
calc.king.kz_min = 1e-12;

% --- band refine (optional) ---
calc.king.band_refine.enable = true;

%% -------------------- fig setup --------------------
fig.unwrap = false;

%% -------------------- save setup --------------------
if SAVE_PNG
    tstr = datestr(datetime('now'), 'mmdd_HHMM');
    a_str = sprintf('%.2fm', source.a);
    m_str = sprintf('m=%d', source.m);
    SAVE_DIR = sprintf('%s_%s__%s', a_str, m_str, tstr);
    if ~exist(SAVE_DIR, 'dir'), mkdir(SAVE_DIR); end
else
    SAVE_DIR = '';
end

%% -------------------- target plane & ds --------------------
z_target = 1.0;

% ds: downsample factor on rho grid (default 32)
ds = 32;

%% ==================== PROFILER: helpers ====================
get_mem_mb = @() local_get_mem_mb();
fmt_mem    = @(x) local_fmt_mem(x);

%% ============================================================
% PREPARE GRID (from make_source_velocity)
% ============================================================
[~, fht_tmp] = make_source_velocity(source, medium, calc); %#ok<ASGLU>
rho_full = fht_tmp.xh(:);
z_full   = fht_tmp.z_ultra(:);

[~, iz] = min(abs(z_full - z_target));
z_use = z_full(iz);

rho_ds = rho_full(1:ds:end);

fprintf('z_target=%.3f, z_use=%.6f, ds=%d\n', z_target, z_use, ds);

%% ============================================================
% Prepare DIM observation grid (phi=0 line): x=rho_ds, y=0, z=z_use
% ============================================================
if ~isfield(calc,'dim') || ~isstruct(calc.dim), calc.dim = struct(); end

calc.dim.z_use = z_use;  % required by DIM block in velocity solver

calc.dim.obs_grid = struct();
calc.dim.obs_grid.x = rho_ds(:).';  % 1 x Nx
calc.dim.obs_grid.y = 0;            % 1 x Ny (=1)
calc.dim.obs_grid.z = z_use;        % scalar plane

% blocks for DIM velocity
calc.dim.block_size     = 20000;
calc.dim.src_block_size = 5000;

%% ============================================================
% Compute velocity field (King + DIM)
%% ============================================================
fprintf('\n==================== VELOCITY (KING + DIM) ====================\n');
mem0 = get_mem_mb(); t0 = tic;

resV = calc_ultrasound_velocity_field(source, medium, calc, 'both');

tAll = toc(t0); mem1 = get_mem_mb();
fprintf('Velocity time (both): %.3f s\n', tAll);
fprintf('Velocity memory: start %s, end %s, delta %s\n', fmt_mem(mem0), fmt_mem(mem1), fmt_mem(mem1-mem0));

%% ============================================================
% Extract King (f1) on same z_use, downsample to rho_ds
%% ============================================================
zK   = resV.king.z(:);
rhoK = resV.king.rho(:);
[~, izK] = min(abs(zK - z_use));

vR_king_full = resV.king.v_rho_f1(:, izK);
vP_king_full = resV.king.v_phi_f1(:, izK);
vZ_king_full = resV.king.v_z_f1(:,   izK);

vR_king = vR_king_full(1:ds:end);
vP_king = vP_king_full(1:ds:end);
vZ_king = vZ_king_full(1:ds:end);

%% ============================================================
% Extract DIM (f1) on line (y=0): Ny=1, Nx=numel(rho_ds)
%% ============================================================
% NOTE: your DIM block packs:
% result.dim.v_rho_f1, v_phi_f1, v_z_f1 as Ny x Nx arrays
vR_dim = squeeze(resV.dim.v_rho_f1(1,:)).';
vP_dim = squeeze(resV.dim.v_phi_f1(1,:)).';
vZ_dim = squeeze(resV.dim.v_z_f1(1,:)).';

%% ============================================================
% Phase & errors
%% ============================================================
if fig.unwrap
    phR_k = unwrap(angle(vR_king)); phR_d = unwrap(angle(vR_dim));
    phP_k = unwrap(angle(vP_king)); phP_d = unwrap(angle(vP_dim));
    phZ_k = unwrap(angle(vZ_king)); phZ_d = unwrap(angle(vZ_dim));
else
    phR_k = angle(vR_king); phR_d = angle(vR_dim);
    phP_k = angle(vP_king); phP_d = angle(vP_dim);
    phZ_k = angle(vZ_king); phZ_d = angle(vZ_dim);
end

eps0 = 1e-12;

errR_log = log10( abs(vR_dim - vR_king) ./ (abs(vR_king) + eps0) );
errP_log = log10( abs(vP_dim - vP_king) ./ (abs(vP_king) + eps0) );
errZ_log = log10( abs(vZ_dim - vZ_king) ./ (abs(vZ_king) + eps0) );

%% ============================================================
% 1D FIGS: compare v_rho / v_phi / v_z  (3 subplots each)
%% ============================================================

% -------------------- v_rho --------------------
figure('Name',sprintf('v_rho: King vs DIM @ z=%.3f m', z_use), 'position',[100 100 1200 900]);

subplot(3,1,1);
plot(rho_ds, abs(vR_king), 'LineWidth',1.5); hold on;
plot(rho_ds, abs(vR_dim),  '--', 'LineWidth',1.5);
grid on;
xlabel('\rho (m)'); ylabel('|v_\rho| (m/s)');
title(sprintf('Magnitude: $v_\\rho$ @ $z=%.3f$ m', z_use), 'Interpreter','latex');
legend('King','DIM','Location','best');

subplot(3,1,2);
plot(rho_ds, phR_k, 'LineWidth',1.5); hold on;
plot(rho_ds, phR_d, '--', 'LineWidth',1.5);
grid on;
xlabel('\rho (m)'); ylabel('Phase (rad)');
title('Phase', 'Interpreter','latex');
legend('King','DIM','Location','best');

subplot(3,1,3);
plot(rho_ds, errR_log, 'LineWidth',1.5);
grid on;
xlabel('\rho (m)');
ylabel('log_{10} relative error');
title('log_{10}(|v_{DIM}-v_{King}| / (|v_{King}|+\epsilon))', 'Interpreter','latex');

local_save_fig_png(gcf, sprintf('vRho_1D_compare3_z%.3fm', z_use));

% -------------------- v_phi --------------------
figure('Name',sprintf('v_phi: King vs DIM @ z=%.3f m', z_use), 'position',[120 120 1200 900]);

subplot(3,1,1);
plot(rho_ds, abs(vP_king), 'LineWidth',1.5); hold on;
plot(rho_ds, abs(vP_dim),  '--', 'LineWidth',1.5);
grid on;
xlabel('\rho (m)'); ylabel('|v_\phi| (m/s)');
title(sprintf('Magnitude: $v_\\phi$ @ $z=%.3f$ m', z_use), 'Interpreter','latex');
legend('King','DIM','Location','best');

subplot(3,1,2);
plot(rho_ds, phP_k, 'LineWidth',1.5); hold on;
plot(rho_ds, phP_d, '--', 'LineWidth',1.5);
grid on;
xlabel('\rho (m)'); ylabel('Phase (rad)');
title('Phase', 'Interpreter','latex');
legend('King','DIM','Location','best');

subplot(3,1,3);
plot(rho_ds, errP_log, 'LineWidth',1.5);
grid on;
xlabel('\rho (m)');
ylabel('log_{10} relative error');
title('log_{10}(|v_{DIM}-v_{King}| / (|v_{King}|+\epsilon))', 'Interpreter','latex');

local_save_fig_png(gcf, sprintf('vPhi_1D_compare3_z%.3fm', z_use));

% -------------------- v_z --------------------
figure('Name',sprintf('v_z: King vs DIM @ z=%.3f m', z_use), 'position',[140 140 1200 900]);

subplot(3,1,1);
plot(rho_ds, abs(vZ_king), 'LineWidth',1.5); hold on;
plot(rho_ds, abs(vZ_dim),  '--', 'LineWidth',1.5);
grid on;
xlabel('\rho (m)'); ylabel('|v_z| (m/s)');
title(sprintf('Magnitude: $v_z$ @ $z=%.3f$ m', z_use), 'Interpreter','latex');
legend('King','DIM','Location','best');

subplot(3,1,2);
plot(rho_ds, phZ_k, 'LineWidth',1.5); hold on;
plot(rho_ds, phZ_d, '--', 'LineWidth',1.5);
grid on;
xlabel('\rho (m)'); ylabel('Phase (rad)');
title('Phase', 'Interpreter','latex');
legend('King','DIM','Location','best');

subplot(3,1,3);
plot(rho_ds, errZ_log, 'LineWidth',1.5);
grid on;
xlabel('\rho (m)');
ylabel('log_{10} relative error');
title('log_{10}(|v_{DIM}-v_{King}| / (|v_{King}|+\epsilon))', 'Interpreter','latex');

local_save_fig_png(gcf, sprintf('vZ_1D_compare3_z%.3fm', z_use));

%% ============================================================
% 2D xOy @ z≈1 m (King only): |v|, phase(vz), and s=v·v
% via phase extension exp(j*m*theta)
%% ============================================================

% ---- view settings ----
r_boundary = 0.30;
theta_fig  = 0:0.01:2*pi;

% ---- m ----
m_use = resV.king.m;

% ---- crop radial samples (NO interp) ----
rho_fig = rho_ds(:).';
idx_rb  = find(rho_fig <= r_boundary, 1, 'last');
if isempty(idx_rb), idx_rb = numel(rho_fig); end

rho_fig = rho_fig(1:idx_rb);
vR_line = vR_king(:).'; vR_line = vR_line(1:idx_rb);
vP_line = vP_king(:).'; vP_line = vP_line(1:idx_rb);
vZ_line = vZ_king(:).'; vZ_line = vZ_line(1:idx_rb);

% ---- polar grid -> Cartesian ----
[TH, R] = meshgrid(theta_fig, rho_fig);
[X, Y]  = pol2cart(TH, R);

% ---- phase extension (complex) ----
Ephi = exp(1i * m_use * TH);

vR_2D = (vR_line(:) * ones(1, numel(theta_fig))) .* Ephi;
vP_2D = (vP_line(:) * ones(1, numel(theta_fig))) .* Ephi;
vZ_2D = (vZ_line(:) * ones(1, numel(theta_fig))) .* Ephi;

% ---- cylindrical -> Cartesian ----
cT = cos(TH); sT = sin(TH);
vX_2D = vR_2D .* cT - vP_2D .* sT;
vY_2D = vR_2D .* sT + vP_2D .* cT;

% ---- |v| ----
Vmag = sqrt(abs(vX_2D).^2 + abs(vY_2D).^2 + abs(vZ_2D).^2);

% ---- phase(vz)/pi ----
VZph = angle(vZ_2D) / pi;

% ---- s = v·v (NO conjugate) ----
S   = vX_2D.^2 + vY_2D.^2 + vZ_2D.^2;
Smag = abs(S);
Sph  = angle(S)/pi;

%% -------------------- FIG: |v| + quiver(Re{vx,vy}) and Re{vz} --------------------
figure('Name',sprintf('King: xOy |v| + quiver @ z=%.3f m, m=%d', z_use, m_use), ...
    'position',[100 100 1400 600]);

subplot(1,2,1);
pcolor(X, Y, Vmag); shading flat
colormap(MyColor('vik'));
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$|v|$ (m/s)';
set(clb,'Fontsize',18);

axis equal
xlim([-1.2*r_boundary 1.2*r_boundary]);
ylim([-1.2*r_boundary 1.2*r_boundary]);
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');
xlabel('$x$ (m)','Interpreter','latex','Fontsize',18);
ylabel('$y$ (m)','Interpreter','latex','Fontsize',18);
title(sprintf('$|v|$ @ $z=%.3f$ m', z_use), 'Interpreter','latex','Fontsize',20);

hold on;
qstep_r = max(1, round(numel(rho_fig)/25));
qstep_t = max(1, round(numel(theta_fig)/60));
Xq  = X(1:qstep_r:end, 1:qstep_t:end);
Yq  = Y(1:qstep_r:end, 1:qstep_t:end);
vXq = real(vX_2D(1:qstep_r:end, 1:qstep_t:end));
vYq = real(vY_2D(1:qstep_r:end, 1:qstep_t:end));

sc = max(Vmag(:));
if sc > 0
    vXq = vXq / sc;
    vYq = vYq / sc;
end
quiver(Xq, Yq, vXq, vYq, 0.6, 'k', 'LineWidth', 1.0);
hold off;

subplot(1,2,2);
pcolor(X, Y, real(vZ_2D)); shading flat
colormap(MyColor('vik'));
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$\Re\{v_z\}$ (m/s)';
set(clb,'Fontsize',18);

axis equal
xlim([-1.2*r_boundary 1.2*r_boundary]);
ylim([-1.2*r_boundary 1.2*r_boundary]);
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');
xlabel('$x$ (m)','Interpreter','latex','Fontsize',18);
ylabel('$y$ (m)','Interpreter','latex','Fontsize',18);
title(sprintf('$\\Re\\{v_z\\}$ @ $z=%.3f$ m', z_use), 'Interpreter','latex','Fontsize',20);

sgtitle(sprintf('King velocity @ $z=%.3f$ m, $m=%d$, $r\\le%.2f$ m', z_use, m_use, r_boundary), ...
    'Interpreter','latex','Fontsize',20);

local_save_fig_png(gcf, sprintf('King_xOy_Vmag_quiver_z%.3fm_m%d', z_use, m_use));

%% -------------------- FIG: phase(vz)/pi --------------------
figure('Name',sprintf('King: xOy phase(vz)/pi @ z=%.3f m, m=%d', z_use, m_use), ...
    'position',[120 120 700 650]);

pcolor(X, Y, VZph); shading flat
colormap('hsv'); clim([-1 1]);
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$\angle v_z/\pi$';
set(clb,'Fontsize',18);

axis equal
xlim([-1.2*r_boundary 1.2*r_boundary]);
ylim([-1.2*r_boundary 1.2*r_boundary]);
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');
xlabel('$x$ (m)','Interpreter','latex','Fontsize',18);
ylabel('$y$ (m)','Interpreter','latex','Fontsize',18);
title(sprintf('$\\angle v_z/\\pi$ @ $z=%.3f$ m', z_use), 'Interpreter','latex','Fontsize',20);

local_save_fig_png(gcf, sprintf('King_xOy_PhaseVz_z%.3fm_m%d', z_use, m_use));

%% -------------------- FIG: inner product S=v·v (magnitude + phase) --------------------
figure('Name',sprintf('King: xOy S=v·v @ z=%.3f m, m=%d', z_use, m_use), ...
    'position',[100 100 1400 600]);

subplot(1,2,1);
pcolor(X, Y, Smag); shading flat
colormap(MyColor('vik'));
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$|s|$, $s=v\cdot v$';
set(clb,'Fontsize',18);

axis equal
xlim([-1.2*r_boundary 1.2*r_boundary]);
ylim([-1.2*r_boundary 1.2*r_boundary]);
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');
xlabel('$x$ (m)','Interpreter','latex','Fontsize',18);
ylabel('$y$ (m)','Interpreter','latex','Fontsize',18);
title('$|v\cdot v|$', 'Interpreter','latex','Fontsize',20);

subplot(1,2,2);
pcolor(X, Y, Sph); shading flat
colormap('hsv'); clim([-1 1]);
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$\angle (v\cdot v)/\pi$';
set(clb,'Fontsize',18);

axis equal
xlim([-1.2*r_boundary 1.2*r_boundary]);
ylim([-1.2*r_boundary 1.2*r_boundary]);
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');
xlabel('$x$ (m)','Interpreter','latex','Fontsize',18);
ylabel('$y$ (m)','Interpreter','latex','Fontsize',18);
title('$\angle (v\cdot v)/\pi$', 'Interpreter','latex','Fontsize',20);

sgtitle(sprintf('King scalar $s=v\\cdot v$ @ $z=%.3f$ m, $m=%d$, $r\\le%.2f$ m', ...
    z_use, m_use, r_boundary), 'Interpreter','latex','Fontsize',20);

local_save_fig_png(gcf, sprintf('King_xOy_vdotv_z%.3fm_m%d', z_use, m_use));

fprintf('\nDONE.\n');

%% ==================== local functions ====================

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
if isnan(x), s = 'NaN';
else,        s = sprintf('%.1f MB', x);
end
end

function local_save_fig_png(fig_handle, base_name)
global SAVE_PNG SAVE_DIR
if isempty(SAVE_PNG) || ~SAVE_PNG, return; end
if isempty(SAVE_DIR) || ~exist(SAVE_DIR,'dir'), return; end

fname = local_sanitize_filename(base_name);
fp = fullfile(SAVE_DIR, [fname, '.png']);
set(fig_handle, 'Color', 'w');

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
