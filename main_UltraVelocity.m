%% ============================================================
% MAIN: King/FHT velocity field (v_rho, v_phi, v_z)
% 1) 1D @ z≈1 m: three figures (each: magnitude + phase)
% 2) 2D @ z≈1 m plane (via phase extension exp(j*m*theta)):
%    - show velocity vector v (quiver of Re{v_x,v_y} with |v| background)
%    - show inner product s = v·v = v_x^2 + v_y^2 + v_z^2 (magnitude + phase)
%
% Dependencies:
%   - make_source_velocity.m (you provided)
%   - calc_ultrasound_velocity_field.m (King/FHT velocity solver; see earlier)
%   - m_FHT.m, solve_kappa0.m
%   - MyColor.m, fontsize.m (your plotting helpers)
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
source.m = 10;
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

% --- DIM source discretization (not used here, but keep consistent) ---
calc.dim.use_freq = 'f2';
calc.dim.dis_coe = 16 * 1;
calc.dim.margin = 1;
calc.dim.src_discretization = 'polar';

% --- King analytic spectrum stability ---
calc.king.gspec_method = 'analytic';
calc.king.eps_kzz = 1e-3;
calc.king.eps_phase = calc.king.eps_kzz;
calc.king.kz_min = 1e-12;

% --- band refine (not used in velocity solver below) ---
calc.king.band_refine.enable = false;

%% -------------------- fig setup --------------------
fig.unwrap = false;

%% -------------------- save setup --------------------
if SAVE_PNG
    tstr = datestr(datetime('now'), 'mmdd_HHMM');
    a_str = sprintf('%.2fm', source.a);
    m_str = sprintf('m=%d', source.m);
    SAVE_DIR = sprintf('%s_%s__%s', a_str, m_str, tstr);
    if ~exist(SAVE_DIR, 'dir'); mkdir(SAVE_DIR); end
else
    SAVE_DIR = '';
end

%% -------------------- target plane --------------------
z_target = 1;
ds = 16 * 2;

%% ============================================================
% PREPARE GRID (from make_source_velocity)
%% ============================================================
[~, fht_tmp] = make_source_velocity(source, medium, calc); %#ok<ASGLU>
rho_full = fht_tmp.xh(:);
z_full   = fht_tmp.z_ultra(:);

[~, iz] = min(abs(z_full - z_target));
z_use = z_full(iz);

rho_ds = rho_full(1:ds:end);

fprintf('z_target=%.3f, z_use=%.6f, ds=%d\n', z_target, z_use, ds);

%% ============================================================
% Compute velocity field (King/FHT)
%% ============================================================
fprintf('\n==================== KING VELOCITY (FHT) ====================\n');

t0 = tic;
resV = calc_ultrasound_velocity_field(source, medium, calc, 'king');
tV = toc(t0);

fprintf('Velocity compute time: %.3f s\n', tV);

% ---- extract f1 velocities at z_use and downsample in rho ----
zK   = resV.king.z(:);
rhoK = resV.king.rho(:);

[~, izK] = min(abs(zK - z_use));

vR_full = resV.king.v_rho_f1(:, izK);   % Nr x 1
vP_full = resV.king.v_phi_f1(:, izK);
vZ_full = resV.king.v_z_f1(:,   izK);

vR = vR_full(1:ds:end);
vP = vP_full(1:ds:end);
vZ = vZ_full(1:ds:end);

% ---- phases ----
if fig.unwrap
    phR = unwrap(angle(vR));
    phP = unwrap(angle(vP));
    phZ = unwrap(angle(vZ));
else
    phR = angle(vR);
    phP = angle(vP);
    phZ = angle(vZ);
end

%% ============================================================
% 1D FIGS: rho-direction @ z≈1 m
% Three figures: (v_rho, v_phi, v_z), each has magnitude+phase
%% ============================================================

local_plot_1D_mag_phase(rho_ds, vR, phR, ...
    sprintf('$|v_\\rho|$ @ $z=%.3f$ m', z_use), ...
    sprintf('$\\angle v_\\rho$ @ $z=%.3f$ m', z_use), ...
    sprintf('vRho_1D_z%.3fm', z_use));

local_plot_1D_mag_phase(rho_ds, vP, phP, ...
    sprintf('$|v_\\varphi|$ @ $z=%.3f$ m', z_use), ...
    sprintf('$\\angle v_\\varphi$ @ $z=%.3f$ m', z_use), ...
    sprintf('vPhi_1D_z%.3fm', z_use));

local_plot_1D_mag_phase(rho_ds, vZ, phZ, ...
    sprintf('$|v_z|$ @ $z=%.3f$ m', z_use), ...
    sprintf('$\\angle v_z$ @ $z=%.3f$ m', z_use), ...
    sprintf('vZ_1D_z%.3fm', z_use));

%% ============================================================
% 2D plane @ z≈1 m (phase extension)
% v(r,theta) = v_line(r) * exp(j*m*theta)
% Convert to (v_x, v_y, v_z), show:
%   - velocity vector (quiver of Re{v_x,v_y} over |v| background)
%   - s = v·v = v_x^2 + v_y^2 + v_z^2 : magnitude and phase
%% ============================================================

% ---- view settings ----
r_boundary = 0.30;           % radius to show (m)
theta_fig  = 0:0.01:2*pi;    % angular grid

% ---- m ----
m_use = resV.king.m;

% ---- crop radial samples (NO interp) ----
rho_fig = rho_ds(:).';
idx_rb = find(rho_fig <= r_boundary, 1, 'last');
if isempty(idx_rb); idx_rb = numel(rho_fig); end

rho_fig = rho_fig(1:idx_rb);
vR_line = vR(:).'; vR_line = vR_line(1:idx_rb);
vP_line = vP(:).'; vP_line = vP_line(1:idx_rb);
vZ_line = vZ(:).'; vZ_line = vZ_line(1:idx_rb);

% ---- polar grid -> Cartesian ----
[TH, R] = meshgrid(theta_fig, rho_fig);
[X, Y]  = pol2cart(TH, R);

% ---- phase extension (complex) ----
Ephi = exp(1i * m_use * TH);

vR_2D = (vR_line(:) * ones(1, numel(theta_fig))) .* Ephi;
vP_2D = (vP_line(:) * ones(1, numel(theta_fig))) .* Ephi;
vZ_2D = (vZ_line(:) * ones(1, numel(theta_fig))) .* Ephi;

% ---- cylindrical -> Cartesian for in-plane components ----
% v_x = v_rho*cosθ - v_phi*sinθ
% v_y = v_rho*sinθ + v_phi*cosθ
cT = cos(TH); sT = sin(TH);
vX_2D = vR_2D .* cT - vP_2D .* sT;
vY_2D = vR_2D .* sT + vP_2D .* cT;

% ---- vector magnitude (use |v| = sqrt(|vx|^2+|vy|^2+|vz|^2)) ----
Vmag = sqrt(abs(vX_2D).^2 + abs(vY_2D).^2 + abs(vZ_2D).^2);

% ---- inner product with itself (NO conjugate): s = v·v ----
S = vX_2D.^2 + vY_2D.^2 + vZ_2D.^2;

Smag = abs(S);
Sph  = angle(S)/pi;   % phase/pi

%% -------------------- FIG: 2D velocity vector + |v| background --------------------
figure('Name',sprintf('2D |v| + quiver @ z=%.3f m, m=%d', z_use, m_use), ...
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

% overlay quiver of Re{v_x, v_y}
hold on;
qstep_r = max(1, round(numel(rho_fig)/25));
qstep_t = max(1, round(numel(theta_fig)/60));
Xq = X(1:qstep_r:end, 1:qstep_t:end);
Yq = Y(1:qstep_r:end, 1:qstep_t:end);

vXq = real(vX_2D(1:qstep_r:end, 1:qstep_t:end));
vYq = real(vY_2D(1:qstep_r:end, 1:qstep_t:end));

% scale arrows for visibility (relative)
scale = max(Vmag(:));
if scale > 0
    vXq = vXq / scale;
    vYq = vYq / scale;
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

sgtitle(sprintf('$xOy$ velocity @ $z=%.3f$ m, $m=%d$, $r\\le%.2f$ m', ...
    z_use, m_use, r_boundary), 'Interpreter','latex','Fontsize',20);

local_save_fig_png(gcf, sprintf('xOy_Vmag_quiver_z%.3fm_m%d', z_use, m_use));

%% -------------------- FIG: inner product S = v·v (magnitude + phase) --------------------
figure('Name',sprintf('2D S=v·v @ z=%.3f m, m=%d', z_use, m_use), ...
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
colormap('hsv');
clim([-1 1]);
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

sgtitle(sprintf('$xOy$ scalar $s=v\\cdot v$ @ $z=%.3f$ m, $m=%d$, $r\\le%.2f$ m', ...
    z_use, m_use, r_boundary), 'Interpreter','latex','Fontsize',20);

local_save_fig_png(gcf, sprintf('xOy_vdotv_z%.3fm_m%d', z_use, m_use));

fprintf('\nDONE.\n');

%% ==================== local functions ====================

function local_plot_1D_mag_phase(rho, v, ph, ttl_mag, ttl_ph, save_name)
global SAVE_PNG SAVE_DIR

mag = abs(v);

figure('Name',save_name,'position',[100 100 1200 650]);

subplot(2,1,1);
plot(rho, mag, 'LineWidth', 1.5);
grid on;
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');
xlabel('$\rho$ (m)','Interpreter','latex','Fontsize',18);
ylabel('Magnitude (m/s)','Interpreter','latex','Fontsize',18);
title(ttl_mag,'Interpreter','latex','Fontsize',20);

subplot(2,1,2);
plot(rho, ph, 'LineWidth', 1.5);
grid on;
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');
xlabel('$\rho$ (m)','Interpreter','latex','Fontsize',18);
ylabel('Phase (rad)','Interpreter','latex','Fontsize',18);
title(ttl_ph,'Interpreter','latex','Fontsize',20);

local_save_fig_png(gcf, save_name);
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
