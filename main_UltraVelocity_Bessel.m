%% ============================================================
% MAIN: DIM-only Bessel beam velocity @ z≈1 m (2D plane)
% - Calls:
%   make_source_velocity.m
%   calc_ultrasound_velocity_field.m   (method='dim')
%
% - Source: multi-ring superposition + optional axicon radial phase
%   v_rho(r) = w_ring(bin(r)) * exp(-j*k*sin(theta)*r)
%   v_n(r,phi) = v_rho(r) * exp(j*n*phi)     (n via source.m_custom)
%
% - Plots (DIM, f2) on z=z_use plane:
%   0) |v| (interp) + quiver(Re{v_x,v_y})      + phase(vz)/pi (flat)
%   1) v_x : |v_x| (interp) + angle(v_x)/pi (flat)
%   2) v_y : |v_y| (interp) + angle(v_y)/pi (flat)
%   3) v_z : |v_z| (interp) + angle(v_z)/pi (flat)
%   4) s=v·v (NO conjugate): |s| (interp) + angle(s)/pi (flat)
%
% Style:
%   - magnitude: colormap(MyColor('vik')), shading interp
%   - phase:     colormap(hsv),           shading flat, range [-1,1]
% ============================================================

clear; clc; close all;

%% -------------------- save figures (GLOBAL control) --------------------
global SAVE_PNG SAVE_DIR
SAVE_PNG = false;
SAVE_DIR = '';

%% -------------------- profiler helpers --------------------
get_mem_mb = @() local_get_mem_mb();
fmt_mem    = @(x) local_fmt_mem(x);

%% -------------------- medium --------------------
medium.c0   = 343;
medium.rho0 = 1.21;
medium.beta = 1.2;
medium.pref = 2e-5;
medium.use_absorp = false;
medium.atten_handle = [];

%% -------------------- source (multi-ring superposition) --------------------
n = 2;                      % n阶贝塞尔（同时也是角向阶数）
source.profile = 'Custom';  % 用 Custom 自定义 v_rho(r)
source.a = 0.1;
source.v0 = 0.172;
source.v_ratio = 1;

% 关键：需要你在 make_source_velocity.m 的 Custom 分支里启用 m_custom
% （见文末“必要改动”说明）
source.m_custom = n;        % vn(x,y)=v_rho(r)*exp(j*n*phi)

source.f1 = 42e3;
source.fa = 2e3;
source.f2 = source.f1 + source.fa;

% --------- 多环参数：边界 r_edges + 复权重 w_ring ----------
Q = 24;
r_edges = linspace(0, source.a, Q+1);

w_ring = zeros(Q,1);
q0 = round(0.7*Q);
w_ring(max(1,q0-1):min(Q,q0+1)) = 1;
w_ring = w_ring / max(abs(w_ring)+eps);

% --------- axicon 径向相位（可选） ----------
theta_deg = 10;
theta = deg2rad(theta_deg);

k1r = 2*pi*source.f1/medium.c0;
k2r = 2*pi*source.f2/medium.c0;
kr1 = k1r*sin(theta);
kr2 = k2r*sin(theta);

source.custom_vrho_handle_1 = @(rho) source.v0 .* local_vrho_rings(rho, r_edges, w_ring) .* exp(-1i*kr1.*rho);
source.custom_vrho_handle_2 = @(rho) (source.v_ratio*source.v0) .* local_vrho_rings(rho, r_edges, w_ring) .* exp(-1i*kr2.*rho);

%% -------------------- calc --------------------
calc = struct();

% FHT fields not used, but required by make_source_velocity() defaults
calc.fht.N_FHT = 32768;
calc.fht.rho_max = 0.25;
calc.fht.Nh_scale = 1.2;
calc.fht.NH_scale = 4;
calc.fht.Nh_v_scale = 1.1;
calc.fht.zu_max = 1.2;
calc.fht.za_max = 1;

% DIM source discretization
calc.dim.use_freq = 'f2';
calc.dim.dis_coe  = 16;
calc.dim.margin   = 1;
calc.dim.src_discretization = 'polar';   % 'cart' or 'polar'

% DIM blocks
calc.dim.block_size     = 20000;
calc.dim.src_block_size = 5000;

%% -------------------- target plane (z≈1 m) --------------------
z_target = 1.0;

% Get consistent z grid from make_source_velocity
[~, fht_tmp] = make_source_velocity(source, medium, calc); %#ok<ASGLU>
z_full = fht_tmp.z_ultra(:);
[~, iz] = min(abs(z_full - z_target));
z_use = z_full(iz);

calc.dim.z_use = z_use;

%% -------------------- observation grid (xOy) --------------------
r_boundary = 0.008;        % plot window
dx_obs = 0.0002;           % observation spacing

x_obs = -r_boundary:dx_obs:r_boundary;
y_obs = -r_boundary:dx_obs:r_boundary;

calc.dim.obs_grid = struct();
calc.dim.obs_grid.x = x_obs;
calc.dim.obs_grid.y = y_obs;
calc.dim.obs_grid.z = z_use;

%% ============================================================
% Compute (DIM ONLY) + timing/memory
%% ============================================================
fprintf('\n==================== DIM VELOCITY ====================\n');
fprintf('z_target=%.3f, z_use=%.6f\n', z_target, z_use);
fprintf('obs grid: Nx=%d, Ny=%d, dx=%.4g m, r_boundary=%.4g m\n', numel(x_obs), numel(y_obs), dx_obs, r_boundary);

mem0 = get_mem_mb(); t0 = tic;
resV = calc_ultrasound_velocity_field(source, medium, calc, 'dim');
tAll = toc(t0); mem1 = get_mem_mb();

fprintf('DIM time: %.3f s\n', tAll);
fprintf('DIM memory: start %s, end %s, delta %s\n', fmt_mem(mem0), fmt_mem(mem1), fmt_mem(mem1-mem0));

%% -------------------- fetch fields (use f2) --------------------
X  = resV.dim.X;
Y  = resV.dim.Y;

vX = resV.dim.v_x_f2;
vY = resV.dim.v_y_f2;
vZ = resV.dim.v_z_f2;

S  = vX.^2 + vY.^2 + vZ.^2;                       % v·v (NO conjugate)
Vmag = sqrt(abs(vX).^2 + abs(vY).^2 + abs(vZ).^2);
VZph = angle(vZ)/pi;

%% -------------------- plot settings --------------------
lim_ph = [-1 1];

lim_vx = [min(abs(vX(:))), max(abs(vX(:)))];
lim_vy = [min(abs(vY(:))), max(abs(vY(:)))];
lim_vz = [min(abs(vZ(:))), max(abs(vZ(:)))];
lim_s  = [min(abs(S(:))),  max(abs(S(:)))];
lim_v  = [min(Vmag(:)),     max(Vmag(:))];

%% ============================================================
% Fig0: |v| + quiver(Re{v_x,v_y})  AND  phase(vz)/pi
%% ============================================================
do_quiver = true;

figure('Name', sprintf('|v| & phase(vz) @ z=%.3f m', z_use), ...
    'position',[100 100 1400 650], 'Color','w');
set(gcf,'Renderer','opengl');

% ---- |v| ----
subplot(1,2,1);
surf(X, Y, zeros(size(X)), Vmag, 'EdgeColor','none'); view(2);
shading interp; colormap(MyColor('vik'));
if all(isfinite(lim_v)) && lim_v(2) > lim_v(1), clim(lim_v); end
clb = colorbar; clb.Title.Interpreter='latex'; clb.Title.String='$|v|$';
set(clb,'Fontsize',18);

axis equal;
xlim([min(X(:)) max(X(:))]); ylim([min(Y(:)) max(Y(:))]);
set(gca,'linewidth',2,'TickLabelInterpreter','latex'); fontsize(gca,22,'points');
xlabel('$x$ (m)','Interpreter','latex'); ylabel('$y$ (m)','Interpreter','latex');
title('$|v|$', 'Interpreter','latex','Fontsize',18);

if do_quiver
    hold on;
    qstep_r = max(1, round(size(X,1)/25));
    qstep_t = max(1, round(size(X,2)/60));
    Xq  = X(1:qstep_r:end, 1:qstep_t:end);
    Yq  = Y(1:qstep_r:end, 1:qstep_t:end);
    vXq = real(vX(1:qstep_r:end, 1:qstep_t:end));
    vYq = real(vY(1:qstep_r:end, 1:qstep_t:end));
    sc  = max(hypot(vXq(:), vYq(:)));
    if sc > 0, vXq = vXq/sc; vYq = vYq/sc; end
    quiver(Xq, Yq, vXq, vYq, 0.6, 'k', 'LineWidth', 1.0);
    hold off;
end

% ---- phase(vz)/pi ----
subplot(1,2,2);
surf(X, Y, zeros(size(X)), VZph, 'EdgeColor','none'); view(2);
shading flat; colormap(hsv); clim(lim_ph);
clb = colorbar; clb.Title.Interpreter='latex'; clb.Title.String='$\angle v_z/\pi$';
set(clb,'Fontsize',18);

axis equal;
xlim([min(X(:)) max(X(:))]); ylim([min(Y(:)) max(Y(:))]);
set(gca,'linewidth',2,'TickLabelInterpreter','latex'); fontsize(gca,22,'points');
xlabel('$x$ (m)','Interpreter','latex'); ylabel('$y$ (m)','Interpreter','latex');
title('$\angle v_z/\pi$', 'Interpreter','latex','Fontsize',18);

sgtitle(sprintf('DIM, $f=%.1f$ kHz, $\\theta=%.1f^\\circ$, $n=%d$', source.f2/1e3, theta_deg, n), ...
    'Interpreter','latex','Fontsize',18);

local_save_fig_png(gcf, sprintf('DIM_Vmag_PhaseVz_z%.3f', z_use));

%% ============================================================
% Fig1..4: v_x, v_y, v_z, s=v·v  (mag/phase)
%% ============================================================
local_plot_mag_phase_cart(X, Y, vX, 'v_x', z_use, source.f2, theta_deg, lim_vx, lim_ph);
local_save_fig_png(gcf, sprintf('DIM_vx_magphase_z%.3f', z_use));

local_plot_mag_phase_cart(X, Y, vY, 'v_y', z_use, source.f2, theta_deg, lim_vy, lim_ph);
local_save_fig_png(gcf, sprintf('DIM_vy_magphase_z%.3f', z_use));

local_plot_mag_phase_cart(X, Y, vZ, 'v_z', z_use, source.f2, theta_deg, lim_vz, lim_ph);
local_save_fig_png(gcf, sprintf('DIM_vz_magphase_z%.3f', z_use));

local_plot_mag_phase_cart(X, Y, S, 's=v\cdot v', z_use, source.f2, theta_deg, lim_s, lim_ph, true);
local_save_fig_png(gcf, sprintf('DIM_s_vdotv_magphase_z%.3f', z_use));


%% ==================== local functions ====================

function local_plot_mag_phase_cart(X, Y, V, nameStr, z_use, f_hz, theta_deg, lim_mag, lim_ph, isS)
if nargin < 10, isS = false; end

figure('Color','w', 'Position',[100 100 1400 650], ...
    'Name',sprintf('%s (DIM) @ z=%.3f m', nameStr, z_use));
set(gcf,'Renderer','opengl');

% ----- magnitude (interp) -----
subplot(1,2,1);
surf(X, Y, zeros(size(X)), abs(V), 'EdgeColor','none'); view(2);
shading interp; colormap(MyColor('vik'));
if all(isfinite(lim_mag)) && lim_mag(2) > lim_mag(1), clim(lim_mag); end
clb = colorbar;
clb.Title.Interpreter = 'latex';
if ~isS
    clb.Title.String = sprintf('$|%s|$', strrep(nameStr,'_','\_'));
else
    clb.Title.String = '$|s|$';
end
set(clb,'Fontsize',18);

axis equal;
xlim([min(X(:)) max(X(:))]); ylim([min(Y(:)) max(Y(:))]);
set(gca,'linewidth',2,'TickLabelInterpreter','latex'); fontsize(gca,22,'points');
xlabel('$x$ (m)','Interpreter','latex'); ylabel('$y$ (m)','Interpreter','latex');
title(sprintf('Magnitude @ $z=%.3f$ m', z_use), 'Interpreter','latex','Fontsize',18);

% ----- phase/pi (flat) -----
subplot(1,2,2);
surf(X, Y, zeros(size(X)), angle(V)/pi, 'EdgeColor','none'); view(2);
shading flat; colormap(hsv);
clim(lim_ph);
clb = colorbar;
clb.Title.Interpreter = 'latex';
if ~isS
    clb.Title.String = sprintf('$\\angle %s/\\pi$', strrep(nameStr,'_','\_'));
else
    clb.Title.String = '$\angle s/\pi$';
end
set(clb,'Fontsize',18);

axis equal;
xlim([min(X(:)) max(X(:))]); ylim([min(Y(:)) max(Y(:))]);
set(gca,'linewidth',2,'TickLabelInterpreter','latex'); fontsize(gca,22,'points');
xlabel('$x$ (m)','Interpreter','latex'); ylabel('$y$ (m)','Interpreter','latex');
title('Phase', 'Interpreter','latex','Fontsize',18);

sgtitle(sprintf('DIM, $f=%.1f$ kHz, $\\theta=%.1f^\\circ$', f_hz/1e3, theta_deg), ...
    'Interpreter','latex','Fontsize',18);
end

function vr = local_vrho_rings(rho, r_edges, w_ring)
rho = abs(rho);
vr = zeros(size(rho));

q = discretize(rho, r_edges);   % 1..Q, NaN outside
mask = ~isnan(q);
if any(mask(:))
    vr(mask) = w_ring(q(mask));
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
