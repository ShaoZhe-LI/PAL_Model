%% ============================================================
% MAIN: DIM-only spiral-source (paper-like) velocity @ z = z_obs
% - Source: circular spiral binary phase grating (0/pi) as vibrating surface
% - Requires: make_source_velocity.m  (your updated version supports
%             source.custom_vn_xy_handle_1/2)
% - Calls:    calc_ultrasound_velocity_field(...,'dim')
%
% Plots (DIM, f2) on z=z_use plane:
%   0) |v| + quiver(Re{v_x,v_y})      + angle(vz)/pi
%   1) v_x : |v_x| + angle(v_x)/pi
%   2) v_y : |v_y| + angle(v_y)/pi
%   3) v_z : |v_z| + angle(v_z)/pi
%   4) s=v·v (NO conj): |s| + angle(s)/pi
%
% Colormap:
%   - magnitude: hsv + shading interp  (consistent with other components)
%   - phase:     hsv + shading interp  (avoid pixel grid; note phase jump artifacts possible)
% ============================================================

clear; clc; close all;

%% -------------------- save figures (GLOBAL control) --------------------
global SAVE_PNG SAVE_DIR
SAVE_PNG = false;

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

%% ============================================================
% SOURCE (paper parameters): circular spiral grating (binary phase)
% From the paper snippet (your page2.png):
%   M = 5 turns, r0 = 10 mm, f = 97 kHz, theta0 = pi/4 => P = sqrt(2)*lambda
%   dP = 0.45 P
% ============================================================
source = struct();
source.profile = 'Custom';

% frequencies (keep PAL-compatible; f2 used in DIM by default)
source.f1 = 97e3;
source.fa = 2e3;
source.f2 = source.f1 + source.fa;

% basic amplitude (surface normal velocity magnitude)
source.v0 = 0.172;
source.v_ratio = 1;

% IMPORTANT: for arbitrary 2D vn(x,y), set m_custom = 0 (do NOT add exp(j*m*phi))
source.m_custom = 0;

% wavelength at f1 (paper uses monochromatic 97 kHz)
lambda1 = medium.c0 / source.f1;

% grating parameters
theta0 = pi/4;              % cone angle in geometric grating design
P  = sqrt(2) * lambda1;      % from paper: theta0=pi/4 => P = sqrt(2)*lambda
dP = 0.45 * P;               % from paper: dP = 0.45 P
M  = 5;                      % turns
r0 = 10e-3;                  % minimum radius (m)

% aperture radius implied by spiral extent
a_spiral = r0 + M * P;
source.a = a_spiral;

% choose handedness: +1 (CCW) or -1 (CW)
handed = +1;

% Define the spiral binary phase grating as vn(x,y):
%   u = rho - r0 - b*phi, b = P/(2*pi)
%   mask = mod(u,P) < dP  => phase 0, else phase pi  (=> +1 / -1)
%   vn = v0 * (+1 or -1) inside spiral annulus; 0 outside
vn_handle = @(X,Y) local_spiral_binary_grating_vn(X,Y, ...
    source.v0, source.a, r0, P, dP, M, handed);

% Use same vn for f1 and f2 (paper is monochromatic; here keep consistent)
source.custom_vn_xy_handle_1 = vn_handle;
source.custom_vn_xy_handle_2 = @(X,Y) (source.v_ratio) * vn_handle(X,Y);

%% -------------------- calc --------------------
calc = struct();

% FHT fields not used, but make_source_velocity needs defaults
calc.fht.N_FHT      = 32768;
calc.fht.rho_max    = 0.25;
calc.fht.Nh_scale   = 1.2;
calc.fht.NH_scale   = 4;
calc.fht.Nh_v_scale = 1.1;
calc.fht.zu_max     = 0.2;
calc.fht.za_max     = 0.2;

% DIM discretization
calc.dim.use_freq = 'f2';
calc.dim.dis_coe  = 16;
calc.dim.margin   = 1;
calc.dim.src_discretization = 'polar';     % 'cart' or 'polar'

% DIM blocks
calc.dim.block_size     = 20000;
calc.dim.src_block_size = 5000;

%% ============================================================
% Target z from paper Eq.(3):
%   z_min = r0 * (P/lambda) * sqrt(1 - lambda^2/P^2)
%   z_max = z_min * (1 + M*P/r0)
%   z_obs = (z_min+z_max)/2
% With P = sqrt(2)*lambda => z_min = r0
% ============================================================
z_min = r0;                  % simplifies under P=sqrt(2)lambda
z_max = z_min * (1 + M*P/r0);
z_target = 0.5 * (z_min + z_max);

% Get consistent z grid from make_source_velocity (your function)
[~, fht_tmp] = make_source_velocity(source, medium, calc);
z_full = fht_tmp.z_ultra(:);
[~, iz] = min(abs(z_full - z_target));
z_use = z_full(iz);
calc.dim.z_use = z_use;

%% ============================================================
% Observation plane range (paper figures use scale bar ~ lambda)
% Use square with full width ~ lambda (so half-width = lambda/2)
% ============================================================
lambda2 = medium.c0 / source.f2;          % for reference
half_width = 0.5 * lambda1;              % match paper scale bar ~ lambda
r_boundary = 1.2 * half_width;           % small margin
dx_obs = half_width / 64;                   % resolution

x_obs = -r_boundary:dx_obs:r_boundary;
y_obs = -r_boundary:dx_obs:r_boundary;

calc.dim.obs_grid = struct();
calc.dim.obs_grid.x = x_obs;
calc.dim.obs_grid.y = y_obs;
calc.dim.obs_grid.z = z_use;

%% ============================================================
% Compute (DIM ONLY) + timing/memory
%% ============================================================
fprintf('\n==================== DIM VELOCITY (SPIRAL SOURCE) ====================\n');
fprintf('f1=%.1f kHz, f2=%.1f kHz, lambda1=%.4g m\n', source.f1/1e3, source.f2/1e3, lambda1);
fprintf('Spiral: M=%d, r0=%.3g m, P=%.3g m, dP=%.3g m, a=%.3g m, handed=%+d\n', ...
    M, r0, P, dP, source.a, handed);
fprintf('z_target=%.6f m, z_use=%.6f m\n', z_target, z_use);
fprintf('obs grid: Nx=%d, Ny=%d, dx=%.4g m, half-width=%.4g m\n', ...
    numel(x_obs), numel(y_obs), dx_obs, r_boundary);

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

%% -------------------- plot settings --------------------
lim_ph = [-1 1];

lim_vx = [min(abs(vX(:))), max(abs(vX(:)))];
lim_vy = [min(abs(vY(:))), max(abs(vY(:)))];
lim_vz = [min(abs(vZ(:))), max(abs(vZ(:)))];
lim_s  = [min(abs(S(:))),  max(abs(S(:)))];
lim_v  = [min(Vmag(:)),    max(Vmag(:))];

%% ============================================================
% Fig0: |v| + quiver(Re{v_x,v_y})  AND  angle(vz)/pi
%% ============================================================
do_quiver = true;

figure('Name', sprintf('|v| & angle(vz) @ z=%.3f m', z_use), ...
    'position',[100 100 1400 650], 'Color','w');
set(gcf,'Renderer','opengl');

% ---- |v| ----
subplot(1,2,1);
surf(X, Y, zeros(size(X)), Vmag, 'EdgeColor','none'); view(2);
shading interp; colormap(hsv);
if lim_v(2) > lim_v(1), clim(lim_v); end
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

% ---- angle(vz)/pi ----
subplot(1,2,2);
surf(X, Y, zeros(size(X)), angle(vZ)/pi, 'EdgeColor','none'); view(2);
shading interp; colormap(hsv); clim(lim_ph);
clb = colorbar; clb.Title.Interpreter='latex'; clb.Title.String='$\angle v_z/\pi$';
set(clb,'Fontsize',18);

axis equal;
xlim([min(X(:)) max(X(:))]); ylim([min(Y(:)) max(Y(:))]);
set(gca,'linewidth',2,'TickLabelInterpreter','latex'); fontsize(gca,22,'points');
xlabel('$x$ (m)','Interpreter','latex'); ylabel('$y$ (m)','Interpreter','latex');
title('$\angle v_z/\pi$', 'Interpreter','latex','Fontsize',18);

sgtitle(sprintf('Spiral source (paper-like), DIM @ $z=%.3f$ m, $f_2=%.1f$ kHz', z_use, source.f2/1e3), ...
    'Interpreter','latex','Fontsize',18);

local_save_fig_png(gcf, sprintf('DIM_Vmag_PhaseVz_spiral_z%.3f', z_use));

%% ============================================================
% Fig1..4: v_x, v_y, v_z, s=v·v  (mag/phase)
%% ============================================================
local_plot_mag_phase_cart(X, Y, vX, 'v_x', z_use, source.f2, lim_vx, lim_ph);
local_save_fig_png(gcf, sprintf('DIM_vx_magphase_spiral_z%.3f', z_use));

local_plot_mag_phase_cart(X, Y, vY, 'v_y', z_use, source.f2, lim_vy, lim_ph);
local_save_fig_png(gcf, sprintf('DIM_vy_magphase_spiral_z%.3f', z_use));

local_plot_mag_phase_cart(X, Y, vZ, 'v_z', z_use, source.f2, lim_vz, lim_ph);
local_save_fig_png(gcf, sprintf('DIM_vz_magphase_spiral_z%.3f', z_use));

local_plot_mag_phase_cart(X, Y, S, 's=v\cdot v', z_use, source.f2, lim_s, lim_ph, true);
local_save_fig_png(gcf, sprintf('DIM_s_vdotv_magphase_spiral_z%.3f', z_use));

%% ==================== local functions ====================

function local_plot_mag_phase_cart(X, Y, V, nameStr, z_use, f_hz, lim_mag, lim_ph, isS)
if nargin < 9, isS = false; end

figure('Color','w', 'Position',[100 100 1400 650], ...
    'Name',sprintf('%s (DIM) @ z=%.3f m', nameStr, z_use));
set(gcf,'Renderer','opengl');

% ----- magnitude -----
subplot(1,2,1);
surf(X, Y, zeros(size(X)), abs(V), 'EdgeColor','none'); view(2);
shading interp; colormap(hsv);
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

% ----- phase/pi -----
subplot(1,2,2);
surf(X, Y, zeros(size(X)), angle(V)/pi, 'EdgeColor','none'); view(2);
shading interp; colormap(hsv); clim(lim_ph);
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

sgtitle(sprintf('DIM, $f=%.1f$ kHz', f_hz/1e3), 'Interpreter','latex','Fontsize',18);
end

function vn = local_spiral_binary_grating_vn(X,Y,v0,a,r0,P,dP,M,handed)
% Binary phase spiral grating (0/pi) on an annulus [r0, r0+M*P]
% handed: +1 CCW, -1 CW (controls spiral chirality)

rho = hypot(X,Y);
phi = atan2(Y,X);                    % [-pi,pi]

vn = zeros(size(rho));

r1 = r0 + M*P;
mask_ap = (rho <= a) & (rho >= r0) & (rho <= r1);

if ~any(mask_ap(:)), return; end

b = P/(2*pi);

% Spiral coordinate (Archimedean): rho = r0 + b*(handed*phi) + k*P
u = rho(mask_ap) - r0 - b*(handed*phi(mask_ap));

t = mod(u, P);                       % [0,P)

% duty cycle: dP/P
on = (t < dP);

% two phase levels: 0 or pi => +1 or -1
g = ones(size(t));
g(~on) = -1;

vn(mask_ap) = v0 * g;
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
