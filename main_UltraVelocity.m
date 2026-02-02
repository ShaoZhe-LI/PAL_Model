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
%   - calc_ultrasound_velocity_field.m   (incl. DIM block)
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
source.m = 1;
source.F = 0.2;

source.f1 = 42e3;
source.fa = 2e3;

%% -------------------- calc parameters --------------------
calc = struct();

% --- FHT / King ---
calc.fht.N_FHT = 32768;
calc.fht.rho_max = 0.2;
calc.fht.Nh_scale = 1.2;
calc.fht.NH_scale = 4;
calc.fht.Nh_v_scale = 1.1;
calc.fht.zu_max = 1.1;
calc.fht.za_max = 1;

% --- DIM source discretization ---
calc.dim.use_freq = 'f2';
calc.dim.dis_coe = 16;
calc.dim.margin = 1;
calc.dim.src_discretization = 'polar'; % 'cart' or 'polar'

% --- King analytic spectrum stability ---
calc.king.gspec_method = 'analytic';
calc.king.eps_kzz   = 1e-3;
calc.king.eps_phase = calc.king.eps_kzz;
calc.king.kz_min    = 1e-12;

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
calc.dim.obs_grid.z = z_use;        % scalar

calc.dim.block_size     = 20000;
calc.dim.src_block_size = 5000;

%% ============================================================
% Compute velocity field (King + DIM)
% ============================================================
fprintf('\n==================== VELOCITY (KING + DIM) ====================\n');
mem0 = get_mem_mb(); t0 = tic;

resV = calc_ultrasound_velocity_field(source, medium, calc, 'both');

tAll = toc(t0); mem1 = get_mem_mb();
fprintf('Velocity time (both): %.3f s\n', tAll);
fprintf('Velocity memory: start %s, end %s, delta %s\n', fmt_mem(mem0), fmt_mem(mem1), fmt_mem(mem1-mem0));

%% ============================================================
% Extract King (f1) on same z_use, downsample to rho_ds
% ============================================================
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
% ============================================================
% enforce consistent indexing with obs_grid
Nx_line = numel(calc.dim.obs_grid.x);
if size(resV.dim.v_rho_f1,1) ~= 1 || size(resV.dim.v_rho_f1,2) ~= Nx_line
    error('DIM output size mismatch: expected (Ny=1, Nx=%d).', Nx_line);
end

vR_dim = reshape(resV.dim.v_rho_f1(1,1:Nx_line), [], 1);
vP_dim = reshape(resV.dim.v_phi_f1(1,1:Nx_line), [], 1);
vZ_dim = reshape(resV.dim.v_z_f1  (1,1:Nx_line), [], 1);

%% ============================================================
% Phase & errors (exclude rho=0 to avoid singular/ill-conditioned points)
% ============================================================
rho_line = rho_ds(:);
mask_ok = rho_line > 0;          % drop rho=0 point
rho_p = rho_line(mask_ok);

vRk = vR_king(mask_ok); vRd = vR_dim(mask_ok);
vPk = vP_king(mask_ok); vPd = vP_dim(mask_ok);
vZk = vZ_king(mask_ok); vZd = vZ_dim(mask_ok);

if fig.unwrap
    phR_k = unwrap(angle(vRk)); phR_d = unwrap(angle(vRd));
    phP_k = unwrap(angle(vPk)); phP_d = unwrap(angle(vPd));
    phZ_k = unwrap(angle(vZk)); phZ_d = unwrap(angle(vZd));
else
    phR_k = angle(vRk); phR_d = angle(vRd);
    phP_k = angle(vPk); phP_d = angle(vPd);
    phZ_k = angle(vZk); phZ_d = angle(vZd);
end

eps0 = 1e-12;
errR_log = log10( abs(vRd - vRk) ./ (abs(vRk) + eps0) );
errP_log = log10( abs(vPd - vPk) ./ (abs(vPk) + eps0) );
errZ_log = log10( abs(vZd - vZk) ./ (abs(vZk) + eps0) );

%% ============================================================
% 1D FIGS: compare v_rho / v_phi / v_z  (3 subplots each) — style aligned
% ============================================================

% -------- helper: unify axes style --------
apply_axes_style = @(ax) set(ax,'LineWidth',2,'TickLabelInterpreter','latex');

% -------------------- v_rho --------------------
figure('Name',sprintf('v_rho: King vs DIM @ z=%.3f m', z_use), 'position',[100 100 1200 900]);

subplot(3,1,1);
plot(rho_p, abs(vRk), 'LineWidth',1.8); hold on;
plot(rho_p, abs(vRd), '--', 'LineWidth',1.8);
grid on; xlabel('\rho (m)','Interpreter','latex'); ylabel('$|v_\rho|$ (m/s)','Interpreter','latex');
title(sprintf('Magnitude: $v_\\rho$ @ $z=%.3f$ m', z_use), 'Interpreter','latex');
legend('King','DIM','Location','best');
apply_axes_style(gca); fontsize(gca,22,'points');

subplot(3,1,2);
plot(rho_p, phR_k, 'LineWidth',1.8); hold on;
plot(rho_p, phR_d, '--', 'LineWidth',1.8);
grid on; xlabel('\rho (m)','Interpreter','latex'); ylabel('Phase (rad)','Interpreter','latex');
title('Phase', 'Interpreter','latex');
legend('King','DIM','Location','best');
apply_axes_style(gca); fontsize(gca,22,'points');

subplot(3,1,3);
plot(rho_p, errR_log, 'LineWidth',1.8);
grid on; xlabel('\rho (m)','Interpreter','latex'); ylabel('$\log_{10}$ relative error','Interpreter','latex');
title('$\log_{10}(|v_{DIM}-v_{King}|/(|v_{King}|+\epsilon))$', 'Interpreter','latex');
apply_axes_style(gca); fontsize(gca,22,'points');

local_save_fig_png(gcf, sprintf('vRho_1D_compare3_z%.3fm', z_use));

% -------------------- v_phi --------------------
figure('Name',sprintf('v_phi: King vs DIM @ z=%.3f m', z_use), 'position',[120 120 1200 900]);

subplot(3,1,1);
plot(rho_p, abs(vPk), 'LineWidth',1.8); hold on;
plot(rho_p, abs(vPd), '--', 'LineWidth',1.8);
grid on; xlabel('\rho (m)','Interpreter','latex'); ylabel('$|v_\phi|$ (m/s)','Interpreter','latex');
title(sprintf('Magnitude: $v_\\phi$ @ $z=%.3f$ m', z_use), 'Interpreter','latex');
legend('King','DIM','Location','best');
apply_axes_style(gca); fontsize(gca,22,'points');

subplot(3,1,2);
plot(rho_p, phP_k, 'LineWidth',1.8); hold on;
plot(rho_p, phP_d, '--', 'LineWidth',1.8);
grid on; xlabel('\rho (m)','Interpreter','latex'); ylabel('Phase (rad)','Interpreter','latex');
title('Phase', 'Interpreter','latex');
legend('King','DIM','Location','best');
apply_axes_style(gca); fontsize(gca,22,'points');

subplot(3,1,3);
plot(rho_p, errP_log, 'LineWidth',1.8);
grid on; xlabel('\rho (m)','Interpreter','latex'); ylabel('$\log_{10}$ relative error','Interpreter','latex');
title('$\log_{10}(|v_{DIM}-v_{King}|/(|v_{King}|+\epsilon))$', 'Interpreter','latex');
apply_axes_style(gca); fontsize(gca,22,'points');

local_save_fig_png(gcf, sprintf('vPhi_1D_compare3_z%.3fm', z_use));

% -------------------- v_z --------------------
figure('Name',sprintf('v_z: King vs DIM @ z=%.3f m', z_use), 'position',[140 140 1200 900]);

subplot(3,1,1);
plot(rho_p, abs(vZk), 'LineWidth',1.8); hold on;
plot(rho_p, abs(vZd), '--', 'LineWidth',1.8);
grid on; xlabel('\rho (m)','Interpreter','latex'); ylabel('$|v_z|$ (m/s)','Interpreter','latex');
title(sprintf('Magnitude: $v_z$ @ $z=%.3f$ m', z_use), 'Interpreter','latex');
legend('King','DIM','Location','best');
apply_axes_style(gca); fontsize(gca,22,'points');

subplot(3,1,2);
plot(rho_p, phZ_k, 'LineWidth',1.8); hold on;
plot(rho_p, phZ_d, '--', 'LineWidth',1.8);
grid on; xlabel('\rho (m)','Interpreter','latex'); ylabel('Phase (rad)','Interpreter','latex');
title('Phase', 'Interpreter','latex');
legend('King','DIM','Location','best');
apply_axes_style(gca); fontsize(gca,22,'points');

subplot(3,1,3);
plot(rho_p, errZ_log, 'LineWidth',1.8);
grid on; xlabel('\rho (m)','Interpreter','latex'); ylabel('$\log_{10}$ relative error','Interpreter','latex');
title('$\log_{10}(|v_{DIM}-v_{King}|/(|v_{King}|+\epsilon))$', 'Interpreter','latex');
apply_axes_style(gca); fontsize(gca,22,'points');

local_save_fig_png(gcf, sprintf('vZ_1D_compare3_z%.3fm', z_use));


%% ============================================================
% 2D xOy @ z≈1 m (King + DIM-phase-extended) — side-by-side compare
% One figure per metric (SPL-style layout):
%   Fig1: |v|  (interp)    King | DIM
%   Fig2: phase(vz)/pi (flat) King | DIM
%   Fig3: s=v·v  |s| (interp) and angle(s)/pi (flat), King row + DIM row
%
% Seam fix: surf(view(2)) + theta-wrap (duplicate first column)
% Magnitude: shading interp
% Phase:     shading flat
% Short titles (avoid truncation)
%% ============================================================

r_boundary = 0.04;

% theta grids (linspace)
Ntheta_k = 721;  theta_k = linspace(0, 2*pi, Ntheta_k);
Ntheta_d = 181;  theta_d = linspace(0, 2*pi, Ntheta_d);

m_use = resV.king.m;

% ---- crop radial samples ----
rho_fig = rho_line(:).';
idx_rb  = find(rho_fig <= r_boundary, 1, 'last');
if isempty(idx_rb), idx_rb = numel(rho_fig); end
rho_fig = rho_fig(1:idx_rb);

% King 1D lines
vRk = vR_king(:).'; vRk = vRk(1:idx_rb);
vPk = vP_king(:).'; vPk = vPk(1:idx_rb);
vZk = vZ_king(:).'; vZk = vZk(1:idx_rb);

% DIM 1D lines
vRd = vR_dim(:).';  vRd = vRd(1:idx_rb);
vPd = vP_dim(:).';  vPd = vPd(1:idx_rb);
vZd = vZ_dim(:).';  vZd = vZd(1:idx_rb);

%% ============================================================
% Build KING 2D (double)
%% ============================================================
[THk, Rk] = meshgrid(theta_k, rho_fig);
[Xk, Yk]  = pol2cart(THk, Rk);
Ephi_k = exp(1i * m_use * THk);

vR2k = (vRk(:) * ones(1, numel(theta_k))) .* Ephi_k;
vP2k = (vPk(:) * ones(1, numel(theta_k))) .* Ephi_k;
vZ2k = (vZk(:) * ones(1, numel(theta_k))) .* Ephi_k;

cTk = cos(THk); sTk = sin(THk);
vX2k = vR2k .* cTk - vP2k .* sTk;
vY2k = vR2k .* sTk + vP2k .* cTk;

Vmag_k = sqrt(abs(vX2k).^2 + abs(vY2k).^2 + abs(vZ2k).^2);
VZph_k = angle(vZ2k) / pi;

S_k    = vX2k.^2 + vY2k.^2 + vZ2k.^2;
Smag_k = abs(S_k);
Sph_k  = angle(S_k) / pi;

% wrap theta to kill seam
Xkw     = [Xk, Xk(:,1)];
Ykw     = [Yk, Yk(:,1)];
Vmag_kw = [Vmag_k, Vmag_k(:,1)];
VZph_kw = [VZph_k, VZph_k(:,1)];
Smag_kw = [Smag_k, Smag_k(:,1)];
Sph_kw  = [Sph_k,  Sph_k(:,1)];

% for quiver
vX_kw = [vX2k, vX2k(:,1)];
vY_kw = [vY2k, vY2k(:,1)];
vZ_kw = [vZ2k, vZ2k(:,1)];

%% ============================================================
% Build DIM 2D (single -> double for plotting)
%% ============================================================
[THd, Rd] = meshgrid(single(theta_d), single(rho_fig));
[Xd, Yd]  = pol2cart(double(THd), double(Rd));

Ephi_d = exp(1i * single(m_use) .* THd);

vR2d = (single(vRd(:)) * ones(1, size(THd,2), 'single')) .* Ephi_d;
vP2d = (single(vPd(:)) * ones(1, size(THd,2), 'single')) .* Ephi_d;
vZ2d = (single(vZd(:)) * ones(1, size(THd,2), 'single')) .* Ephi_d;

cTd = cos(THd); sTd = sin(THd);
vX2d = vR2d .* cTd - vP2d .* sTd;
vY2d = vR2d .* sTd + vP2d .* cTd;

Vmag_d = sqrt(abs(vX2d).^2 + abs(vY2d).^2 + abs(vZ2d).^2);
VZph_d = angle(vZ2d) / pi;

S_d    = vX2d.^2 + vY2d.^2 + vZ2d.^2;
Smag_d = abs(S_d);
Sph_d  = angle(S_d) / pi;

% wrap theta to kill seam
Xdw     = [Xd, Xd(:,1)];
Ydw     = [Yd, Yd(:,1)];
Vmag_dw = [double(Vmag_d), double(Vmag_d(:,1))];
VZph_dw = [double(VZph_d), double(VZph_d(:,1))];
Smag_dw = [double(Smag_d), double(Smag_d(:,1))];
Sph_dw  = [double(Sph_d),  double(Sph_d(:,1))];

vX_dw = [double(vX2d), double(vX2d(:,1))];
vY_dw = [double(vY2d), double(vY2d(:,1))];
vZ_dw = [double(vZ2d), double(vZ2d(:,1))];

%% ============================================================
% Unified limits (so King/DIM comparable)
%% ============================================================
lim_Vmag = [min([Vmag_kw(:); Vmag_dw(:)]), max([Vmag_kw(:); Vmag_dw(:)])];
lim_Smag = [min([Smag_kw(:); Smag_dw(:)]), max([Smag_kw(:); Smag_dw(:)])];
lim_ph   = [-1 1];

%% ============================================================
% Fig1: |v| + quiver (King | DIM)
%% ============================================================
figure('Name',sprintf('|v| (z=%.2f)', z_use), 'position',[100 100 1400 650]);
set(gcf,'Renderer','opengl');

% ---- King ----
subplot(1,2,1);
surf(Xkw, Ykw, zeros(size(Xkw)), Vmag_kw, 'EdgeColor','none'); view(2);
shading interp; colormap(MyColor('vik')); clim(lim_Vmag);
clb = colorbar; clb.Title.Interpreter='latex'; clb.Title.String='$|v|$'; set(clb,'Fontsize',18);
axis equal; xlim(1.2*r_boundary*[-1 1]); ylim(1.2*r_boundary*[-1 1]);
set(gca,'linewidth',2,'TickLabelInterpreter','latex'); fontsize(gca,22,'points');
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
title('King','Interpreter','latex','Fontsize',18);

hold on;
qstep_r = max(1, round(size(Xkw,1)/25));
qstep_t = max(1, round(size(Xkw,2)/60));
Xq  = Xkw(1:qstep_r:end, 1:qstep_t:end);
Yq  = Ykw(1:qstep_r:end, 1:qstep_t:end);
vXq = real(vX_kw(1:qstep_r:end, 1:qstep_t:end));
vYq = real(vY_kw(1:qstep_r:end, 1:qstep_t:end));
sc  = max(hypot(vXq(:), vYq(:)));
if sc > 0, vXq = vXq/sc; vYq = vYq/sc; end
quiver(Xq, Yq, vXq, vYq, 0.6, 'k', 'LineWidth', 1.0);
hold off;

% ---- DIM ----
subplot(1,2,2);
surf(Xdw, Ydw, zeros(size(Xdw)), Vmag_dw, 'EdgeColor','none'); view(2);
shading interp; colormap(MyColor('vik')); clim(lim_Vmag);
clb = colorbar; clb.Title.Interpreter='latex'; clb.Title.String='$|v|$'; set(clb,'Fontsize',18);
axis equal; xlim(1.2*r_boundary*[-1 1]); ylim(1.2*r_boundary*[-1 1]);
set(gca,'linewidth',2,'TickLabelInterpreter','latex'); fontsize(gca,22,'points');
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
title('DIM','Interpreter','latex','Fontsize',18);

hold on;
qstep_r = max(1, round(size(Xdw,1)/25));
qstep_t = max(1, round(size(Xdw,2)/60));
Xq  = Xdw(1:qstep_r:end, 1:qstep_t:end);
Yq  = Ydw(1:qstep_r:end, 1:qstep_t:end);
vXq = real(vX_dw(1:qstep_r:end, 1:qstep_t:end));
vYq = real(vY_dw(1:qstep_r:end, 1:qstep_t:end));
sc  = max(hypot(vXq(:), vYq(:)));
if sc > 0, vXq = vXq/sc; vYq = vYq/sc; end
quiver(Xq, Yq, vXq, vYq, 0.6, 'k', 'LineWidth', 1.0);
hold off;

sgtitle(sprintf('$|v|$, $m=%d$, $r\\le%.2f$', m_use, r_boundary), ...
    'Interpreter','latex','Fontsize',18);

local_save_fig_png(gcf, sprintf('cmp2D_Vmag_z%.2f_m%d', z_use, m_use));

%% ============================================================
% Fig2: phase(vz)/pi (King | DIM)
%% ============================================================
figure('Name',sprintf('phase(vz) (z=%.2f)', z_use), 'position',[120 120 1400 650]);
set(gcf,'Renderer','opengl');

subplot(1,2,1);
surf(Xkw, Ykw, zeros(size(Xkw)), VZph_kw, 'EdgeColor','none'); view(2);
shading flat; colormap('hsv'); clim(lim_ph);
clb = colorbar; clb.Title.Interpreter='latex'; clb.Title.String='$\angle v_z/\pi$'; set(clb,'Fontsize',18);
axis equal; xlim(1.2*r_boundary*[-1 1]); ylim(1.2*r_boundary*[-1 1]);
set(gca,'linewidth',2,'TickLabelInterpreter','latex'); fontsize(gca,22,'points');
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
title('King','Interpreter','latex','Fontsize',18);

subplot(1,2,2);
surf(Xdw, Ydw, zeros(size(Xdw)), VZph_dw, 'EdgeColor','none'); view(2);
shading flat; colormap('hsv'); clim(lim_ph);
clb = colorbar; clb.Title.Interpreter='latex'; clb.Title.String='$\angle v_z/\pi$'; set(clb,'Fontsize',18);
axis equal; xlim(1.2*r_boundary*[-1 1]); ylim(1.2*r_boundary*[-1 1]);
set(gca,'linewidth',2,'TickLabelInterpreter','latex'); fontsize(gca,22,'points');
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
title('DIM','Interpreter','latex','Fontsize',18);

sgtitle(sprintf('$\\angle v_z/\\pi$, $m=%d$, $r\\le%.2f$', m_use, r_boundary), ...
    'Interpreter','latex','Fontsize',18);

local_save_fig_png(gcf, sprintf('cmp2D_PhaseVz_z%.2f_m%d', z_use, m_use));

%% ============================================================
% Fig3: s=v·v  (2x2): |s| and angle(s)/pi, King row + DIM row
%% ============================================================
figure('Name',sprintf('s=v·v (z=%.2f)', z_use), 'position',[140 140 1400 900]);
set(gcf,'Renderer','opengl');

% ---- King |s| ----
subplot(2,2,1);
surf(Xkw, Ykw, zeros(size(Xkw)), Smag_kw, 'EdgeColor','none'); view(2);
shading interp; colormap(MyColor('vik')); clim(lim_Smag);
clb = colorbar; clb.Title.Interpreter='latex'; clb.Title.String='$|s|$'; set(clb,'Fontsize',18);
axis equal; xlim(1.2*r_boundary*[-1 1]); ylim(1.2*r_boundary*[-1 1]);
set(gca,'linewidth',2,'TickLabelInterpreter','latex'); fontsize(gca,22,'points');
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
title('King $|s|$','Interpreter','latex','Fontsize',18);

% ---- King angle(s)/pi ----
subplot(2,2,2);
surf(Xkw, Ykw, zeros(size(Xkw)), Sph_kw, 'EdgeColor','none'); view(2);
shading flat; colormap('hsv'); clim(lim_ph);
clb = colorbar; clb.Title.Interpreter='latex'; clb.Title.String='$\angle s/\pi$'; set(clb,'Fontsize',18);
axis equal; xlim(1.2*r_boundary*[-1 1]); ylim(1.2*r_boundary*[-1 1]);
set(gca,'linewidth',2,'TickLabelInterpreter','latex'); fontsize(gca,22,'points');
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
title('King $\angle s/\pi$','Interpreter','latex','Fontsize',18);

% ---- DIM |s| ----
subplot(2,2,3);
surf(Xdw, Ydw, zeros(size(Xdw)), Smag_dw, 'EdgeColor','none'); view(2);
shading interp; colormap(MyColor('vik')); clim(lim_Smag);
clb = colorbar; clb.Title.Interpreter='latex'; clb.Title.String='$|s|$'; set(clb,'Fontsize',18);
axis equal; xlim(1.2*r_boundary*[-1 1]); ylim(1.2*r_boundary*[-1 1]);
set(gca,'linewidth',2,'TickLabelInterpreter','latex'); fontsize(gca,22,'points');
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
title('DIM $|s|$','Interpreter','latex','Fontsize',18);

% ---- DIM angle(s)/pi ----
subplot(2,2,4);
surf(Xdw, Ydw, zeros(size(Xdw)), Sph_dw, 'EdgeColor','none'); view(2);
shading flat; colormap('hsv'); clim(lim_ph);
clb = colorbar; clb.Title.Interpreter='latex'; clb.Title.String='$\angle s/\pi$'; set(clb,'Fontsize',18);
axis equal; xlim(1.2*r_boundary*[-1 1]); ylim(1.2*r_boundary*[-1 1]);
set(gca,'linewidth',2,'TickLabelInterpreter','latex'); fontsize(gca,22,'points');
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
title('DIM $\angle s/\pi$','Interpreter','latex','Fontsize',18);

sgtitle(sprintf('$s=v\\cdot v$, $m=%d$, $r\\le%.2f$', m_use, r_boundary), ...
    'Interpreter','latex','Fontsize',18);

local_save_fig_png(gcf, sprintf('cmp2D_vdotv_z%.2f_m%d', z_use, m_use));

%% ---- free DIM big vars (optional) ----
clear vR2d vP2d vZ2d vX2d vY2d Vmag_d VZph_d S_d Smag_d Sph_d THd Rd Ephi_d cTd sTd
clear Xd Yd Xdw Ydw Vmag_dw VZph_dw Smag_dw Sph_dw vX_dw vY_dw vZ_dw














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
