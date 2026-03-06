%% ============================================================
% MAIN: DIM-only spiral-source (paper-like) velocity @ z = z_obs
% - Source: circular spiral binary phase grating (0/pi) as vibrating surface
% - Requires: make_source_velocity.m  (supports source.custom_vn_xy_handle_1/2)
% - Calls:    calc_ultrasound_velocity_field(...,'dim')
%
% Figures:
%   Fig0: |v| + quiver(Re{vx,vy})   AND   angle(vz)/pi
%   Fig1: vx  (|vx| + angle(vx)/pi)
%   Fig2: vy  (|vy| + angle(vy)/pi)
%   Fig3: vz  (|vz| + angle(vz)/pi)
%   Fig4: s=v·v (NO conj) (|s| + angle(s)/pi)
%
% + (ADDED) Fig5-like (NOW IN 3D):
%   - detect C-points from S=v·v phase vortices (winding number)
%   - compute major-axis direction (3D) and fix sign by instantaneous velocity at phase phi_sign
%   - build two closed contours (C1: enclose 1 point, C2: enclose selected points)
%   - plot 3D bi-vectors (double-ended arrows) using quiver3 on the plane z=z_use
%
% Colormap:
%   - magnitude: hsv + shading interp
%   - phase:     hsv + shading flat
% ============================================================

% 本代码箭头方向为长轴方向与初始时刻振速方向确定
% 圈内每存在一个C点则发生一次翻转，奇数个则整体呈现翻转，偶数个则整体不翻转

clear; clc; close all;

%% -------------------- save switches --------------------
global SAVE_PNG SAVE_MAT SAVE_DIR
SAVE_PNG = true;        % save figures
SAVE_MAT = SAVE_PNG;    % save .mat

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
% SOURCE (paper parameters)
%% ============================================================
source = struct();
source.profile = 'Custom';

% frequencies (PAL-compatible; DIM uses f2 by default)
source.f1 = 97e3;
source.fa = 2e3;
source.f2 = source.f1 + source.fa;

% amplitude (surface normal velocity magnitude)
source.v0 = 0.172;
source.v_ratio = 1;

% IMPORTANT: for arbitrary 2D vn(x,y), set m_custom = 0
source.m_custom = 0;

% wavelength at f1
lambda1 = medium.c0 / source.f1;

% grating parameters
theta0 = pi/4; %#ok<NASGU>
P  = sqrt(2) * lambda1;
dP = 0.45 * P;
M  = 5;
r0 = 10e-3;

% aperture radius implied by spiral extent
source.a = r0 + M * P;

% handedness: +1 CCW, -1 CW
handed = +1;

% ===== Bessel/OAM order (YOU CHANGE THIS) =====
l = 1;     % <<<<<< 改这里：1/2/3/...  贝塞尔阶数（拓扑荷）
r_boundary_coef = 1.8;  % 2D图像范围
r_small_k = 4;          % 小圆
r_large_k = 16;         % 大圆

extra = struct();
extra.l = l;
extra.r_boundary_coef = r_boundary_coef;
extra.r_small_k = r_small_k;
extra.r_large_k = r_large_k;

% Define vn(x,y)
vn_handle = @(X,Y) local_spiral_binary_grating_vn(X,Y, ...
    source.v0, source.a, r0, P, dP, M, handed, l);

source.custom_vn_xy_handle_1 = vn_handle;
source.custom_vn_xy_handle_2 = @(X,Y) source.v_ratio * vn_handle(X,Y);

%% -------------------- calc --------------------
calc = struct();

% FHT not used, but make_source_velocity needs defaults
calc.fht.N_FHT      = 32768;
calc.fht.rho_max    = 0.25;
calc.fht.Nh_scale   = 1.2;
calc.fht.NH_scale   = 4;
calc.fht.Nh_v_scale = 1.1;
calc.fht.zu_max     = 0.2;
calc.fht.za_max     = 0.2;

% DIM
calc.dim.use_freq = 'f2';
calc.dim.dis_coe  = 16;
calc.dim.margin   = 1;
calc.dim.src_discretization = 'polar';   % 'cart' or 'polar'
calc.dim.block_size     = 20000;
calc.dim.src_block_size = 5000;

%% -------------------- fig (for runinfo compatibility) --------------------
fig = struct();
fig.unwrap = false;

%% -------------------- save setup --------------------
if (SAVE_PNG || SAVE_MAT)
    tstr  = datestr(datetime('now'), 'mmdd_HHMM');
    f_str = sprintf('f1=%.0fk', source.f1/1e3);
    l_str = sprintf('l=%d', l);
    a_str = sprintf('a=%.1fmm', source.a*1e3);
    SAVE_DIR = sprintf('LiteratureReview__%s__%s__%s__%s', f_str, a_str, l_str, tstr);

    if ~exist(SAVE_DIR, 'dir'), mkdir(SAVE_DIR); end
    local_write_runinfo_txt(SAVE_DIR, medium, source, calc, fig, extra);
else
    SAVE_DIR = '';
end

%% ============================================================
% Target z from paper Eq.(3)
%% ============================================================
z_min = r0;
z_max = z_min * (1 + M*P/r0);
z_target = 0.5 * (z_min + z_max);

% Get consistent z grid
[~, fht_tmp] = make_source_velocity(source, medium, calc);
z_full = fht_tmp.z_ultra(:);
[~, iz] = min(abs(z_full - z_target));
z_use = z_full(iz);
calc.dim.z_use = z_use;

%% ============================================================
% Observation plane range (paper-like: scale ~ lambda)
%% ============================================================
half_width = 0.5 * lambda1;
r_boundary = r_boundary_coef * half_width;
dx_obs = half_width / 32;

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
fprintf('order l=%d, handed=%+d\n', l, handed);
fprintf('f1=%.1f kHz, f2=%.1f kHz, lambda1=%.4g m\n', source.f1/1e3, source.f2/1e3, lambda1);
fprintf('Spiral: M=%d, r0=%.3g m, P=%.3g m, dP=%.3g m, a=%.3g m\n', M, r0, P, dP, source.a);
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

S    = vX.^2 + vY.^2 + vZ.^2;                 % v·v (NO conjugate)
Vmag = sqrt(abs(vX).^2 + abs(vY).^2 + abs(vZ).^2);

%% -------------------- save data (.mat) --------------------
if SAVE_MAT && ~isempty(SAVE_DIR)
    data_file = fullfile(SAVE_DIR, 'resV_spiral_DIM.mat');

    meta = struct();
    meta.time_str   = datestr(datetime('now'), 'yyyy-mm-dd HH:MM:SS');
    meta.l          = l;
    meta.handed     = handed;
    meta.M          = M;
    meta.r0         = r0;
    meta.P          = P;
    meta.dP         = dP;
    meta.z_target   = z_target;
    meta.z_use      = z_use;
    meta.dx_obs     = dx_obs;
    meta.r_boundary = r_boundary;
    meta.tAll       = tAll;
    meta.mem0_MB    = mem0;
    meta.mem1_MB    = mem1;

    save(data_file, ...
        'medium','source','calc','meta', ...
        'resV','X','Y','vX','vY','vZ','S','Vmag', ...
        '-v7.3');

    fprintf('Saved MAT: %s\n', data_file);
end

%% ============================================================
% (ADDED) C-point finding + Fig5-like bi-vector plot (3D)
%% ============================================================
do_fig5_like = true;

if do_fig5_like
    fprintf('\n==================== C-POINT DETECTION (from S=v·v) ====================\n');

    Smag = abs(S);
    phS  = angle(S);

    phi_sign = 0;   % u(t0)=Re{v*exp(-j*phi_sign)} to set sign

    % ---- detection options ----
    copt = struct();
    copt.max_candidates = 25;
    copt.smag_rel_th    = 0.1;
    copt.min_sep_pix    = 4;
    copt.ring_radius    = 2;

    Cpts = local_find_Cpoints_from_S(X, Y, Smag, phS, copt);  % [x,y,q,|S|]

    fprintf('Candidates (x,y,charge,|S|): %d\n', size(Cpts,1));
    for k = 1:size(Cpts,1)
        fprintf('  #%02d: x=%.6g, y=%.6g, q=%+d, |S|=%.3g\n', ...
            k, Cpts(k,1), Cpts(k,2), Cpts(k,3), Cpts(k,4));
    end

    % ---- robust loop charge + pick set ----
    Cnz_all = Cpts(Cpts(:,3)~=0, :);

    dxy = mean_grid_step(X,Y);
    r_small = r_small_k * dxy;
    r_large = r_large_k * dxy;

    if ~isempty(Cnz_all)
        [~,ii0] = min(Cnz_all(:,4));
        p0 = Cnz_all(ii0,1:2);
    else
        p0 = [mean(X(:)), mean(Y(:))];
    end

    rGamma = max(10*dxy, 0.35*max(range(X(:)), range(Y(:))));
    qGamma = local_charge_along_circle(p0, rGamma, @(x,y) local_interp_phase_S(X,Y,phS,x,y));
    qExp = 2*l;

    fprintf('\n[Fig5-like] Total charge on loop Γ: qGamma=%+d, expected ~ %+d (2*l)\n', qGamma, qExp);

    Cnz_in = Cnz_all;
    if ~isempty(Cnz_in)
        rr = hypot(Cnz_in(:,1)-p0(1), Cnz_in(:,2)-p0(2));
        Cnz_in = Cnz_in(rr <= rGamma, :);
    end
    sum_in = sum(Cnz_in(:,3));
    fprintf('[Fig5-like] Inside Γ: %d nonzero-charge candidates, sum(q)=%+d\n', size(Cnz_in,1), sum_in);

    % pick set
    picked = zeros(0,4);    % [x y q |S|]
    if ~isempty(Cnz_in)
        [~,ii] = sort(Cnz_in(:,4),'ascend');
        Cnz_in = Cnz_in(ii,:);

        qsum = 0;
        for k = 1:min(20,size(Cnz_in,1))
            qk = Cnz_in(k,3);

            % keep common low orders
            if ~ismember(abs(qk), [1 2 4])
                continue;
            end

            picked(end+1,:) = Cnz_in(k,:); %#ok<AGROW>
            qsum = qsum + qk;

            if qsum == qGamma
                break;
            end
        end
        fprintf('[Fig5-like] Picked %d points, sum(q)=%+d (target qGamma=%+d)\n', size(picked,1), qsum, qGamma);
    end

    if ~isempty(picked)
        p1 = picked(1,1:2);
    else
        p1 = p0;
    end

    % contours
    C1 = local_make_circle(p1, r_small, 70);
    if ~isempty(picked) && size(picked,1) >= 2
        pc = mean(picked(:,1:2), 1);
        rr = max(vecnorm(picked(:,1:2) - pc, 2, 2)) + r_small;
        C2 = local_make_circle(pc, max(rr, r_large), 90);
    else
        C2 = local_make_circle(p1, r_large, 90);
    end

    % ---- interpolants for v ----
    is_meshgrid = ...
        size(X,1) > 1 && size(X,2) > 1 && ...
        all(abs(diff(X(:,1))) < 1e-12) && ...
        all(abs(diff(Y(1,:))) < 1e-12);

    if is_meshgrid
        Xn  = X.';   Yn  = Y.';
        vXn = vX.';  vYn = vY.';  vZn = vZ.';
        Fvx = griddedInterpolant(Xn, Yn, vXn, 'linear', 'none');
        Fvy = griddedInterpolant(Xn, Yn, vYn, 'linear', 'none');
        Fvz = griddedInterpolant(Xn, Yn, vZn, 'linear', 'none');
    else
        Fvx = griddedInterpolant(X, Y, vX, 'linear', 'none');
        Fvy = griddedInterpolant(X, Y, vY, 'linear', 'none');
        Fvz = griddedInterpolant(X, Y, vZ, 'linear', 'none');
    end

    % ---- major-axis (SIGNED, 3D) on contours ----
    E1 = local_major_axis_directors_on_contour(C1, Fvx, Fvy, Fvz, phi_sign);
    E2 = local_major_axis_directors_on_contour(C2, Fvx, Fvy, Fvz, phi_sign);

    E1 = local_norm_rows(E1);
    E2 = local_norm_rows(E2);

    % ---- plot Fig5-like in 3D ----
    col_pos = [0 0.45 1.00];
    col_neg = [1.00 0.35 0.35];

    qscale3 = 0.25;     % arrow length scale (3D)
    z0 = z_use;         % all points lie on plane z=z_use

    figure('Color','w','Position',[120 120 1200 560], ...
        'Name',sprintf('Fig5-like bi-vectors (3D) @ z=%.4g m', z_use));
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    % --- left: C1 ---
    nexttile; hold on; box on; grid on;
    plot3(C1(:,1), C1(:,2), z0*ones(size(C1,1),1), 'k', 'LineWidth', 1.0);
    local_plot_bivectors3d(C1, E1, z0, qscale3, col_pos, col_neg);
    plot3(p1(1), p1(2), z0, 'k.', 'MarkerSize', 18);
    title('$C_1$ (3D bi-vectors)','Interpreter','latex');
    xlabel('$x$ (m)','Interpreter','latex'); ylabel('$y$ (m)','Interpreter','latex'); zlabel('$z$ (m)','Interpreter','latex');
    axis equal;
    view(35,25);

    % --- right: C2 ---
    nexttile; hold on; box on; grid on;
    plot3(C2(:,1), C2(:,2), z0*ones(size(C2,1),1), 'k', 'LineWidth', 1.0);
    % show C1 range on C2 panel
    plot3(C1(:,1), C1(:,2), z0*ones(size(C1,1),1), '--', 'LineWidth', 1.6, 'Color', [0 0.7 0]);
    local_plot_bivectors3d(C2, E2, z0, qscale3, col_pos, col_neg);
    if ~isempty(picked)
        plot3(picked(:,1), picked(:,2), z0*ones(size(picked,1),1), 'k.', 'MarkerSize', 18);
    else
        plot3(p1(1), p1(2), z0, 'k.', 'MarkerSize', 18);
    end
    title('$C_2$ (3D bi-vectors)','Interpreter','latex');
    xlabel('$x$ (m)','Interpreter','latex'); ylabel('$y$ (m)','Interpreter','latex'); zlabel('$z$ (m)','Interpreter','latex');
    axis equal;
    view(35,25);

    local_save_fig_png(gcf, sprintf('DIM_Fig5_bivector3D_spiral_z%.4f_l%d', z_use, l));
    local_save_fig_fig(gcf, sprintf('DIM_Fig5_bivector3D_spiral_z%.4f_l%d', z_use, l));
end

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
shading flat; colormap(hsv); clim(lim_ph);
clb = colorbar; clb.Title.Interpreter='latex'; clb.Title.String='$\angle v_z/\pi$';
set(clb,'Fontsize',18);

axis equal;
xlim([min(X(:)) max(X(:))]); ylim([min(Y(:)) max(Y(:))]);
set(gca,'linewidth',2,'TickLabelInterpreter','latex'); fontsize(gca,22,'points');
xlabel('$x$ (m)','Interpreter','latex'); ylabel('$y$ (m)','Interpreter','latex');
title('$\angle v_z/\pi$', 'Interpreter','latex','Fontsize',18);

sgtitle(sprintf('Spiral source, DIM @ $z=%.3f$ m, $f_2=%.1f$ kHz, $l=%d$', ...
    z_use, source.f2/1e3, l), 'Interpreter','latex','Fontsize',18);

local_save_fig_png(gcf, sprintf('DIM_Vmag_PhaseVz_spiral_z%.4f_l%d', z_use, l));

%% ============================================================
% Fig1..4: v_x, v_y, v_z, s=v·v  (mag/phase)
%% ============================================================
local_plot_mag_phase_cart(X, Y, vX, 'v_x', z_use, source.f2, lim_vx, lim_ph);
local_save_fig_png(gcf, sprintf('DIM_vx_magphase_spiral_z%.4f_l%d', z_use, l));

local_plot_mag_phase_cart(X, Y, vY, 'v_y', z_use, source.f2, lim_vy, lim_ph);
local_save_fig_png(gcf, sprintf('DIM_vy_magphase_spiral_z%.4f_l%d', z_use, l));

local_plot_mag_phase_cart(X, Y, vZ, 'v_z', z_use, source.f2, lim_vz, lim_ph);
local_save_fig_png(gcf, sprintf('DIM_vz_magphase_spiral_z%.4f_l%d', z_use, l));

local_plot_mag_phase_cart(X, Y, S, 's=v\cdot v', z_use, source.f2, lim_s, lim_ph, true);
local_save_fig_png(gcf, sprintf('DIM_s_vdotv_magphase_spiral_z%.4f_l%d', z_use, l));

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
shading flat; colormap(hsv); clim(lim_ph);
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

sgtitle(sprintf('DIM, $f=%.1f$ kHz, $z=%.3f$ m', f_hz/1e3, z_use), ...
    'Interpreter','latex','Fontsize',18);
end

function vn = local_spiral_binary_grating_vn(X,Y,v0,a,r0,P,dP,M,handed,l)
rho = hypot(X,Y);
phi = atan2(Y,X);

vn = zeros(size(rho));

r1 = r0 + M*P;
mask_ap = (rho <= a) & (rho >= r0) & (rho <= r1);
if ~any(mask_ap(:)), return; end

b = P/(2*pi);
phi_eff = l * phi(mask_ap);

u = rho(mask_ap) - r0 - b*(handed*phi_eff);
t = mod(u, P);

on = (t < dP);
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


function local_save_fig_fig(fig_handle, base_name)
global SAVE_PNG SAVE_DIR
if isempty(SAVE_PNG) || ~SAVE_PNG, return; end
if isempty(SAVE_DIR) || ~exist(SAVE_DIR,'dir'), return; end

fname = local_sanitize_filename(base_name);
fp = fullfile(SAVE_DIR, [fname, '.fig']);
set(fig_handle, 'Color', 'w');
try
    savefig(fig_handle, fp);                % <-- 用 savefig 保存 .fig
catch
    % 兼容旧版本 MATLAB
    try
        saveas(fig_handle, fp);
    catch
        warning('Failed to save FIG: %s', fp);
    end
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

function local_write_runinfo_txt(save_dir, medium, source, calc, fig, extra)
fp = fullfile(save_dir, 'run_info.txt');
fid = fopen(fp, 'w');
if fid < 0
    warning('Cannot create run_info.txt in %s', save_dir);
    return;
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, '===== RUN INFO =====\n');
fprintf(fid, 'Time: %s\n\n', datestr(datetime('now'), 'yyyy-mm-dd HH:MM:SS'));

fprintf(fid, '%% medium\n');
fprintf(fid, 'medium.c0=%.15g; medium.rho0=%.15g; medium.beta=%.15g; medium.pref=%.15g;\n', ...
    medium.c0, medium.rho0, medium.beta, medium.pref);
fprintf(fid, 'medium.use_absorp=%s;\n\n', local_bool_str(medium.use_absorp));

fprintf(fid, '%% source\n');
fprintf(fid, 'source.profile=''%s''; source.a=%.15g; source.v0=%.15g; source.v_ratio=%.15g;\n', ...
    char(string(source.profile)), source.a, source.v0, source.v_ratio);
if isfield(source,'m_custom')
    fprintf(fid, 'source.m_custom=%d;\n', source.m_custom);
end
fprintf(fid, 'source.f1=%.15g; source.f2=%.15g; source.fa=%.15g;\n\n', source.f1, source.f2, source.fa);

fprintf(fid, '%% calc.fht\n');
cf = calc.fht;
fprintf(fid, 'N_FHT=%d; rho_max=%.15g; Nh_scale=%.15g; NH_scale=%.15g; Nh_v_scale=%.15g; zu_max=%.15g; za_max=%.15g;\n\n', ...
    cf.N_FHT, cf.rho_max, cf.Nh_scale, cf.NH_scale, cf.Nh_v_scale, cf.zu_max, cf.za_max);

fprintf(fid, '%% calc.dim\n');
cd = calc.dim;
fprintf(fid, 'use_freq=''%s''; dis_coe=%.15g; margin=%.15g; src_discretization=''%s'';\n', ...
    char(string(cd.use_freq)), cd.dis_coe, cd.margin, char(string(cd.src_discretization)));

fprintf(fid, '\n%% fig\n');
fprintf(fid, 'fig.unwrap=%s;\n', local_bool_str(fig.unwrap));

fprintf(fid, '\n%% extra (for reproduction)\n');
fprintf(fid, 'l=%d;\n', extra.l);
fprintf(fid, 'r_boundary_coef=%.15g;\n', extra.r_boundary_coef);
fprintf(fid, 'r_small_k=%d;  %% r_small = r_small_k * dxy\n', extra.r_small_k);
fprintf(fid, 'r_large_k=%d;  %% r_large = r_large_k * dxy\n', extra.r_large_k);

fprintf(fid, '\n===== END =====\n');
end

function s = local_bool_str(x)
if logical(x), s = 'true'; else, s = 'false'; end
end

function dxy = mean_grid_step(X,Y)
dx = abs(X(1,2)-X(1,1));
dy = abs(Y(2,1)-Y(1,1));
dxy = mean([dx,dy]);
end

function Cpts = local_find_Cpoints_from_S(X, Y, Smag, phS, opts)
if ~isfield(opts,'max_candidates'), opts.max_candidates = 20; end
if ~isfield(opts,'min_sep_pix'),    opts.min_sep_pix    = 4;  end
if ~isfield(opts,'ring_radius'),    opts.ring_radius    = 2;  end
if ~isfield(opts,'smag_rel_th'),    opts.smag_rel_th    = 0.06; end

medS = median(Smag(:), 'omitnan');
th   = opts.smag_rel_th * medS;

A = Smag;
A(~isfinite(A)) = inf;

[vals, idx] = sort(A(:), 'ascend');
idx = idx(vals < th);
if isempty(idx)
    Cpts = zeros(0,4); return;
end
idx = idx(1:min(numel(idx), 5*opts.max_candidates));

[ny,nx] = size(A);
[iy,ix] = ind2sub([ny,nx], idx);

supp = false(size(ix));
keep_list = [];

for k = 1:numel(ix)
    if numel(keep_list) >= opts.max_candidates, break; end
    if ~supp(k)
        keep_list(end+1) = k; %#ok<AGROW>
        dx = ix - ix(k);
        dy = iy - iy(k);
        near = (dx.^2 + dy.^2) <= opts.min_sep_pix^2;
        supp(near) = true;
    end
end

Cpts = zeros(numel(keep_list),4);
for kk = 1:numel(keep_list)
    k = keep_list(kk);
    q = local_phase_winding_charge(phS, ix(k), iy(k), opts.ring_radius);
    Cpts(kk,:) = [X(iy(k),ix(k)), Y(iy(k),ix(k)), q, Smag(iy(k),ix(k))];
end
end

function q = local_phase_winding_charge(ph, ix, iy, r)
[ny,nx] = size(ph);
if ix-r < 2 || ix+r > nx-1 || iy-r < 2 || iy+r > ny-1
    q = 0; return;
end

pts = [];
for x = ix-r:ix+r,       pts(end+1,:) = [x,   iy-r]; end %#ok<AGROW>
for y = iy-r+1:iy+r,     pts(end+1,:) = [ix+r,y  ]; end %#ok<AGROW>
for x = ix+r-1:-1:ix-r,  pts(end+1,:) = [x,   iy+r]; end %#ok<AGROW>
for y = iy+r-1:-1:iy-r+1,pts(end+1,:) = [ix-r,y  ]; end %#ok<AGROW>

ang = zeros(size(pts,1),1);
for k = 1:size(pts,1)
    ang(k) = ph(pts(k,2), pts(k,1));
end

d = diff([ang; ang(1)]);
d = atan2(sin(d), cos(d));
w = sum(d);
q = round(w/(2*pi));
end

function C = local_make_circle(center_xy, radius, Nsamp)
tt = linspace(0,2*pi,Nsamp+1);   % include endpoint
C  = [center_xy(1)+radius*cos(tt(:)), center_xy(2)+radius*sin(tt(:))];
end

function E = local_major_axis_directors_on_contour(C, Fvx, Fvy, Fvz, phi_sign)
C_use = C;
if size(C,1) >= 2 && norm(C(end,:) - C(1,:)) < 1e-12
    C_use = C(1:end-1,:);
end

Ns = size(C_use,1);
E  = nan(Ns,3);
for i = 1:Ns
    x = C_use(i,1); y = C_use(i,2);
    vx = Fvx(x,y); vy = Fvy(x,y); vz = Fvz(x,y);
    if any(~isfinite([vx,vy,vz])), continue; end
    E(i,:) = local_major_axis_vector_signed([vx;vy;vz], phi_sign).';
end
end

function e = local_major_axis_director(v)
a = real(v(:));
b = imag(v(:));
M = a*a.' + b*b.';
[V,D] = eig(M);
[~,id] = max(diag(D));
e = V(:,id);
n = norm(e);
if n > 0, e = e/n; end
end

function e = local_major_axis_vector_signed(v, phi_sign)
e = local_major_axis_director(v);
u0 = real(v(:) * exp(-1j*phi_sign));
if norm(u0) > 1e-12
    if dot(e, u0) < 0
        e = -e;
    end
end
end

function local_plot_bivectors3d(C, E, z0, qscale, col_pos, col_neg)
% Plot bi-vectors (±E) along a closed contour C on plane z=z0
% C may be explicitly closed (last==first); E is computed on C(1:end-1,:)

C_use = C;
if size(C,1) >= 2 && norm(C(end,:) - C(1,:)) < 1e-12
    C_use = C(1:end-1,:);
end

NsC = size(C_use,1);
NsE = size(E,1);
Ns  = min(NsC, NsE);

x  = C_use(1:Ns,1);
y  = C_use(1:Ns,2);
z  = z0*ones(Ns,1);

ex = E(1:Ns,1);
ey = E(1:Ns,2);
ez = E(1:Ns,3);

ok = all(isfinite([x,y,ex,ey,ez]),2);
x=x(ok); y=y(ok); z=z(ok);
ex=ex(ok); ey=ey(ok); ez=ez(ok);

quiver3(x,y,z,  ex, ey, ez, qscale, 'Color',col_pos,'LineWidth',1.2,'MaxHeadSize',0.9);
quiver3(x,y,z, -ex,-ey,-ez, qscale, 'Color',col_neg,'LineWidth',1.2,'MaxHeadSize',0.9);
end

function q = local_charge_along_circle(center_xy, radius, phase_fun)
Ns = 200;
tt = linspace(0,2*pi,Ns+1); tt(end)=[];
x  = center_xy(1) + radius*cos(tt);
y  = center_xy(2) + radius*sin(tt);

ang = zeros(Ns,1);
for k = 1:Ns
    ang(k) = phase_fun(x(k), y(k));
    if ~isfinite(ang(k)), ang(k) = 0; end
end

d = diff([ang; ang(1)]);
d = atan2(sin(d), cos(d));
w = sum(d);
q = round(w/(2*pi));
end

function ph = local_interp_phase_S(X,Y,phS,xq,yq)
persistent Fph is_init
if isempty(is_init) || isempty(Fph)
    is_init = true;

    is_meshgrid = ...
        size(X,1) > 1 && size(X,2) > 1 && ...
        all(abs(diff(X(:,1))) < 1e-12) && ...
        all(abs(diff(Y(1,:))) < 1e-12);

    if is_meshgrid
        Xn = X.'; Yn = Y.'; phn = phS.';
        Fph = griddedInterpolant(Xn, Yn, phn, 'linear', 'none');
    else
        Fph = griddedInterpolant(X, Y, phS, 'linear', 'none');
    end
end
ph = Fph(xq, yq);
end

function A = local_norm_rows(A)
nn = sqrt(sum(A.^2,2));
A = A ./ max(nn, eps);
end