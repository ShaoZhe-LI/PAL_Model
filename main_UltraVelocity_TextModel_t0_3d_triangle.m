%% ============================================================
% MAIN: DIM-only spiral-source (paper-like) velocity @ z = z_obs
% - Source: circular spiral binary phase grating (0/pi) as vibrating surface
% - Requires: make_source_velocity.m  (supports source.custom_vn_xy_handle_1/2)
% - Calls:    calc_ultrasound_velocity_field(..., obs_grid_v, 'dim')
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
% + (ADDED) Pressure field:
%   - DIM ultrasonic pressure p at the same z_use
%   - pressure observation plane is 3x larger in side length than velocity plane
%   - plot |p| and angle(p)/pi and save
%
% Colormap:
%   - magnitude: vik + shading interp
%   - phase:     hsv + shading flat
%
% NOTE:
%   - All PHYSICS/CALCULATION still use SI unit: meter (m)
%   - All PLOTTING axes are displayed in millimeter (mm)
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
source.f1 = 95e3;
source.fa = 2e3;
source.f2 = source.f1 + source.fa;

% amplitude (surface normal velocity magnitude)
source.v0 = 0.172;
source.v_ratio = 1;

% IMPORTANT: for arbitrary 2D vn(x,y), set m_custom = 0
source.m_custom = 0;

% source mode:
%   'spiral_abnormal' : 原始写法，径向等宽、法向不等宽
%   'spiral_normal'   : 近似按螺旋局部法向等宽
%   'triangle_normal' : 三角形
source.mode = 'triangle_normal';

% wavelength at f1/f2
lambda1 = medium.c0 / source.f2;

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
r_boundary_coef = 1.2;  % 2D图像范围（速度场）
r_small_k = 4;          % 小圆
r_large_k = 16;         % 大圆
qk_vec = 1:10;

%% ============================================================
% Observation plane range (paper-like: scale ~ lambda)
%% ============================================================

% ---- range mode ----
use_custom_obs_range = true;   % true: 用自定义范围; false: 保持原来的对称 r_boundary 范围

% ---- custom ranges (unit: m) ----
xlim_v_custom = [-2.7e-3,  4.9e-3];
ylim_v_custom = [-5.1e-3,  2.4e-3];

xlim_p_custom = [-2.6e-3, 11.8e-3];
ylim_p_custom = [-9.0e-3,  5.1e-3];

if strcmp(source.mode,'triangle_normal')
    source.f1   = 68e3;
    source.fa   = 2e3;
    source.f2   = source.f1 + source.fa;

    theta0 = pi/4;
    lambda1 = medium.c0 / source.f2;
    P  = sqrt(2) * lambda1;
    dP = 0.45 * P;

    M  = 8;
    r0 = 22e-3;

    source.a = 1.1*(r0 + 2*M*P + dP);
    handed = -1;
    l = -1;
end

extra = struct();
extra.l = l;
extra.r_boundary_coef = r_boundary_coef;
extra.r_small_k = r_small_k;
extra.r_large_k = r_large_k;

% Define vn(x,y)
vn_handle = @(X,Y) local_spiral_binary_grating_vn(X,Y, ...
    source.v0, source.a, r0, P, dP, M, handed, l, source.mode);

source.custom_vn_xy_handle_1 = vn_handle;
source.custom_vn_xy_handle_2 = @(X,Y) source.v_ratio * vn_handle(X,Y);

%% -------------------- detection options --------------------
copt = struct();
copt.max_candidates = 80;
copt.smag_rel_th    = 0.2;
copt.min_sep_pix    = 2;
copt.ring_radius    = 1;

%% -------------------- preview source pattern --------------------
src_plot_N = 801;
src_margin = 1.15;

x_src = linspace(-src_margin*source.a, src_margin*source.a, src_plot_N);
y_src = linspace(-src_margin*source.a, src_margin*source.a, src_plot_N);
[Xs, Ys] = meshgrid(x_src, y_src);

Vsrc = vn_handle(Xs, Ys);
mask_disk = hypot(Xs, Ys) <= source.a;
Vsrc(~mask_disk) = NaN;

Xsrc_mm = 1e3 * Xs;
Ysrc_mm = 1e3 * Ys;

fig_source = figure('Color','w', 'Position',[120 120 760 680], ...
    'Name', sprintf('Source pattern: %s', char(string(source.mode))));
ax = axes; %#ok<NASGU>
surf(Xsrc_mm, Ysrc_mm, zeros(size(Xs)), Vsrc, 'EdgeColor','none');
view(2);
shading flat;
axis equal;
xlim(1e3 * src_margin * [-source.a, source.a]);
ylim(1e3 * src_margin * [-source.a, source.a]);
xlabel('$x$ (mm)','Interpreter','latex');
ylabel('$y$ (mm)','Interpreter','latex');
title(sprintf('Source pattern: %s', strrep(char(string(source.mode)),'_','\_')), ...
    'Interpreter','none');

colormap(parula);
clb = colorbar;
clb.Title.String = '$v_n$';
clb.Title.Interpreter = 'latex';

set(gca,'LineWidth',1.5,'TickLabelInterpreter','latex');
fontsize(gca,20,'points');

hold on;
tt = linspace(0,2*pi,600);
plot(1e3*source.a*cos(tt), 1e3*source.a*sin(tt), 'k-', 'LineWidth', 1.0);
hold off;

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

% --- DIM pressure method for calc_ultrasound_field ---
calc.dim.method = 'rayleigh';   % 'rayleigh' or 'asm'
calc.dim.use_parallel = true;
calc.dim.num_workers  = 4;

%% -------------------- fig (for runinfo compatibility) --------------------
fig = struct();
fig.unwrap = false;

%% -------------------- save setup --------------------
if (SAVE_PNG || SAVE_MAT)
    tstr  = datestr(datetime('now'), 'mmdd_HHMM');
    f_str = sprintf('f=%.0fk', source.f2/1e3);
    l_str = sprintf('l=%d', l);
    a_str = sprintf('a=%.1fmm', source.a*1e3);
    SAVE_DIR = sprintf('LiteratureReview__%s__%s__%s__%s', f_str, a_str, l_str, tstr);

    if ~exist(SAVE_DIR, 'dir'), mkdir(SAVE_DIR); end
    local_write_runinfo_txt(SAVE_DIR, medium, source, calc, fig, extra, copt);
else
    SAVE_DIR = '';
end

%% ============================================================
% Target z from paper Eq.(3)
%% ============================================================
z_min = r0;
z_max = z_min * (1 + M*P/r0);
z_target = 0.5 * (z_min + z_max);
z_use = z_target;

%% ============================================================
% Observation plane range (paper-like: scale ~ lambda)
%% ============================================================
half_width = 0.5 * lambda1;
r_boundary = r_boundary_coef * half_width;
dx_obs = half_width / 32;

if use_custom_obs_range
    % ---- custom velocity-field observation grid ----
    x_obs = xlim_v_custom(1):dx_obs:xlim_v_custom(2);
    y_obs = ylim_v_custom(1):dx_obs:ylim_v_custom(2);

    % ---- custom pressure-field observation grid ----
    dx_obs_p = dx_obs;   % keep same spacing as before
    x_obs_p = xlim_p_custom(1):dx_obs_p:xlim_p_custom(2);
    y_obs_p = ylim_p_custom(1):dx_obs_p:ylim_p_custom(2);

    % keep variables for compatibility / logging
    r_boundary_p = NaN;
    press_range_scale = NaN;
else
    % ---- original symmetric ranges ----
    x_obs = -r_boundary:dx_obs:r_boundary;
    y_obs = -r_boundary:dx_obs:r_boundary;

    press_range_scale = 3;
    r_boundary_p = press_range_scale * r_boundary;
    dx_obs_p     = dx_obs;   % keep same sampling step

    x_obs_p = -r_boundary_p:dx_obs_p:r_boundary_p;
    y_obs_p = -r_boundary_p:dx_obs_p:r_boundary_p;
end

% velocity-field observation grid
obs_grid_v = struct();
obs_grid_v.dim = struct();
obs_grid_v.dim.x = x_obs;
obs_grid_v.dim.y = y_obs;
obs_grid_v.dim.z = z_use;
obs_grid_v.dim.block_size     = calc.dim.block_size;
obs_grid_v.dim.src_block_size = calc.dim.src_block_size;

% pressure-field observation grid
obs_grid_p = struct();
obs_grid_p.dim = struct();
obs_grid_p.dim.x = x_obs_p;
obs_grid_p.dim.y = y_obs_p;
obs_grid_p.dim.z = z_use;
obs_grid_p.dim.block_size     = calc.dim.block_size;
obs_grid_p.dim.src_block_size = calc.dim.src_block_size;

% ---- record observation-grid info into extra ----
extra.use_custom_obs_range = use_custom_obs_range;
extra.z_use   = z_use;

extra.xlim_v  = [x_obs(1),   x_obs(end)];
extra.ylim_v  = [y_obs(1),   y_obs(end)];
extra.xlim_p  = [x_obs_p(1), x_obs_p(end)];
extra.ylim_p  = [y_obs_p(1), y_obs_p(end)];

extra.Nx_v = numel(x_obs);
extra.Ny_v = numel(y_obs);
extra.Nx_p = numel(x_obs_p);
extra.Ny_p = numel(y_obs_p);

if numel(x_obs) >= 2,   extra.dx_obs_v = x_obs(2)   - x_obs(1);   else, extra.dx_obs_v = NaN; end
if numel(y_obs) >= 2,   extra.dy_obs_v = y_obs(2)   - y_obs(1);   else, extra.dy_obs_v = NaN; end
if numel(x_obs_p) >= 2, extra.dx_obs_p = x_obs_p(2) - x_obs_p(1); else, extra.dx_obs_p = NaN; end
if numel(y_obs_p) >= 2, extra.dy_obs_p = y_obs_p(2) - y_obs_p(1); else, extra.dy_obs_p = NaN; end

%% ============================================================
% Compute velocity field (DIM ONLY) + timing/memory
%% ============================================================
fprintf('\n==================== DIM VELOCITY (SPIRAL SOURCE) ====================\n');
fprintf('order l=%d, handed=%+d\n', l, handed);
fprintf('f1=%.1f kHz, f2=%.1f kHz, lambda1=%.4g m\n', source.f1/1e3, source.f2/1e3, lambda1);
fprintf('Spiral: M=%d, r0=%.3g m, P=%.3g m, dP=%.3g m, a=%.3g m\n', M, r0, P, dP, source.a);
fprintf('z_target=%.6f m, z_use=%.6f m\n', z_target, z_use);
if use_custom_obs_range
    fprintf('velocity obs grid: Nx=%d, Ny=%d, dx=%.4g m, x=[%.4g, %.4g] m, y=[%.4g, %.4g] m\n', ...
        numel(x_obs), numel(y_obs), dx_obs, x_obs(1), x_obs(end), y_obs(1), y_obs(end));
else
    fprintf('velocity obs grid: Nx=%d, Ny=%d, dx=%.4g m, half-width=%.4g m\n', ...
        numel(x_obs), numel(y_obs), dx_obs, r_boundary);
end

mem0 = get_mem_mb(); t0 = tic;
resV = calc_ultrasound_velocity_field(source, medium, calc, obs_grid_v, 'dim');
tAll = toc(t0); mem1 = get_mem_mb();

fprintf('DIM velocity time: %.3f s\n', tAll);
fprintf('DIM velocity memory: start %s, end %s, delta %s\n', fmt_mem(mem0), fmt_mem(mem1), fmt_mem(mem1-mem0));

%% ============================================================
% Compute ultrasonic pressure field (DIM ONLY, same z_use)
%% ============================================================
fprintf('\n==================== DIM PRESSURE (SPIRAL SOURCE) ====================\n');
if use_custom_obs_range
    fprintf('pressure obs grid: Nx=%d, Ny=%d, dx=%.4g m, x=[%.4g, %.4g] m, y=[%.4g, %.4g] m\n', ...
        numel(x_obs_p), numel(y_obs_p), dx_obs_p, x_obs_p(1), x_obs_p(end), y_obs_p(1), y_obs_p(end));
else
    fprintf('pressure obs grid: Nx=%d, Ny=%d, dx=%.4g m, half-width=%.4g m\n', ...
        numel(x_obs_p), numel(y_obs_p), dx_obs_p, r_boundary_p);
end

memP0 = get_mem_mb(); tP0 = tic;
resP = calc_ultrasound_field(source, medium, calc, obs_grid_p, 'dim');
tP = toc(tP0); memP1 = get_mem_mb();

fprintf('DIM pressure time: %.3f s\n', tP);
fprintf('DIM pressure memory: start %s, end %s, delta %s\n', ...
    fmt_mem(memP0), fmt_mem(memP1), fmt_mem(memP1-memP0));

%% -------------------- fetch fields (use f2) --------------------
X  = resV.dim.X;
Y  = resV.dim.Y;

% velocity DIM output is Ny x Nx x Nz; here Nz=1
if ndims(resV.dim.v_x_f2) == 3
    vX = resV.dim.v_x_f2(:,:,1);
    vY = resV.dim.v_y_f2(:,:,1);
    vZ = resV.dim.v_z_f2(:,:,1);
else
    vX = resV.dim.v_x_f2;
    vY = resV.dim.v_y_f2;
    vZ = resV.dim.v_z_f2;
end

% ultrasonic pressure field at f2
if isfield(resP.dim,'X') && isfield(resP.dim,'Y')
    Xp = resP.dim.X;
    Yp = resP.dim.Y;
else
    [Xp, Yp] = meshgrid(resP.dim.x, resP.dim.y);
end

if ndims(resP.dim.p_f2) == 3
    pU = resP.dim.p_f2(:,:,1);
else
    pU = resP.dim.p_f2;
end

Xpmm = 1e3 * Xp;
Ypmm = 1e3 * Yp;

Xmm = 1e3 * X;
Ymm = 1e3 * Y;
z0mm = 1e3 * z_use;

fprintf('Velocity X range = [%.3f, %.3f] mm\n', min(Xmm(:)), max(Xmm(:)));
fprintf('Velocity Y range = [%.3f, %.3f] mm\n', min(Ymm(:)), max(Ymm(:)));
fprintf('Pressure X range = [%.3f, %.3f] mm\n', min(Xpmm(:)), max(Xpmm(:)));
fprintf('Pressure Y range = [%.3f, %.3f] mm\n', min(Ypmm(:)), max(Ypmm(:)));

S    = vX.^2 + vY.^2 + vZ.^2;
Vmag = sqrt(abs(vX).^2 + abs(vY).^2 + abs(vZ).^2);

%% -------------------- save data (.mat) --------------------
if SAVE_MAT && ~isempty(SAVE_DIR)
    data_file = fullfile(SAVE_DIR, 'resV_spiral_DIM.mat');

    meta = struct();
    meta.time_str     = datestr(datetime('now'), 'yyyy-mm-dd HH:MM:SS');
    meta.l            = l;
    meta.handed       = handed;
    meta.M            = M;
    meta.r0           = r0;
    meta.P            = P;
    meta.dP           = dP;
    meta.z_target     = z_target;
    meta.z_use        = z_use;

    meta.dx_obs_v     = dx_obs;
    meta.r_boundary_v = r_boundary;

    meta.press_range_scale = press_range_scale;
    meta.dx_obs_p     = dx_obs_p;
    meta.r_boundary_p = r_boundary_p;

    meta.t_vel        = tAll;
    meta.mem0_vel_MB  = mem0;
    meta.mem1_vel_MB  = mem1;

    meta.t_pres       = tP;
    meta.mem0_pres_MB = memP0;
    meta.mem1_pres_MB = memP1;

    save(data_file, ...
        'medium','source','calc','obs_grid_v','obs_grid_p','meta', ...
        'resV','resP', ...
        'X','Y','vX','vY','vZ','S','Vmag', ...
        'Xp','Yp','pU', ...
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

    phi_sign = 0;

    Cpts = local_find_Cpoints_from_S(X, Y, Smag, phS, copt);

    fprintf('Candidates (x,y,charge,|S|): %d\n', size(Cpts,1));
    for k = 1:size(Cpts,1)
        fprintf('  #%02d: x=%.6g, y=%.6g, q=%+d, |S|=%.3g\n', ...
            k, Cpts(k,1), Cpts(k,2), Cpts(k,3), Cpts(k,4));
    end

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

    picked = zeros(0,4);
    if ~isempty(Cnz_in)
        [~,ii] = sort(Cnz_in(:,4),'ascend');
        Cnz_in = Cnz_in(ii,:);

        qsum = 0;
        for k = 1:min(20,size(Cnz_in,1))
            qk = Cnz_in(k,3);

            if ~ismember(abs(qk), qk_vec)
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

    C1 = local_make_circle(p1, r_small, 70);
    if ~isempty(picked) && size(picked,1) >= 2
        pc = mean(picked(:,1:2), 1);
        rr = max(vecnorm(picked(:,1:2) - pc, 2, 2)) + r_small;
        C2 = local_make_circle(pc, max(rr, r_large), 90);
    else
        C2 = local_make_circle(p1, r_large, 90);
    end

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

    E1 = local_major_axis_directors_on_contour(C1, Fvx, Fvy, Fvz, phi_sign);
    E2 = local_major_axis_directors_on_contour(C2, Fvx, Fvy, Fvz, phi_sign);

    E1 = local_norm_rows(E1);
    E2 = local_norm_rows(E2);

    col_pos = [0 0.45 1.00];
    col_neg = [1.00 0.35 0.35];

    qscale3 = 0.25;
    z0 = z0mm;

    figure('Color','w','Position',[120 120 1200 560], ...
        'Name',sprintf('Fig5-like bi-vectors (3D) @ z=%.4g mm', z0mm));
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    nexttile; hold on; box on; grid on;
    plot3(1e3*C1(:,1), 1e3*C1(:,2), z0*ones(size(C1,1),1), 'k', 'LineWidth', 1.0);
    local_plot_bivectors3d(C1, E1, z0, qscale3, col_pos, col_neg);
    plot3(1e3*p1(1), 1e3*p1(2), z0, 'k.', 'MarkerSize', 18);
    title('$C_1$ (3D bi-vectors)','Interpreter','latex');
    xlabel('$x$ (mm)','Interpreter','latex');
    ylabel('$y$ (mm)','Interpreter','latex');
    zlabel('$z$ (mm)','Interpreter','latex');
    axis equal;
    view(35,25);

    nexttile; hold on; box on; grid on;
    plot3(1e3*C2(:,1), 1e3*C2(:,2), z0*ones(size(C2,1),1), 'k', 'LineWidth', 1.0);
    plot3(1e3*C1(:,1), 1e3*C1(:,2), z0*ones(size(C1,1),1), '--', 'LineWidth', 1.6, 'Color', [0 0.7 0]);
    local_plot_bivectors3d(C2, E2, z0, qscale3, col_pos, col_neg);
    if ~isempty(picked)
        plot3(1e3*picked(:,1), 1e3*picked(:,2), z0*ones(size(picked,1),1), 'k.', 'MarkerSize', 18);
    else
        plot3(1e3*p1(1), 1e3*p1(2), z0, 'k.', 'MarkerSize', 18);
    end
    title('$C_2$ (3D bi-vectors)','Interpreter','latex');
    xlabel('$x$ (mm)','Interpreter','latex');
    ylabel('$y$ (mm)','Interpreter','latex');
    zlabel('$z$ (mm)','Interpreter','latex');
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
lim_p  = [min(abs(pU(:))), max(abs(pU(:)))];

%% ============================================================
% Fig0: |v| + quiver(Re{v_x,v_y})  AND  angle(vz)/pi
%% ============================================================
do_quiver = true;

figure('Name', sprintf('|v| & angle(vz) @ z=%.3f mm', z0mm), ...
    'position',[100 100 1400 650], 'Color','w');
set(gcf,'Renderer','opengl');

ax1 = subplot(1,2,1);
surf(Xmm, Ymm, zeros(size(X)), Vmag, 'EdgeColor','none'); view(2);
shading interp; colormap(ax1, MyColor('vik'));
if lim_v(2) > lim_v(1), clim(lim_v); end
clb = colorbar; clb.Title.Interpreter='latex'; clb.Title.String='$|v|$';
set(clb,'Fontsize',18);

axis equal;
xlim([min(Xmm(:)) max(Xmm(:))]); ylim([min(Ymm(:)) max(Ymm(:))]);
set(gca,'linewidth',2,'TickLabelInterpreter','latex'); fontsize(gca,22,'points');
xlabel('$x$ (mm)','Interpreter','latex'); ylabel('$y$ (mm)','Interpreter','latex');
title('$|v|$', 'Interpreter','latex','Fontsize',18);

if do_quiver
    hold on;
    qstep_r = max(1, round(size(X,1)/25));
    qstep_t = max(1, round(size(X,2)/60));
    Xq  = Xmm(1:qstep_r:end, 1:qstep_t:end);
    Yq  = Ymm(1:qstep_r:end, 1:qstep_t:end);
    vXq = real(vX(1:qstep_r:end, 1:qstep_t:end));
    vYq = real(vY(1:qstep_r:end, 1:qstep_t:end));
    sc  = max(hypot(vXq(:), vYq(:)));
    if sc > 0
        vXq = vXq/sc;
        vYq = vYq/sc;
    end
    quiver(Xq, Yq, vXq, vYq, 0.6, 'k', 'LineWidth', 1.0);
    hold off;
end

ax2 = subplot(1,2,2);
surf(Xmm, Ymm, zeros(size(X)), angle(vZ)/pi, 'EdgeColor','none'); view(2);
shading flat; colormap(ax2, hsv); clim(lim_ph);
clb = colorbar; clb.Title.Interpreter='latex'; clb.Title.String='$\angle v_z/\pi$';
set(clb,'Fontsize',18);

axis equal;
xlim([min(Xmm(:)) max(Xmm(:))]); ylim([min(Ymm(:)) max(Ymm(:))]);
set(gca,'linewidth',2,'TickLabelInterpreter','latex'); fontsize(gca,22,'points');
xlabel('$x$ (mm)','Interpreter','latex'); ylabel('$y$ (mm)','Interpreter','latex');
title('$\angle v_z/\pi$', 'Interpreter','latex','Fontsize',18);

sgtitle(sprintf('Spiral source, DIM @ $z=%.3f$ mm, $f_2=%.1f$ kHz, $l=%d$', ...
    z0mm, source.f2/1e3, l), 'Interpreter','latex','Fontsize',18);

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

local_save_fig_png(fig_source, sprintf('SourcePattern_%s', char(string(source.mode))));

%% ============================================================
% Pressure field (paper Fig.6-like): |p| and angle(p)/pi
%% ============================================================
local_plot_pressure_mag_phase_cart(Xp, Yp, pU, z_use, source.f2, lim_p, lim_ph);
local_save_fig_png(gcf, sprintf('DIM_p_magphase_spiral_z%.4f_l%d', z_use, l));
local_save_fig_fig(gcf, sprintf('DIM_p_magphase_spiral_z%.4f_l%d', z_use, l));

%% ============================================================
% Fig-like: oblique 3D view of velocity-direction distribution
% - arrows lie on the observation plane z=z_use
% - arrow direction: in-plane components of normalized instantaneous velocity
% - arrow color    : n_z
%% ============================================================
do_plane_dirvec = true;

if do_plane_dirvec
    % ---------------- user params ----------------
    phi_t_dir   = 0;        % instantaneous phase
    qstep_dir   = 1;        % denser arrows
    qscale_dir  = 0.42;     % arrow length
    head_dir    = 0.95;     % arrow head size
    linew_dir   = 0.9;      % arrow line width

    % ---------------- user params ----------------
    phi_t_dir   = 0;        % instantaneous phase
    qstep_dir   = 3;        % 抽稀一点：1最密，2/3更容易看方向，也更省fig
    qscale_dir  = 1.15;     % 箭头整体放大（原来0.42太短）
    head_dir    = 1.45;     % 箭头头部更大
    linew_dir   = 1.1;     % 箭头线更粗

    use_mask_r  = false;    % circular crop if needed
    mask_r_coef = 0.95;

    zplane_mm   = 1e3 * z_use;
    z_lift_mm   = 0.00;     % arrow base lift above plane, mm
    view_az_el  = [42, 28]; % oblique view, like paper style

    % ---------------- instantaneous real velocity ----------------
    ux = real(vX .* exp(-1j*phi_t_dir));
    uy = real(vY .* exp(-1j*phi_t_dir));
    uz = real(vZ .* exp(-1j*phi_t_dir));

    umag = sqrt(ux.^2 + uy.^2 + uz.^2);
    mask_valid = umag > 1e-12;

    nx = nan(size(ux));
    ny = nan(size(uy));
    nz = nan(size(uz));

    nx(mask_valid) = ux(mask_valid) ./ umag(mask_valid);
    ny(mask_valid) = uy(mask_valid) ./ umag(mask_valid);
    nz(mask_valid) = uz(mask_valid) ./ umag(mask_valid);

    if use_mask_r
        rr = hypot(X, Y);
        mask_valid = mask_valid & (rr <= mask_r_coef * max(rr(:)));
        nx(~mask_valid) = nan;
        ny(~mask_valid) = nan;
        nz(~mask_valid) = nan;
    end

    % ---------------- downsample ----------------
    Xq  = X(1:qstep_dir:end, 1:qstep_dir:end);
    Yq  = Y(1:qstep_dir:end, 1:qstep_dir:end);
    nxq = nx(1:qstep_dir:end, 1:qstep_dir:end);
    nyq = ny(1:qstep_dir:end, 1:qstep_dir:end);
    nzq = nz(1:qstep_dir:end, 1:qstep_dir:end);

    ok = isfinite(Xq) & isfinite(Yq) & isfinite(nxq) & isfinite(nyq) & isfinite(nzq);

    xq = Xq(ok);
    yq = Yq(ok);
    uq = nxq(ok);
    vq = nyq(ok);
    cq = nzq(ok);

    xq_mm = 1e3 * xq;
    yq_mm = 1e3 * yq;
    zq_mm = (zplane_mm + z_lift_mm) * ones(size(xq_mm));

    % ---------------- colormap from n_z only ----------------
    cmap_dir = turbo(256);
    cidx = round(1 + (size(cmap_dir,1)-1) * (cq + 1)/2);
    cidx = max(1, min(size(cmap_dir,1), cidx));

    % ---------------- figure ----------------
    figure('Color','w', 'Position',[120 120 1200 820], ...
        'Name', sprintf('Oblique velocity-direction map @ z=%.3f mm', zplane_mm));
    hold on; box on;

    % ---- faint observation plane ----
    surf(1e3*X, 1e3*Y, zplane_mm*ones(size(X)), ...
        zeros(size(X)), ...
        'EdgeColor','none', ...
        'FaceColor',[0.94 0.94 0.94], ...
        'FaceAlpha',0.20);

    % ---- colored arrows (color only by n_z) ----
    for k = 1:numel(xq_mm)
        quiver3(xq_mm(k), yq_mm(k), zq_mm(k), ...
            uq(k), vq(k), 0, qscale_dir, ...
            'Color', cmap_dir(cidx(k),:), ...
            'LineWidth', linew_dir, ...
            'MaxHeadSize', head_dir);
    end

    % ---- axes / view ----
    axis equal;
    xlim(1e3*[min(X(:)) max(X(:))]);
    ylim(1e3*[min(Y(:)) max(Y(:))]);

    xr = 1e3 * range(X(:));
    yr = 1e3 * range(Y(:));
    zr = max(xr, yr);
    zlim(zplane_mm + [-0.18 0.18] * zr);

    view(view_az_el(1), view_az_el(2));

    set(gca,'LineWidth',1.5,'TickLabelInterpreter','latex');
    fontsize(gca,20,'points');
    xlabel('$x$ (mm)','Interpreter','latex');
    ylabel('$y$ (mm)','Interpreter','latex');
    zlabel('$z$ (mm)','Interpreter','latex');

    title(sprintf('Normalized velocity-direction field at $z=%.3f$ mm', zplane_mm), ...
        'Interpreter','latex','FontSize',18);

    % ---- colorbar for n_z ----
    colormap(cmap_dir);
    clim([-1 1]);
    cb = colorbar;
    cb.Ticks = [-1 1];
    cb.TickLabels = {'-max','max'};
    cb.Label.String = '$n_z$';
    cb.Label.Interpreter = 'latex';
    cb.TickLabelInterpreter = 'latex';

    grid on;

    local_save_fig_png(gcf, sprintf('DIM_plane_dirvec_oblique_z%.4f_l%d', z_use, l));
    local_save_fig_fig(gcf, sprintf('DIM_plane_dirvec_oblique_z%.4f_l%d', z_use, l));
end


%% ============================================================
% Fig.8(a)(b)-like THEORETICAL skyrmion plots from ideal 3-plane-wave model
% Paper style:
%   - fixed t and z
%   - use Eq.(2)-type 3-wave interference, N=3, Phi_j=0
%   - hue  -> longitude
%   - light/dark -> latitude
%   - horizontal rectangular cells
%% ============================================================
%% ============================================================
% Fig.8(a)(b)-like THEORETICAL skyrmion plots from ideal 3-plane-wave model
%% ============================================================
do_fig8_theory = true;

if do_fig8_theory
    k0 = real(2*pi*source.f2/medium.c0);
    theta_pw = pi/4;

    lambda_pw = medium.c0 / source.f2;
    L_half = 1.5 * lambda_pw;
    dx_pw  = lambda_pw / 80;

    x_pw = -L_half:dx_pw:L_half;
    y_pw = -L_half:dx_pw:L_half;
    [Xpw, Ypw] = meshgrid(x_pw, y_pw);
    Zpw = zeros(size(Xpw));

    phi_list = [0, 2*pi/3, 4*pi/3];

    x_shift = 0.18 * lambda_pw;
    y_shift = -0.10 * lambda_pw;

    Phi_list = zeros(1,3);
    for j = 1:3
        phi_j = phi_list(j);
        kxj = k0 * sin(theta_pw) * cos(phi_j);
        kyj = k0 * sin(theta_pw) * sin(phi_j);
        Phi_list(j) = -(kxj*x_shift + kyj*y_shift);
    end

    vx_pw = zeros(size(Xpw));
    vy_pw = zeros(size(Xpw));
    vz_pw = zeros(size(Xpw));

    for j = 1:3
        phi_j = phi_list(j);
        Phi_j = Phi_list(j);

        kxj = k0 * sin(theta_pw) * cos(phi_j);
        kyj = k0 * sin(theta_pw) * sin(phi_j);
        kzj = k0 * cos(theta_pw);

        phase_j = kxj*Xpw + kyj*Ypw + kzj*Zpw + Phi_j;
        ej = exp(1j * phase_j);

        exj = sin(theta_pw) * cos(phi_j);
        eyj = sin(theta_pw) * sin(phi_j);
        ezj = cos(theta_pw);

        vx_pw = vx_pw + exj * ej;
        vy_pw = vy_pw + eyj * ej;
        vz_pw = vz_pw + ezj * ej;
    end

    ux_pw = real(vx_pw);
    uy_pw = real(vy_pw);
    uz_pw = real(vz_pw);

    umag_pw = sqrt(ux_pw.^2 + uy_pw.^2 + uz_pw.^2);
    mask_pw = umag_pw > 1e-12;

    nx_pw = nan(size(ux_pw));
    ny_pw = nan(size(uy_pw));
    nz_pw = nan(size(uz_pw));

    nx_pw(mask_pw) = ux_pw(mask_pw) ./ umag_pw(mask_pw);
    ny_pw(mask_pw) = uy_pw(mask_pw) ./ umag_pw(mask_pw);
    nz_pw(mask_pw) = uz_pw(mask_pw) ./ umag_pw(mask_pw);

    az_pw  = atan2(ny_pw, nx_pw);
    lat_pw = asin(max(-1, min(1, nz_pw)));
    RGB_pw = local_dirfield_to_rgb_paper(az_pw, lat_pw);

    % optional horizontal guide lines (only for visual reference)
    row_h = 0.90 * lambda_pw;
    y_edges_pw = (-2:2) * row_h;
    y_edges_pw = y_edges_pw(y_edges_pw >= -L_half & y_edges_pw <= L_half);

    % ---- white box / selected cell: fully continuous control in x and y ----
    % adjusted as requested:
    %   - slightly left
    %   - noticeably lower
    %   - overall larger
    xc_pw     = 0.18 * lambda_pw;   % slightly more right
    yc_pw     = -0.08 * lambda_pw;  % slightly more up
    cell_w_pw = 1.55 * lambda_pw;   % wider
    cell_h_pw = 1.55 * lambda_pw;   % taller

    x1c = xc_pw - 0.5 * cell_w_pw;
    x2c = xc_pw + 0.5 * cell_w_pw;
    y1c = yc_pw - 0.5 * cell_h_pw;
    y2c = yc_pw + 0.5 * cell_h_pw;

    mask_cell_pw = (Xpw >= x1c) & (Xpw <= x2c) & ...
        (Ypw >= y1c) & (Ypw <= y2c) & mask_pw;

    fig8_th = figure('Color','w','Position',[120 120 1320 620], ...
        'Name', 'Fig8-like THEORY from ideal 3-plane-wave model');
    tl = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    % ---------------- (a) direction map ----------------
    axa = nexttile(tl,1);
    image(axa, ...
        [min(Xpw(:))/lambda_pw, max(Xpw(:))/lambda_pw], ...
        [min(Ypw(:))/lambda_pw, max(Ypw(:))/lambda_pw], ...
        RGB_pw);
    set(axa,'YDir','normal');
    hold(axa,'on');
    box(axa,'on');

    for iy = 2:numel(y_edges_pw)-1
        plot(axa, [min(Xpw(:)) max(Xpw(:))]/lambda_pw, ...
            [y_edges_pw(iy) y_edges_pw(iy)]/lambda_pw, ...
            'w-', 'LineWidth', 1.0);
    end

    rectangle(axa, 'Position', ...
        [x1c/lambda_pw, y1c/lambda_pw, (x2c-x1c)/lambda_pw, (y2c-y1c)/lambda_pw], ...
        'EdgeColor','w', 'LineWidth', 1.8);

    axis(axa,'equal');
    xlim(axa, [min(Xpw(:)) max(Xpw(:))]/lambda_pw);
    ylim(axa, [min(Ypw(:)) max(Ypw(:))]/lambda_pw);
    set(axa,'LineWidth',1.5,'TickLabelInterpreter','latex');
    fontsize(axa,20,'points');
    xlabel(axa,'$x/\lambda$','Interpreter','latex');
    ylabel(axa,'$y/\lambda$','Interpreter','latex');
    title(axa,'(a)','Interpreter','latex','FontSize',18);

    % ---------------- (b) sphere coverage ----------------
    axb = nexttile(tl,2);
    hold(axb,'on');
    box(axb,'on');

    [Xsph, Ysph, Zsph] = sphere(160);

    surf(axb, Xsph, Ysph, Zsph, ...
        'FaceColor',[0.92 0.92 0.92], ...
        'EdgeColor','none', ...
        'FaceAlpha',0.42);

    material(axb, 'dull');
    lighting(axb, 'gouraud');
    camlight(axb, 'headlight');
    camlight(axb, 'right');

    nx_cell_pw = nx_pw(mask_cell_pw);
    ny_cell_pw = ny_pw(mask_cell_pw);
    nz_cell_pw = nz_pw(mask_cell_pw);

    az_cell_pw  = atan2(ny_cell_pw, nx_cell_pw);
    lat_cell_pw = asin(max(-1, min(1, nz_cell_pw)));
    rgb_cell_pw = local_dirpts_to_rgb_paper(az_cell_pw, lat_cell_pw);

    % slightly sparser points on the sphere
    idx_pick = 1:2:numel(nx_cell_pw);

    nx_q  = nx_cell_pw(idx_pick);
    ny_q  = ny_cell_pw(idx_pick);
    nz_q  = nz_cell_pw(idx_pick);
    rgb_q = rgb_cell_pw(idx_pick,:);

    % lift points slightly above sphere to avoid z-fighting
    r_pts = 1.015;
    px = r_pts * nx_q;
    py = r_pts * ny_q;
    pz = r_pts * nz_q;

    scatter3(axb, px, py, pz, 18, rgb_q, ...
        'filled', ...
        'MarkerEdgeColor',[0.35 0.35 0.35], ...
        'LineWidth',0.25);

    qL = 1.28;
    quiver3(axb, 0,0,0, qL,0,0, 0, 'k', 'LineWidth',1.2, 'MaxHeadSize',0.18);
    quiver3(axb, 0,0,0, 0,qL,0, 0, 'k', 'LineWidth',1.2, 'MaxHeadSize',0.18);
    quiver3(axb, 0,0,0, 0,0,qL, 0, 'k', 'LineWidth',1.2, 'MaxHeadSize',0.18);

    text(axb, qL+0.05, 0, 0, '$x$', 'Interpreter','latex', 'FontSize',18);
    text(axb, 0, qL+0.05, 0, '$y$', 'Interpreter','latex', 'FontSize',18);
    text(axb, 0, 0, qL+0.05, '$z$', 'Interpreter','latex', 'FontSize',18);

    axis(axb,'equal');
    xlim(axb,[-1.35 1.35]); ylim(axb,[-1.35 1.35]); zlim(axb,[-1.35 1.35]);
    view(axb,28,18);
    set(axb,'LineWidth',1.5,'TickLabelInterpreter','latex');
    fontsize(axb,18,'points');
    axis(axb,'off');
    title(axb,'(b)','Interpreter','latex','FontSize',18);

    local_save_fig_png(gcf, sprintf('THEORY_Fig8ab_skyrmion_3pw_f%.1fk', source.f2/1e3));
    local_save_fig_fig(gcf, sprintf('THEORY_Fig8ab_skyrmion_3pw_f%.1fk', source.f2/1e3));
end

%% ============================================================
% Fig.8-like palette legend: color/brightness vs angular coordinates
%% ============================================================
do_fig8_palette = true;

if do_fig8_palette
    n_az  = 401;
    n_lat = 201;

    az_vec  = linspace(-pi, pi, n_az);
    lat_vec = linspace(-pi/2, pi/2, n_lat);

    [AZ, LAT] = meshgrid(az_vec, lat_vec);

    RGB_leg = local_dirfield_to_rgb_paper(AZ, LAT);

    figure('Color','w', 'Position',[140 140 900 420], ...
        'Name','Fig8-like direction-color legend');

    image(az_vec/pi, lat_vec/pi, RGB_leg);
    set(gca,'YDir','normal');
    axis tight;

    set(gca,'LineWidth',1.5,'TickLabelInterpreter','latex');
    fontsize(gca,18,'points');

    xlabel('azimuth / $\pi$', 'Interpreter','latex');
    ylabel('latitude / $\pi$', 'Interpreter','latex');

    xticks([-1 -0.5 0 0.5 1]);
    yticks([-0.5 -0.25 0 0.25 0.5]);

    title('Direction-color legend', 'Interpreter','latex', 'FontSize',18);

    local_save_fig_png(gcf, 'THEORY_Fig8_palette_legend');
    local_save_fig_fig(gcf, 'THEORY_Fig8_palette_legend');
end


%% ==================== local functions ====================

function local_plot_mag_phase_cart(X, Y, V, nameStr, z_use, f_hz, lim_mag, lim_ph, isS)
if nargin < 9, isS = false; end

Xmm = 1e3 * X;
Ymm = 1e3 * Y;
zmm = 1e3 * z_use;

figure('Color','w', 'Position',[100 100 1400 650], ...
    'Name',sprintf('%s (DIM) @ z=%.3f mm', nameStr, zmm));
set(gcf,'Renderer','opengl');

ax1 = subplot(1,2,1);
surf(Xmm, Ymm, zeros(size(X)), abs(V), 'EdgeColor','none'); view(2);
shading interp; colormap(ax1, MyColor('vik'));
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
xlim([min(Xmm(:)) max(Xmm(:))]); ylim([min(Ymm(:)) max(Ymm(:))]);
set(gca,'linewidth',2,'TickLabelInterpreter','latex'); fontsize(gca,22,'points');
xlabel('$x$ (mm)','Interpreter','latex'); ylabel('$y$ (mm)','Interpreter','latex');
title(sprintf('Magnitude @ $z=%.3f$ mm', zmm), 'Interpreter','latex','Fontsize',18);

ax2 = subplot(1,2,2);
surf(Xmm, Ymm, zeros(size(X)), angle(V)/pi, 'EdgeColor','none'); view(2);
shading flat; colormap(ax2, hsv); clim(lim_ph);
clb = colorbar;
clb.Title.Interpreter = 'latex';
if ~isS
    clb.Title.String = sprintf('$\\angle %s/\\pi$', strrep(nameStr,'_','\_'));
else
    clb.Title.String = '$\angle s/\pi$';
end
set(clb,'Fontsize',18);

axis equal;
xlim([min(Xmm(:)) max(Xmm(:))]); ylim([min(Ymm(:)) max(Ymm(:))]);
set(gca,'linewidth',2,'TickLabelInterpreter','latex'); fontsize(gca,22,'points');
xlabel('$x$ (mm)','Interpreter','latex'); ylabel('$y$ (mm)','Interpreter','latex');
title('Phase', 'Interpreter','latex','Fontsize',18);

sgtitle(sprintf('DIM, $f=%.1f$ kHz, $z=%.3f$ mm', f_hz/1e3, zmm), ...
    'Interpreter','latex','Fontsize',18);
end

function vn = local_spiral_binary_grating_vn(X,Y,v0,a,r0,P,dP,M,handed,l,mode)
if nargin < 10 || isempty(mode)
    mode = 'spiral_abnormal';
end
mode = lower(string(mode));

rho = hypot(X,Y);
phi = atan2(Y,X);

vn = zeros(size(rho));

r1 = r0 + M*P;

mask_ap = (rho <= a) & (rho >= r0) & (rho <= r1);
if ~any(mask_ap(:))
    return;
end

b = P/(2*pi);

switch mode
    case "spiral_abnormal"
        rho_m = rho(mask_ap);
        phi_m = phi(mask_ap);

        phi_eff = l * phi_m;

        u = rho_m - r0 - b*(handed*phi_eff);
        t = mod(u, P);

        on = (t < dP);
        g = ones(size(t));
        g(~on) = -1;

        vn(mask_ap) = v0 * g;

    case "spiral_normal"
        xq = X(mask_ap);
        yq = Y(mask_ap);

        psi_max = 2*pi*M;
        npsi = max(4000, ceil(psi_max / (2*pi) * 1200));
        psi = linspace(0, psi_max, npsi);

        rho_c = r0 + b*psi;
        theta_c = psi / (handed * l);

        xc = rho_c .* cos(theta_c);
        yc = rho_c .* sin(theta_c);

        keep = (xc.^2 + yc.^2) <= a^2;
        xc = xc(keep);
        yc = yc(keep);

        dn_min2 = inf(size(xq));

        blk = 4000;
        for ib = 1:blk:numel(xq)
            ie = min(ib + blk - 1, numel(xq));

            Xb = xq(ib:ie);
            Yb = yq(ib:ie);

            dx = Xb(:) - xc(:).';
            dy = Yb(:) - yc(:).';
            d2 = dx.^2 + dy.^2;

            dn_min2(ib:ie) = min(d2, [], 2);
        end

        dn_min = sqrt(dn_min2);
        on = dn_min <= (dP/2);

        val = -v0 * ones(size(xq));
        val(on) = +v0;

        vn(mask_ap) = val;

    case "triangle_normal"
        vn = zeros(size(X));

        pts = local_make_triangle_spiral_centerline(r0, M, P, handed);
        if size(pts,1) < 6
            return;
        end

        segP0 = pts(1:end-1,:);
        segP1 = pts(2:end,:);

        seg0_outer = segP0(1:3,:);
        seg1_outer = segP1(1:3,:);

        ref_outer = mean(pts(1:3,:), 1);
        outer_poly = local_triangle_from_three_segments(seg0_outer, seg1_outer, -dP/2, ref_outer);

        seg0_inner = segP0(end-2:end,:);
        seg1_inner = segP1(end-2:end,:);

        ref_inner = mean(pts(end-3:end,:), 1);
        inner_poly = local_triangle_from_three_segments(seg0_inner, seg1_inner, +dP/2, ref_inner);

        [in_outer, on_outer] = inpolygon(X, Y, outer_poly(:,1), outer_poly(:,2));
        inside_outer = in_outer | on_outer;

        [in_inner, on_inner_poly] = inpolygon(X, Y, inner_poly(:,1), inner_poly(:,2));
        inside_inner = in_inner | on_inner_poly;

        mask_ap = (rho <= a) & inside_outer;
        if ~any(mask_ap(:))
            return;
        end

        xq = X(mask_ap);
        yq = Y(mask_ap);

        inside_inner_q = inside_inner(mask_ap);

        d2_min = inf(size(xq));

        nSeg  = size(segP0,1);
        blk_q = 2000;
        blk_s = 100;

        for iq = 1:blk_q:numel(xq)
            jq = min(iq + blk_q - 1, numel(xq));

            Xb = xq(iq:jq);
            Yb = yq(iq:jq);

            d2_blk = inf(numel(Xb),1);

            for is = 1:blk_s:nSeg
                js = min(is + blk_s - 1, nSeg);

                P0 = segP0(is:js,:);
                P1 = segP1(is:js,:);

                d2_now = local_point_to_segments_min_d2(Xb, Yb, P0, P1);
                d2_blk = min(d2_blk, d2_now);
            end

            d2_min(iq:jq) = d2_blk;
        end

        on_strip = sqrt(d2_min) <= (dP/2);

        val = -v0 * ones(size(xq));
        val(on_strip) = +v0;
        val(inside_inner_q) = 0;

        vn(mask_ap) = val;

    otherwise
        error('Unknown source.mode: %s', char(mode));
end
end

function poly = local_triangle_from_three_segments(segP0, segP1, offset_dist, ref_point)
nvec = zeros(3,2);
cval = zeros(3,1);

for k = 1:3
    p0 = segP0(k,:);
    p1 = segP1(k,:);

    e = p1 - p0;
    Le = norm(e);
    if Le < eps
        error('Degenerate segment encountered in local_triangle_from_three_segments.');
    end

    n = [-e(2), e(1)] / Le;

    if dot(n, ref_point - p0) < 0
        n = -n;
    end

    nvec(k,:) = n;
    cval(k)   = dot(n, p0) + offset_dist;
end

poly = zeros(3,2);
for k = 1:3
    k2 = mod(k,3) + 1;
    poly(k,:) = local_intersect_two_lines(nvec(k,:), cval(k), nvec(k2,:), cval(k2));
end
end

function p = local_intersect_two_lines(n1, c1, n2, c2)
A = [n1(:).'; n2(:).'];
b = [c1; c2];

if abs(det(A)) < 1e-12
    error('Lines are nearly parallel in local_intersect_two_lines.');
end

p = (A \ b).';
end

function pts = local_make_triangle_spiral_centerline(r0, M, P, handed)
if nargin < 4 || isempty(handed)
    handed = -1;
end

r_outer = r0 + 2*M*P;
L0 = sqrt(3)*r_outer;
dL = 2*sqrt(3)*P/3;
Nseg = 3*M;

if handed < 0
    angles_deg = [0, 120, -120];
else
    angles_deg = [0, -120, 120];
end

pts = zeros(Nseg+1, 2);
pts(1,:) = [0, 0];

nkeep = 1;
for k = 1:Nseg
    Lk = L0 - (k-1)*dL;
    if Lk <= 0
        break;
    end

    ang = angles_deg(mod(k-1,3)+1);
    dir = [cosd(ang), sind(ang)];
    pts(k+1,:) = pts(k,:) + Lk * dir;
    nkeep = k+1;
end

pts = pts(1:nkeep,:);

if size(pts,1) >= 3
    pc = mean(pts(1:3,:), 1);
    pts = pts - pc;
end
end

function d2_min = local_point_to_segments_min_d2(Xq, Yq, P0, P1)
Xq = Xq(:);
Yq = Yq(:);

x0 = P0(:,1).';
y0 = P0(:,2).';
x1 = P1(:,1).';
y1 = P1(:,2).';

vx = x1 - x0;
vy = y1 - y0;
vv = vx.^2 + vy.^2;
vv(vv < eps) = eps;

wx = Xq - x0;
wy = Yq - y0;

t = (wx.*vx + wy.*vy) ./ vv;
t = max(0, min(1, t));

xp = x0 + t.*vx;
yp = y0 + t.*vy;

d2 = (Xq - xp).^2 + (Yq - yp).^2;
d2_min = min(d2, [], 2);
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
    savefig(fig_handle, fp);
catch
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

function local_write_runinfo_txt(save_dir, medium, source, calc, fig, extra, copt)
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
fprintf(fid, 'source.profile=''%s'';\n', char(string(source.profile)));
if isfield(source,'mode')
    fprintf(fid, 'source.mode=''%s'';\n', char(string(source.mode)));
end
fprintf(fid, 'source.a=%.15g; source.v0=%.15g; source.v_ratio=%.15g;\n', ...
    source.a, source.v0, source.v_ratio);
if isfield(source,'m_custom')
    fprintf(fid, 'source.m_custom=%d;\n', source.m_custom);
end
fprintf(fid, 'source.f1=%.15g; source.f2=%.15g; source.fa=%.15g;\n\n', ...
    source.f1, source.f2, source.fa);

fprintf(fid, '%% calc.fht\n');
cf = calc.fht;
fprintf(fid, 'N_FHT=%d; rho_max=%.15g; Nh_scale=%.15g; NH_scale=%.15g; Nh_v_scale=%.15g; zu_max=%.15g; za_max=%.15g;\n\n', ...
    cf.N_FHT, cf.rho_max, cf.Nh_scale, cf.NH_scale, cf.Nh_v_scale, cf.zu_max, cf.za_max);

fprintf(fid, '%% calc.dim\n');
cd = calc.dim;
fprintf(fid, 'use_freq=''%s''; dis_coe=%.15g; margin=%.15g; src_discretization=''%s'';\n', ...
    char(string(cd.use_freq)), cd.dis_coe, cd.margin, char(string(cd.src_discretization)));
fprintf(fid, 'dim_method=''%s''; use_parallel=%s; num_workers=%d;\n\n', ...
    char(string(cd.method)), local_bool_str(cd.use_parallel), cd.num_workers);

fprintf(fid, '\n%% fig\n');
fprintf(fid, 'fig.unwrap=%s;\n', local_bool_str(fig.unwrap));

fprintf(fid, '\n%% extra (for reproduction)\n');
fprintf(fid, 'l=%d;\n', extra.l);
fprintf(fid, 'r_boundary_coef=%.15g;\n', extra.r_boundary_coef);
fprintf(fid, 'r_small_k=%d;  %% r_small = r_small_k * dxy\n', extra.r_small_k);
fprintf(fid, 'r_large_k=%d;  %% r_large = r_large_k * dxy\n', extra.r_large_k);

if isfield(extra,'use_custom_obs_range') && ~isempty(extra.use_custom_obs_range)
    fprintf(fid, 'use_custom_obs_range=%s;\n', local_bool_str(extra.use_custom_obs_range));
end
if isfield(extra,'z_use') && ~isempty(extra.z_use)
    fprintf(fid, 'z_use=%.15g;\n', extra.z_use);
end

if isfield(extra,'xlim_v') && numel(extra.xlim_v)==2
    fprintf(fid, 'xlim_v=[%.15g, %.15g];\n', extra.xlim_v(1), extra.xlim_v(2));
end
if isfield(extra,'ylim_v') && numel(extra.ylim_v)==2
    fprintf(fid, 'ylim_v=[%.15g, %.15g];\n', extra.ylim_v(1), extra.ylim_v(2));
end
if isfield(extra,'xlim_p') && numel(extra.xlim_p)==2
    fprintf(fid, 'xlim_p=[%.15g, %.15g];\n', extra.xlim_p(1), extra.xlim_p(2));
end
if isfield(extra,'ylim_p') && numel(extra.ylim_p)==2
    fprintf(fid, 'ylim_p=[%.15g, %.15g];\n', extra.ylim_p(1), extra.ylim_p(2));
end

if isfield(extra,'dx_obs_v') && ~isempty(extra.dx_obs_v)
    fprintf(fid, 'dx_obs_v=%.15g;\n', extra.dx_obs_v);
end
if isfield(extra,'dy_obs_v') && ~isempty(extra.dy_obs_v)
    fprintf(fid, 'dy_obs_v=%.15g;\n', extra.dy_obs_v);
end
if isfield(extra,'dx_obs_p') && ~isempty(extra.dx_obs_p)
    fprintf(fid, 'dx_obs_p=%.15g;\n', extra.dx_obs_p);
end
if isfield(extra,'dy_obs_p') && ~isempty(extra.dy_obs_p)
    fprintf(fid, 'dy_obs_p=%.15g;\n', extra.dy_obs_p);
end

if isfield(extra,'Nx_v') && ~isempty(extra.Nx_v)
    fprintf(fid, 'Nx_v=%d;\n', extra.Nx_v);
end
if isfield(extra,'Ny_v') && ~isempty(extra.Ny_v)
    fprintf(fid, 'Ny_v=%d;\n', extra.Ny_v);
end
if isfield(extra,'Nx_p') && ~isempty(extra.Nx_p)
    fprintf(fid, 'Nx_p=%d;\n', extra.Nx_p);
end
if isfield(extra,'Ny_p') && ~isempty(extra.Ny_p)
    fprintf(fid, 'Ny_p=%d;\n', extra.Ny_p);
end

fprintf(fid, '\n%% copt (for finding q)\n');
fprintf(fid, 'max_candidates=%d;\n', copt.max_candidates);
fprintf(fid, 'smag_rel_th=%.15g;\n', copt.smag_rel_th);
fprintf(fid, 'min_sep_pix=%.15g;\n', copt.min_sep_pix);
fprintf(fid, 'ring_radius=%.15g;\n', copt.ring_radius);

fprintf(fid, '\n===== END =====\n');
end

function s = local_bool_str(x)
if logical(x)
    s = 'true';
else
    s = 'false';
end
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
tt = linspace(0,2*pi,Nsamp+1);
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
C_use = C;
if size(C,1) >= 2 && norm(C(end,:) - C(1,:)) < 1e-12
    C_use = C(1:end-1,:);
end

NsC = size(C_use,1);
NsE = size(E,1);
Ns  = min(NsC, NsE);

x  = 1e3 * C_use(1:Ns,1);
y  = 1e3 * C_use(1:Ns,2);
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
sz = size(xq);
xq = xq(:);
yq = yq(:);

is_meshgrid = ...
    size(X,1) > 1 && size(X,2) > 1 && ...
    all(abs(diff(X(:,1))) < 1e-12) && ...
    all(abs(diff(Y(1,:))) < 1e-12);

if is_meshgrid
    Xn  = X.';
    Yn  = Y.';
    phn = phS.';
    Fph = griddedInterpolant(Xn, Yn, phn, 'linear', 'none');
else
    Fph = griddedInterpolant(X, Y, phS, 'linear', 'none');
end

ph = Fph(xq, yq);
ph = reshape(ph, sz);
end

function A = local_norm_rows(A)
nn = sqrt(sum(A.^2,2));
A = A ./ max(nn, eps);
end

function local_plot_pressure_mag_phase_cart(X, Y, P, z_use, f_hz, lim_mag, lim_ph)
Xmm = 1e3 * X;
Ymm = 1e3 * Y;
zmm = 1e3 * z_use;

figure('Color','w', 'Position',[100 100 1400 650], ...
    'Name', sprintf('Pressure (DIM) @ z=%.3f mm', zmm));
set(gcf,'Renderer','opengl');

ax1 = subplot(1,2,1);
surf(Xmm, Ymm, zeros(size(X)), abs(P), 'EdgeColor','none'); view(2);
shading interp; colormap(ax1, MyColor('vik'));
if all(isfinite(lim_mag)) && lim_mag(2) > lim_mag(1)
    clim(lim_mag);
end
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$|p|$';
set(clb,'Fontsize',18);

axis equal;
xlim([min(Xmm(:)) max(Xmm(:))]);
ylim([min(Ymm(:)) max(Ymm(:))]);
set(gca,'linewidth',2,'TickLabelInterpreter','latex');
fontsize(gca,22,'points');
xlabel('$x$ (mm)','Interpreter','latex');
ylabel('$y$ (mm)','Interpreter','latex');
title(sprintf('Magnitude @ $z=%.3f$ mm', zmm), ...
    'Interpreter','latex','Fontsize',18);

ax2 = subplot(1,2,2);
surf(Xmm, Ymm, zeros(size(X)), angle(P)/pi, 'EdgeColor','none'); view(2);
shading flat; colormap(ax2, hsv); clim(lim_ph);
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$\angle p/\pi$';
set(clb,'Fontsize',18);

axis equal;
xlim([min(Xmm(:)) max(Xmm(:))]);
ylim([min(Ymm(:)) max(Ymm(:))]);
set(gca,'linewidth',2,'TickLabelInterpreter','latex');
fontsize(gca,22,'points');
xlabel('$x$ (mm)','Interpreter','latex');
ylabel('$y$ (mm)','Interpreter','latex');
title('Phase', 'Interpreter','latex','Fontsize',18);

sgtitle(sprintf('Ultrasonic pressure (DIM), $f=%.1f$ kHz, $z=%.3f$ mm', ...
    f_hz/1e3, zmm), 'Interpreter','latex','Fontsize',18);
end


% ---------------------- Fig 8 ----------------------

function RGB = local_dirfield_to_rgb_paper(az, lat)
h = mod(az + pi, 2*pi) / (2*pi);

% desired behavior:
% upper half: brighter color -> white
% lower half: darker color  -> black
gamma_w    = 1.40;   % white appears mainly near the top
gamma_b    = 6.6;   % black appears mainly near the bottom
gamma_s    = 0.72;   % color band width

sat_gain   = 0.96;
rgb_eq = hsv2rgb(cat(3, h, sat_gain*ones(size(h)), ones(size(h))));
rgb_eq = 1.08 * rgb_eq;

t  = (2/pi) * lat;   % latitude normalized to [-1,1]
tp = max(t, 0);      % upper half
tn = max(-t, 0);     % lower half

% keep color over a broad region
eq = max(0, 1 - abs(t).^gamma_s);

% upper colors intrinsically brighter, lower colors intrinsically darker
bright_fac = 0.88 + 0.22 * tp;   % upper half -> brighter
dark_fac   = 0.58 + 0.42 * tn;   % lower half -> darker

% white / black weights
w_top = min(1, 1.35 * tp.^gamma_w);
b_bot = min(1, 0.82 * tn.^gamma_b);

RGB = zeros([size(h), 3]);
for k = 1:3
    ch = rgb_eq(:,:,k);

    % start from color with different upper/lower brightness
    val = eq .* ch;
    val(t >= 0) = 1.35 * val(t >= 0);
    val(t >= 0) = bright_fac(t >= 0) .* val(t >= 0);
    val(t <  0) = (1 - 0.28*dark_fac(t < 0)) .* val(t < 0);

    % blend to white at the top
    val = val .* (1 - 0.88*w_top) + w_top .* 1;

    % blend to black at the bottom
    val = val .* (1 - 0.72*b_bot);

    RGB(:,:,k) = val;
end

bad = ~isfinite(az) | ~isfinite(lat);
for k = 1:3
    tmp = RGB(:,:,k);
    tmp(bad) = 1;
    RGB(:,:,k) = tmp;
end

RGB = min(max(RGB,0),1);
end

function rgb = local_dirpts_to_rgb_paper(az, lat)
h = mod(az + pi, 2*pi) / (2*pi);

gamma_w    = 1.40;
gamma_b    = 6.6;
gamma_s    = 0.72;

sat_gain   = 0.96;
rgb_eq = hsv2rgb([h(:), sat_gain*ones(numel(h),1), ones(numel(h),1)]);
rgb_eq = 1.08 * rgb_eq;

t  = (2/pi) * lat(:);
tp = max(t, 0);
tn = max(-t, 0);

eq = max(0, 1 - abs(t).^gamma_s);

bright_fac = 0.88 + 0.22 * tp;
dark_fac   = 0.58 + 0.42 * tn;

w_top = min(1, 1.35 * tp.^gamma_w);
b_bot = min(1, 0.82 * tn.^gamma_b);

rgb = zeros(numel(h),3);
for k = 1:3
    ch = rgb_eq(:,k);

    val = eq .* ch;
    val(t >= 0) = 1.35 * val(t >= 0);
    val(t >= 0) = bright_fac(t >= 0) .* val(t >= 0);
    val(t <  0) = (1 - 0.28*dark_fac(t < 0)) .* val(t < 0);

    val = val .* (1 - 0.88*w_top) + w_top .* 1;
    val = val .* (1 - 0.72*b_bot);

    rgb(:,k) = val;
end

rgb = min(max(rgb,0),1);
end