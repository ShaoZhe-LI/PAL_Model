%% ============================================================
% MAIN: f1 ultrasound only (2D only)
% - King: compute on its own FHT grid
% - DIM : use the same FHT grid, but rho is downsampled
% - z range: 0 ~ 1 m
%
% Outputs:
%   1) xOz amplitude / phase (King vs DIM)
%   2) radial amplitude / phase at z = 1 m (King vs DIM)
%   3) xOz real / imag (King vs DIM)
%   4) radial real / imag at z = 1 m (King vs DIM)
%   5) timing summary for both King and DIM

% 测试两种方法不同网格的超声结果
clear; clc; close all;



%%
global SAVE_PNG SAVE_DIR
SAVE_PNG = false;
SAVE_DIR = '';

%% -------------------- medium --------------------
medium.c0 = 343;
medium.rho0 = 1.21;
medium.beta = 1.2;
medium.pref = 2e-5;
medium.use_absorp = true;
medium.atten_handle = @(f) AbsorpAttenCoef(f);

%% -------------------- source --------------------
source.profile = 'Vortex-m';
source.a       = 0.1;
source.v0      = 0.172;
source.v_ratio = 1;
source.m       = 3;
source.F       = 0.2;

source.f1 = 11.8e3;
source.fa = 2e3;
source.f2 = source.f1 + source.fa;

%% -------------------- calc --------------------
calc = struct();

% --- FHT / King ---
calc.fht.N_FHT      = 32768;
calc.fht.rho_max    = 1.2;
calc.fht.Nh_scale   = 1.2;
calc.fht.NH_scale   = 4;
calc.fht.Nh_v_scale = 1.1;
calc.fht.zu_max     = 1.0;
calc.fht.za_max     = 1.0;
calc.fht.delta = medium.c0 / source.f1 / 2;

% --- DIM source discretization ---
calc.dim.method             = 'rayleigh';   % 'rayleigh' or 'asm'
calc.dim.use_freq           = 'f1';
calc.dim.dis_coe            = 16;
calc.dim.margin             = 1;
calc.dim.src_discretization = 'polar';
calc.dim.use_parallel       = true;
calc.dim.num_workers        = 4;

% --- King analytic spectrum stability ---
calc.king.gspec_method       = 'analytic';
calc.king.eps_phase          = 1e-3;
calc.king.kz_min             = 1e-12;
calc.king.band_refine.enable = true;

% --- ASM settings (used only if calc.dim.method='asm') ---
calc.asm.pad_factor = 16;
calc.asm.kzz_eps    = 1e-12;

%% -------------------- display / sampling setup --------------------
z_target = 1.0;    % radial profile at z = 1 m
ds_rho   = 16;     % rho downsample factor for DIM, 32

fig.unwrap_phase_1d    = true;    % radial phase curve
fig.mask_low_amp_phase = false;   % mask low-amplitude phase
fig.phase_amp_ratio    = 1e-3;    % threshold ratio

%% ============================================================
% PREPARE FHT GRID FIRST
%% ============================================================
[~, fht_tmp] = make_source_velocity(source, medium, calc);

rho_full = fht_tmp.xh(:);         % FHT rho grid
z_full   = fht_tmp.z_ultra(:);    % FHT z grid

% crop to requested range
rho_full = rho_full(rho_full <= calc.fht.rho_max);
z_full   = z_full(z_full <= 1.0);

% nearest z to target
[~, iz_view_full] = min(abs(z_full - z_target));
z_view = z_full(iz_view_full);

% downsample rho for DIM
rho_ds = rho_full(1:ds_rho:end);
if rho_ds(end) ~= rho_full(end)
    rho_ds = [rho_ds; rho_full(end)];
end

fprintf('FHT full rho points : %d\n', numel(rho_full));
fprintf('DIM rho points      : %d (downsample factor ~ %d)\n', numel(rho_ds), ds_rho);
fprintf('z points            : %d, z in [0, %.6f] m\n', numel(z_full), z_full(end));
fprintf('z_view              : %.6f m\n', z_view);

%% ============================================================
% KING: compute f1 on its own FHT grid
%% ============================================================
fprintf('\n==== KING f1 field ====\n');
tic;
resK = calc_ultrasound_field(source, medium, calc, [], 'king');
tKing = toc;

rhoK_all = resK.king.rho(:);
zK_all   = resK.king.z(:).';

id_rK = find(rhoK_all <= calc.fht.rho_max, 1, 'last');
id_zK = find(zK_all   <= 1.0,            1, 'last');
if isempty(id_rK), id_rK = numel(rhoK_all); end
if isempty(id_zK), id_zK = numel(zK_all);   end

rhoK = rhoK_all(1:id_rK);
zK   = zK_all(1:id_zK);
pK_rz = resK.king.p_f1(1:id_rK, 1:id_zK);   % Nr x Nz

[~, izK] = min(abs(zK - z_view));
z_view_K = zK(izK);
pK_r = pK_rz(:, izK);

NKing = numel(rhoK) * numel(zK);
tKing_per_point = tKing / NKing;
tKing_per_z     = tKing / numel(zK);

%% ============================================================
% DIM: use FHT-aligned grid, but rho is downsampled
%% ============================================================
fprintf('\n==== DIM f1 field ====\n');

obs_grid = struct();
obs_grid.dim.x = rho_ds(:).';     % FHT rho grid with downsampling
obs_grid.dim.y = 0;
obs_grid.dim.z = z_full(:).';     % FHT z grid, cropped to 0~1 m
obs_grid.dim.block_size = 2e5;

% Optional alternative uniform sampling:
dxyz = calc.fht.delta / 4;
obs_grid.dim.x = 0:dxyz:calc.fht.rho_max;
obs_grid.dim.y = 0;
obs_grid.dim.z = 0:calc.fht.delta:calc.fht.zu_max;

Nx_dim = numel(obs_grid.dim.x);
Ny_dim = numel(obs_grid.dim.y);
Nz_dim = numel(obs_grid.dim.z);
N_dim_obs = Nx_dim * Ny_dim * Nz_dim;

tic;
resD = calc_ultrasound_field(source, medium, calc, obs_grid, 'dim');
tDIM = toc;

tDIM_per_point = tDIM / N_dim_obs;
tDIM_per_z     = tDIM / Nz_dim;

if strcmpi(calc.dim.method,'rayleigh')
    % Ny x Nx x Nz, and here Ny = 1
    pD_rz = squeeze(resD.dim.p_f1(1,:,:));   % Nx x Nz
    rhoD  = resD.dim.x(:);
    zD    = resD.dim.z(:).';
else
    xA = resD.dim.x(:);
    yA = resD.dim.y(:);
    zA = resD.dim.z(:).';
    [~, iy0] = min(abs(yA - 0));
    pD_rz = squeeze(resD.dim.p_f1(iy0,:,:)); % Nx x Nz
    rhoD  = xA;
    zD    = zA;
end

[~, izD] = min(abs(zD - z_view));
z_view_D = zD(izD);
pD_r = pD_rz(:, izD);

%% ============================================================
% Build xOz maps
%% ============================================================
% mirror positive rho to signed x
xK_xoz = [-flipud(rhoK); rhoK];
pK_xoz = [flipud(pK_rz); pK_rz];

xD_xoz = [-flipud(rhoD); rhoD];
pD_xoz = [flipud(pD_rz); pD_rz];

AK_xoz  = abs(pK_xoz);
AD_xoz  = abs(pD_xoz);
PhK_xoz = angle(pK_xoz);
PhD_xoz = angle(pD_xoz);

ReK_xoz = real(pK_xoz);
ImK_xoz = imag(pK_xoz);
ReD_xoz = real(pD_xoz);
ImD_xoz = imag(pD_xoz);

% optional mask low-amplitude phase
if fig.mask_low_amp_phase
    thK_xoz = max(AK_xoz(:)) * fig.phase_amp_ratio;
    thD_xoz = max(AD_xoz(:)) * fig.phase_amp_ratio;
    PhK_xoz(AK_xoz < thK_xoz) = NaN;
    PhD_xoz(AD_xoz < thD_xoz) = NaN;
end

%% ============================================================
% Radial phase for 1D plot
%% ============================================================
if fig.unwrap_phase_1d
    phK_r = unwrap(angle(pK_r));
    phD_r = unwrap(angle(pD_r));
else
    phK_r = angle(pK_r);
    phD_r = angle(pD_r);
end

if fig.mask_low_amp_phase
    thK_1d = max(abs(pK_r)) * fig.phase_amp_ratio;
    thD_1d = max(abs(pD_r)) * fig.phase_amp_ratio;
    phK_r(abs(pK_r) < thK_1d) = NaN;
    phD_r(abs(pD_r) < thD_1d) = NaN;
end

ReK_r = real(pK_r);
ImK_r = imag(pK_r);
ReD_r = real(pD_r);
ImD_r = imag(pD_r);

%% ============================================================
% Color limits
%% ============================================================
amp_lim_xoz = [0, max([AK_xoz(:); AD_xoz(:)])];
ph_lim      = [-pi, pi];

re_lim_xoz = [min([ReK_xoz(:); ReD_xoz(:)]), max([ReK_xoz(:); ReD_xoz(:)])];
im_lim_xoz = [min([ImK_xoz(:); ImD_xoz(:)]), max([ImK_xoz(:); ImD_xoz(:)])];

%% ============================================================
% Plot 1: xOz amplitude / phase
%% ============================================================
figure('Name','xOz |p| and phase: King vs DIM', 'Position',[80 60 1500 900]);

subplot(2,2,1);
pcolor(zK, xK_xoz, AK_xoz); shading interp;
colormap(gca, parula(256));
clim(amp_lim_xoz); colorbar;
xlabel('z (m)'); ylabel('x (m)');
title('King: |p_{f1}| on xOz');
xlim([0 zK(end)]); ylim([-calc.fht.rho_max calc.fht.rho_max]);
pbaspect([zK(end), 2*calc.fht.rho_max, 1]);

subplot(2,2,2);
pcolor(zK, xK_xoz, PhK_xoz); shading interp;
colormap(gca, hsv(256));
clim(ph_lim); colorbar;
xlabel('z (m)'); ylabel('x (m)');
title('King: \angle p_{f1} on xOz');
xlim([0 zK(end)]); ylim([-calc.fht.rho_max calc.fht.rho_max]);
pbaspect([zK(end), 2*calc.fht.rho_max, 1]);

subplot(2,2,3);
pcolor(zD, xD_xoz, AD_xoz); shading interp;
colormap(gca, parula(256));
clim(amp_lim_xoz); colorbar;
xlabel('z (m)'); ylabel('x (m)');
title(sprintf('DIM-%s: |p_{f1}| on xOz', upper(calc.dim.method)));
xlim([0 zD(end)]); ylim([-calc.fht.rho_max calc.fht.rho_max]);
pbaspect([zD(end), 2*calc.fht.rho_max, 1]);

subplot(2,2,4);
pcolor(zD, xD_xoz, PhD_xoz); shading interp;
colormap(gca, hsv(256));
clim(ph_lim); colorbar;
xlabel('z (m)'); ylabel('x (m)');
title(sprintf('DIM-%s: \\angle p_{f1} on xOz', upper(calc.dim.method)));
xlim([0 zD(end)]); ylim([-calc.fht.rho_max calc.fht.rho_max]);
pbaspect([zD(end), 2*calc.fht.rho_max, 1]);

%% ============================================================
% Plot 2: radial amplitude / phase at z = z_view
%% ============================================================
figure('Name','Radial profile at z=1m', 'Position',[150 120 1200 480]);

subplot(1,2,1);
plot(rhoK, abs(pK_r), 'LineWidth', 1.8); hold on;
plot(rhoD, abs(pD_r), '--', 'LineWidth', 1.8);
grid on;
xlabel('\rho (m)');
ylabel('|p_{f1}|');
title(sprintf('Radial amplitude at z = %.6f m', z_view));
legend('King','DIM','Location','best');
xlim([0 calc.fht.rho_max]);

subplot(1,2,2);
plot(rhoK, phK_r, 'LineWidth', 1.8); hold on;
plot(rhoD, phD_r, '--', 'LineWidth', 1.8);
grid on;
xlabel('\rho (m)');
ylabel('Phase (rad)');
title(sprintf('Radial phase at z = %.6f m', z_view));
legend('King','DIM','Location','best');
xlim([0 calc.fht.rho_max]);

%% ============================================================
% Plot 3: xOz real / imag
%% ============================================================
figure('Name','xOz real and imag: King vs DIM', 'Position',[100 70 1500 900]);

subplot(2,2,1);
pcolor(zK, xK_xoz, ReK_xoz); shading interp;
colormap(gca, parula(256));
clim(re_lim_xoz); colorbar;
xlabel('z (m)'); ylabel('x (m)');
title('King: Re\{p_{f1}\} on xOz');
xlim([0 zK(end)]); ylim([-calc.fht.rho_max calc.fht.rho_max]);
pbaspect([zK(end), 2*calc.fht.rho_max, 1]);

subplot(2,2,2);
pcolor(zK, xK_xoz, ImK_xoz); shading interp;
colormap(gca, parula(256));
clim(im_lim_xoz); colorbar;
xlabel('z (m)'); ylabel('x (m)');
title('King: Im\{p_{f1}\} on xOz');
xlim([0 zK(end)]); ylim([-calc.fht.rho_max calc.fht.rho_max]);
pbaspect([zK(end), 2*calc.fht.rho_max, 1]);

subplot(2,2,3);
pcolor(zD, xD_xoz, ReD_xoz); shading interp;
colormap(gca, parula(256));
clim(re_lim_xoz); colorbar;
xlabel('z (m)'); ylabel('x (m)');
title(sprintf('DIM-%s: Re\\{p_{f1}\\} on xOz', upper(calc.dim.method)));
xlim([0 zD(end)]); ylim([-calc.fht.rho_max calc.fht.rho_max]);
pbaspect([zD(end), 2*calc.fht.rho_max, 1]);

subplot(2,2,4);
pcolor(zD, xD_xoz, ImD_xoz); shading interp;
colormap(gca, parula(256));
clim(im_lim_xoz); colorbar;
xlabel('z (m)'); ylabel('x (m)');
title(sprintf('DIM-%s: Im\\{p_{f1}\\} on xOz', upper(calc.dim.method)));
xlim([0 zD(end)]); ylim([-calc.fht.rho_max calc.fht.rho_max]);
pbaspect([zD(end), 2*calc.fht.rho_max, 1]);

%% ============================================================
% Plot 4: radial real / imag at z = z_view
%% ============================================================
figure('Name','Radial real and imag at z=1m', 'Position',[170 130 1200 480]);

subplot(1,2,1);
plot(rhoK, ReK_r, 'LineWidth', 1.8); hold on;
plot(rhoD, ReD_r, '--', 'LineWidth', 1.8);
grid on;
xlabel('\rho (m)');
ylabel('Re\{p_{f1}\}');
title(sprintf('Radial real part at z = %.6f m', z_view));
legend('King','DIM','Location','best');
xlim([0 calc.fht.rho_max]);

subplot(1,2,2);
plot(rhoK, ImK_r, 'LineWidth', 1.8); hold on;
plot(rhoD, ImD_r, '--', 'LineWidth', 1.8);
grid on;
xlabel('\rho (m)');
ylabel('Im\{p_{f1}\}');
title(sprintf('Radial imaginary part at z = %.6f m', z_view));
legend('King','DIM','Location','best');
xlim([0 calc.fht.rho_max]);

%% ============================================================
% Timing summary
%% ============================================================
fprintf('\nDone.\n');
fprintf('King z-view = %.6f m\n', z_view_K);
fprintf('DIM  z-view = %.6f m\n', z_view_D);

fprintf('\n==== Timing summary ====\n');
fprintf('King total time            : %.6f s\n', tKing);
fprintf('King grid points           : %d\n', NKing);
fprintf('King average time / point  : %.6e s/point\n', tKing_per_point);
fprintf('King average time / point  : %.6f ms/point\n', tKing_per_point * 1e3);
fprintf('King average time / z-slice: %.6f s/slice\n', tKing_per_z);

fprintf('DIM total time             : %.6f s\n', tDIM);
fprintf('DIM observation points     : %d\n', N_dim_obs);
fprintf('DIM average time / point   : %.6e s/point\n', tDIM_per_point);
fprintf('DIM average time / point   : %.6f ms/point\n', tDIM_per_point * 1e3);
fprintf('DIM average time / z-slice : %.6f s/slice\n', tDIM_per_z);