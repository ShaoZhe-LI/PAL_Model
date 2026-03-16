%% ============================================================
% MAIN: f1 ultrasound only (2D only)
% - King: compute on its own FHT grid
% - DIM : use the same FHT grid, but rho is downsampled
%         OR use uniform sampling
% - z range: 0 ~ 1 m
%
% Outputs:
%   1) xOz amplitude / phase (King vs DIM)
%   2) radial amplitude / phase at z = 1 m (King vs DIM)
%   3) xOz real / imag (King vs DIM)
%   4) radial real / imag at z = 1 m (King vs DIM)
%   5) timing summary for both King and DIM
%   6) save all figures to a result folder
%   7) save all calculation parameters to a txt file
%
% DIM sampling modes:
%   - 'fht_downsample' : use FHT-aligned rho/z, with rho downsampled
%   - 'uniform'        : use uniform x/z sampling
%
% 文件夹结构：
%   ./Audio_test/
%       └─ 12.2kHz_32768_FHTDownsample_20260316_1009/
%             ├─ Fig1_xOz_amp_phase.png/.fig
%             ├─ Fig2_radial_amp_phase.png/.fig
%             ├─ Fig3_xOz_real_imag.png/.fig
%             ├─ Fig4_radial_real_imag.png/.fig
%             └─ params_and_summary.txt
%  本代码用于调试不同采样下dim与fht方法的差别
clear; clc; close all;

%% -------------------- save control --------------------
global SAVE_PNG SAVE_DIR
SAVE_PNG = false;  % 'true' or 'false'
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
source.m       = 0;
source.F       = 0.2;

source.f1 = 12e3;
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
calc.fht.delta      = medium.c0 / source.f1 / 2;

% --- DIM source discretization ---
calc.dim.method             = 'rayleigh';   % 'rayleigh' or 'asm'
calc.dim.use_freq           = 'f1';
calc.dim.dis_coe            = 16;
calc.dim.margin             = 1;
calc.dim.src_discretization = 'polar';
calc.dim.use_parallel       = true;
calc.dim.num_workers        = 4;

% --- DIM observation-grid sampling mode ---
calc.dim.grid_mode = 'fht_downsample';   % 'fht_downsample' or 'uniform'

% used when calc.dim.grid_mode = 'fht_downsample'
calc.dim.ds_rho = 32;   % rho downsample factor for DIM, 32For20kHz/16ForMore

% used when calc.dim.grid_mode = 'uniform'
calc.dim.uniform_dx = calc.fht.delta / 2;   % /2
calc.dim.uniform_dz = calc.fht.delta / 1;   % /1

% observation block size
calc.dim.block_size = 2e5;

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
z_full   = z_full(z_full <= calc.fht.zu_max);

% nearest z to target
[~, iz_view_full] = min(abs(z_full - z_target));
z_view = z_full(iz_view_full);

%% ============================================================
% BUILD RESULT FOLDER
%% ============================================================
root_save_dir = fullfile(pwd, 'Audio_test');
if ~exist(root_save_dir, 'dir')
    mkdir(root_save_dir);
end

f1_kHz = source.f1 / 1e3;

switch lower(calc.dim.grid_mode)
    case 'fht_downsample'
        grid_mode_tag = 'FHTDownsample';
    case 'uniform'
        grid_mode_tag = 'Uniform';
    otherwise
        grid_mode_tag = calc.dim.grid_mode;
end

time_tag = datestr(now, 'yyyymmdd_HHMM');
result_folder_name = sprintf('%.1fkHz_%d_%s_%s', ...
    f1_kHz, calc.fht.N_FHT, grid_mode_tag, time_tag);

SAVE_DIR = fullfile(root_save_dir, result_folder_name);
if ~exist(SAVE_DIR, 'dir')
    mkdir(SAVE_DIR);
end

fprintf('Root folder   : %s\n', root_save_dir);
fprintf('Result folder : %s\n', SAVE_DIR);

%% ============================================================
% BUILD DIM OBSERVATION GRID
%% ============================================================
switch lower(calc.dim.grid_mode)
    case 'fht_downsample'
        ds_rho = calc.dim.ds_rho;

        rho_ds = rho_full(1:ds_rho:end);
        if rho_ds(end) ~= rho_full(end)
            rho_ds = [rho_ds; rho_full(end)];
        end

        obs_grid = struct();
        obs_grid.dim.x = rho_ds(:).';     % FHT rho grid with downsampling
        obs_grid.dim.y = 0;
        obs_grid.dim.z = z_full(:).';     % FHT z grid, cropped to 0~zu_max
        obs_grid.dim.block_size = calc.dim.block_size;

        fprintf('DIM grid mode         : %s\n', calc.dim.grid_mode);
        fprintf('FHT full rho points   : %d\n', numel(rho_full));
        fprintf('DIM rho points        : %d (downsample factor = %d)\n', numel(rho_ds), ds_rho);
        fprintf('z points              : %d, z in [0, %.6f] m\n', numel(z_full), z_full(end));
        fprintf('z_view                : %.6f m\n', z_view);

    case 'uniform'
        dxyz_x = calc.dim.uniform_dx;
        dxyz_z = calc.dim.uniform_dz;

        obs_grid = struct();
        obs_grid.dim.x = 0:dxyz_x:calc.fht.rho_max;
        obs_grid.dim.y = 0;
        obs_grid.dim.z = 0:dxyz_z:calc.fht.zu_max;
        obs_grid.dim.block_size = calc.dim.block_size;

        fprintf('DIM grid mode         : %s\n', calc.dim.grid_mode);
        fprintf('Uniform dx            : %.6e m\n', dxyz_x);
        fprintf('Uniform dz            : %.6e m\n', dxyz_z);
        fprintf('DIM x points          : %d\n', numel(obs_grid.dim.x));
        fprintf('DIM z points          : %d\n', numel(obs_grid.dim.z));
        fprintf('z_view                : %.6f m\n', z_view);

    otherwise
        error('Unknown calc.dim.grid_mode: %s. Use ''fht_downsample'' or ''uniform''.', calc.dim.grid_mode);
end

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
id_zK = find(zK_all   <= calc.fht.zu_max,  1, 'last');
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
% DIM: use selected observation grid
%% ============================================================
fprintf('\n==== DIM f1 field ====\n');

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
fig1 = figure('Name','xOz |p| and phase: King vs DIM', 'Position',[80 60 1500 900]);

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

save_current_figure(fig1, SAVE_DIR, 'Fig1_xOz_amp_phase');

%% ============================================================
% Plot 2: radial amplitude / phase at z = z_view
%% ============================================================
fig2 = figure('Name','Radial profile at z=1m', 'Position',[150 120 1200 480]);

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

save_current_figure(fig2, SAVE_DIR, 'Fig2_radial_amp_phase');

%% ============================================================
% Plot 3: xOz real / imag
%% ============================================================
fig3 = figure('Name','xOz real and imag: King vs DIM', 'Position',[100 70 1500 900]);

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

save_current_figure(fig3, SAVE_DIR, 'Fig3_xOz_real_imag');

%% ============================================================
% Plot 4: radial real / imag at z = z_view
%% ============================================================
fig4 = figure('Name','Radial real and imag at z=1m', 'Position',[170 130 1200 480]);

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

save_current_figure(fig4, SAVE_DIR, 'Fig4_radial_real_imag');

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

%% ============================================================
% Save parameters and summary to txt
%% ============================================================
param_txt_file = fullfile(SAVE_DIR, 'params_and_summary.txt');
write_params_and_summary_txt(param_txt_file, ...
    medium, source, calc, fig, z_target, ...
    z_view, z_view_K, z_view_D, ...
    rho_full, z_full, obs_grid, ...
    tKing, NKing, tKing_per_point, tKing_per_z, ...
    tDIM, N_dim_obs, tDIM_per_point, tDIM_per_z);

fprintf('\nSaved figures and parameter txt to:\n%s\n', SAVE_DIR);

%% ============================================================
% Local functions
%% ============================================================
function save_current_figure(fig_handle, save_dir, file_stem)
    if ~ishandle(fig_handle), return; end
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    png_file = fullfile(save_dir, [file_stem, '.png']);
    fig_file = fullfile(save_dir, [file_stem, '.fig']);

    try
        exportgraphics(fig_handle, png_file, 'Resolution', 300);
    catch
        saveas(fig_handle, png_file);
    end

    try
        savefig(fig_handle, fig_file);
    catch
        warning('Failed to save .fig file: %s', fig_file);
    end
end

function write_params_and_summary_txt(txt_file, ...
    medium, source, calc, fig, z_target, ...
    z_view, z_view_K, z_view_D, ...
    rho_full, z_full, obs_grid, ...
    tKing, NKing, tKing_per_point, tKing_per_z, ...
    tDIM, N_dim_obs, tDIM_per_point, tDIM_per_z)

    fid = fopen(txt_file, 'w');
    if fid < 0
        error('Cannot open txt file for writing: %s', txt_file);
    end

    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    % ===== wavelength / delta relations =====
    lambda_u  = medium.c0 / source.f1;
    lambda_a  = medium.c0 / source.fa;
    delta_val = calc.fht.delta;

    delta_over_lambda_u = delta_val / lambda_u;
    delta_over_lambda_a = delta_val / lambda_a;

    % ===== uniform spacing / delta relations =====
    uniform_dx_over_delta = calc.dim.uniform_dx / delta_val;
    uniform_dz_over_delta = calc.dim.uniform_dz / delta_val;

    fprintf(fid, '============================================================\n');
    fprintf(fid, 'Calculation parameters and summary\n');
    fprintf(fid, 'Generated time: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '============================================================\n\n');

    %% -------------------- medium --------------------
    fprintf(fid, '[medium]\n');
    fprintf(fid, 'c0               = %.16g\n', medium.c0);
    fprintf(fid, 'rho0             = %.16g\n', medium.rho0);
    fprintf(fid, 'beta             = %.16g\n', medium.beta);
    fprintf(fid, 'pref             = %.16g\n', medium.pref);
    fprintf(fid, 'use_absorp       = %d\n', medium.use_absorp);
    fprintf(fid, 'atten_handle     = %s\n\n', func2str(medium.atten_handle));

    %% -------------------- source --------------------
    fprintf(fid, '[source]\n');
    fprintf(fid, 'profile          = %s\n', source.profile);
    fprintf(fid, 'a                = %.16g\n', source.a);
    fprintf(fid, 'v0               = %.16g\n', source.v0);
    fprintf(fid, 'v_ratio          = %.16g\n', source.v_ratio);
    fprintf(fid, 'm                = %.16g\n', source.m);
    fprintf(fid, 'F                = %.16g\n', source.F);
    fprintf(fid, 'f1               = %.16g\n', source.f1);
    fprintf(fid, 'fa               = %.16g\n', source.fa);
    fprintf(fid, 'f2               = %.16g\n', source.f2);
    fprintf(fid, 'lambda_u (=c0/f1)= %.16g\n', lambda_u);
    fprintf(fid, 'lambda_a (=c0/fa)= %.16g\n\n', lambda_a);

    %% -------------------- calc.fht --------------------
    fprintf(fid, '[calc.fht]\n');
    fprintf(fid, 'N_FHT            = %d\n', calc.fht.N_FHT);
    fprintf(fid, 'rho_max          = %.16g\n', calc.fht.rho_max);
    fprintf(fid, 'Nh_scale         = %.16g\n', calc.fht.Nh_scale);
    fprintf(fid, 'NH_scale         = %.16g\n', calc.fht.NH_scale);
    fprintf(fid, 'Nh_v_scale       = %.16g\n', calc.fht.Nh_v_scale);
    fprintf(fid, 'zu_max           = %.16g\n', calc.fht.zu_max);
    fprintf(fid, 'za_max           = %.16g\n', calc.fht.za_max);
    fprintf(fid, 'delta            = %.16g\n', delta_val);
    fprintf(fid, 'delta/lambda_u   = %.16g\n', delta_over_lambda_u);
    fprintf(fid, 'delta/lambda_a   = %.16g\n', delta_over_lambda_a);
    fprintf(fid, 'lambda_u/delta   = %.16g\n', lambda_u / delta_val);
    fprintf(fid, 'lambda_a/delta   = %.16g\n\n', lambda_a / delta_val);

    %% -------------------- calc.dim --------------------
    fprintf(fid, '[calc.dim]\n');
    fprintf(fid, 'method               = %s\n', calc.dim.method);
    fprintf(fid, 'use_freq             = %s\n', calc.dim.use_freq);
    fprintf(fid, 'dis_coe              = %.16g\n', calc.dim.dis_coe);
    fprintf(fid, 'margin               = %.16g\n', calc.dim.margin);
    fprintf(fid, 'src_discretization   = %s\n', calc.dim.src_discretization);
    fprintf(fid, 'use_parallel         = %d\n', calc.dim.use_parallel);
    fprintf(fid, 'num_workers          = %d\n', calc.dim.num_workers);
    fprintf(fid, 'grid_mode            = %s\n', calc.dim.grid_mode);
    fprintf(fid, 'ds_rho               = %.16g\n', calc.dim.ds_rho);
    fprintf(fid, 'uniform_dx           = %.16g\n', calc.dim.uniform_dx);
    fprintf(fid, 'uniform_dx/delta     = %.16g\n', uniform_dx_over_delta);
    fprintf(fid, 'delta/uniform_dx     = %.16g\n', delta_val / calc.dim.uniform_dx);
    fprintf(fid, 'uniform_dz           = %.16g\n', calc.dim.uniform_dz);
    fprintf(fid, 'uniform_dz/delta     = %.16g\n', uniform_dz_over_delta);
    fprintf(fid, 'delta/uniform_dz     = %.16g\n', delta_val / calc.dim.uniform_dz);
    fprintf(fid, 'block_size           = %.16g\n\n', calc.dim.block_size);

    %% -------------------- calc.king --------------------
    fprintf(fid, '[calc.king]\n');
    fprintf(fid, 'gspec_method       = %s\n', calc.king.gspec_method);
    fprintf(fid, 'eps_phase          = %.16g\n', calc.king.eps_phase);
    fprintf(fid, 'kz_min             = %.16g\n', calc.king.kz_min);
    fprintf(fid, 'band_refine.enable = %d\n\n', calc.king.band_refine.enable);

    %% -------------------- calc.asm --------------------
    fprintf(fid, '[calc.asm]\n');
    fprintf(fid, 'pad_factor         = %.16g\n', calc.asm.pad_factor);
    fprintf(fid, 'kzz_eps            = %.16g\n\n', calc.asm.kzz_eps);

    %% -------------------- figure settings --------------------
    fprintf(fid, '[fig]\n');
    fprintf(fid, 'unwrap_phase_1d    = %d\n', fig.unwrap_phase_1d);
    fprintf(fid, 'mask_low_amp_phase = %d\n', fig.mask_low_amp_phase);
    fprintf(fid, 'phase_amp_ratio    = %.16g\n', fig.phase_amp_ratio);
    fprintf(fid, 'z_target           = %.16g\n', z_target);
    fprintf(fid, 'z_view_used        = %.16g\n', z_view);
    fprintf(fid, 'King z_view        = %.16g\n', z_view_K);
    fprintf(fid, 'DIM  z_view        = %.16g\n\n', z_view_D);

    %% -------------------- prepared FHT grid --------------------
    fprintf(fid, '[prepared FHT grid]\n');
    fprintf(fid, 'rho_full points    = %d\n', numel(rho_full));
    fprintf(fid, 'z_full points      = %d\n', numel(z_full));
    if ~isempty(rho_full)
        fprintf(fid, 'rho_full min/max   = %.16g, %.16g\n', rho_full(1), rho_full(end));
    end
    if ~isempty(z_full)
        fprintf(fid, 'z_full min/max     = %.16g, %.16g\n', z_full(1), z_full(end));
    end
    fprintf(fid, '\n');

    %% -------------------- DIM observation grid --------------------
    fprintf(fid, '[DIM observation grid]\n');
    fprintf(fid, 'Nx                = %d\n', numel(obs_grid.dim.x));
    fprintf(fid, 'Ny                = %d\n', numel(obs_grid.dim.y));
    fprintf(fid, 'Nz                = %d\n', numel(obs_grid.dim.z));
    if ~isempty(obs_grid.dim.x)
        fprintf(fid, 'x min/max         = %.16g, %.16g\n', obs_grid.dim.x(1), obs_grid.dim.x(end));
    end
    if ~isempty(obs_grid.dim.z)
        fprintf(fid, 'z min/max         = %.16g, %.16g\n', obs_grid.dim.z(1), obs_grid.dim.z(end));
    end
    fprintf(fid, 'block_size        = %.16g\n\n', obs_grid.dim.block_size);

    %% -------------------- timing summary --------------------
    fprintf(fid, '[timing summary]\n');
    fprintf(fid, 'King total time            = %.16g s\n', tKing);
    fprintf(fid, 'King grid points           = %d\n', NKing);
    fprintf(fid, 'King average time / point  = %.16g s/point\n', tKing_per_point);
    fprintf(fid, 'King average time / point  = %.16g ms/point\n', tKing_per_point * 1e3);
    fprintf(fid, 'King average time / z      = %.16g s/slice\n', tKing_per_z);
    fprintf(fid, '\n');
    fprintf(fid, 'DIM total time             = %.16g s\n', tDIM);
    fprintf(fid, 'DIM observation points     = %d\n', N_dim_obs);
    fprintf(fid, 'DIM average time / point   = %.16g s/point\n', tDIM_per_point);
    fprintf(fid, 'DIM average time / point   = %.16g ms/point\n', tDIM_per_point * 1e3);
    fprintf(fid, 'DIM average time / z       = %.16g s/slice\n', tDIM_per_z);
    fprintf(fid, '\n');

    fprintf(fid, '============================================================\n');
end