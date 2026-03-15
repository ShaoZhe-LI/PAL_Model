%% ============================================================
% MAIN: ultrasound + audio
%
% Ultrasound:
%   - King on its own FHT grid
%   - DIM on the same z grid and downsampled rho grid
%
% Audio:
%   - King: pa_W using your m_FHT-based workflow
%   - DIM : only one radial line at z = z_goal by direct triple integral
%
% Notes:
%   1) pa_W is used as the target audio field
%   2) DIM audio is only evaluated on one line z = z_goal
%   3) King audio uses analytic Green-function spectrum
%   4) FHT and m_FHT are treated as identical, so only m_FHT is used here
%% ============================================================
clear; clc; close all;

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

source.f1 = 42e3;
source.fa = 2e3;
source.f2 = source.f1 + source.fa;

% separate orders for f1 and f2
if ~isfield(source,'m1'), source.m1 = source.m; end
if ~isfield(source,'m2'), source.m2 = source.m; end

%% -------------------- calc --------------------
calc = struct();

% --- FHT / King ---
calc.fht.N_FHT      = 32768;
calc.fht.rho_max    = 0.8;
calc.fht.Nh_scale   = 1.2;
calc.fht.NH_scale   = 4;
calc.fht.Nh_v_scale = 1.1;
calc.fht.zu_max     = 1.0;   % ultrasound z-range
calc.fht.za_max     = 1.0;   % audio z-range
% optional external delta, if your make_source_velocity supports it:
% calc.fht.delta    = medium.c0 / source.f1 / 2;

% --- DIM source discretization ---
calc.dim.method             = 'rayleigh';   % 'rayleigh' or 'asm'
calc.dim.use_freq           = 'f2';
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

% --- ASM settings ---
calc.asm.pad_factor = 16;
calc.asm.kzz_eps    = 1e-12;

%% -------------------- display / comparison setup --------------------
z_target = 1.0;    % ultrasound 1D profile
z_goal   = 1.0;    % audio line
ds_rho   = 32;     % rho downsample factor for DIM

fig.unwrap_phase_1d    = true;
fig.mask_low_amp_phase = false;
fig.phase_amp_ratio    = 1e-3;

%% -------------------- audio setup --------------------
calc.audio.enable           = true;
calc.audio.z_goal           = z_goal;
calc.audio.Nphi_dim         = 181;      % azimuth samples for DIM triple integral
calc.audio.use_parallel_dim = false;    % optional parfor over obs points

%% ============================================================
% PREPARE FHT GRID FIRST
%% ============================================================
[~, fht_tmp] = make_source_velocity(source, medium, calc);

rho_full = fht_tmp.xh(:);
z_full   = fht_tmp.z_ultra(:);

rho_full = rho_full(rho_full <= calc.fht.rho_max);
z_full   = z_full(z_full <= calc.fht.zu_max);

[~, iz_view_full] = min(abs(z_full - z_target));
z_view = z_full(iz_view_full);

rho_ds = rho_full(1:ds_rho:end);
if rho_ds(end) ~= rho_full(end)
    rho_ds = [rho_ds; rho_full(end)];
end

fprintf('FHT full rho points : %d\n', numel(rho_full));
fprintf('DIM rho points      : %d (downsample factor ~ %d)\n', numel(rho_ds), ds_rho);
fprintf('z points            : %d, z in [0, %.6f] m\n', numel(z_full), z_full(end));
fprintf('z_view              : %.6f m\n', z_view);

%% ============================================================
% KING: ultrasound f1/f2
%% ============================================================
fprintf('\n==== KING ultrasound field ====\n');
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

pK1_rz = resK.king.p_f1(1:id_rK, 1:id_zK);
pK2_rz = resK.king.p_f2(1:id_rK, 1:id_zK);

% for ultrasound plots
pK_rz = pK1_rz;

[~, izK] = min(abs(zK - z_view));
z_view_K = zK(izK);
pK_r = pK1_rz(:, izK);

NKing = numel(rhoK) * numel(zK);
tKing_per_point = tKing / NKing;
tKing_per_z     = tKing / numel(zK);

%% ============================================================
% DIM: ultrasound f1/f2
%% ============================================================
fprintf('\n==== DIM ultrasound field ====\n');

obs_grid = struct();
obs_grid.dim.x = rho_ds(:).';
obs_grid.dim.y = 0;
obs_grid.dim.z = z_full(:).';
obs_grid.dim.block_size = 2e5;

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
    pD1_rz = squeeze(resD.dim.p_f1(1,:,:));   % Nx x Nz
    pD2_rz = squeeze(resD.dim.p_f2(1,:,:));   % Nx x Nz
    rhoD   = resD.dim.x(:);
    zD     = resD.dim.z(:).';
else
    xA = resD.dim.x(:);
    yA = resD.dim.y(:);
    zA = resD.dim.z(:).';
    [~, iy0] = min(abs(yA - 0));
    pD1_rz = squeeze(resD.dim.p_f1(iy0,:,:));
    pD2_rz = squeeze(resD.dim.p_f2(iy0,:,:));
    rhoD   = xA;
    zD     = zA;
end

% for ultrasound plots
pD_rz = pD1_rz;

[~, izD] = min(abs(zD - z_view));
z_view_D = zD(izD);
pD_r = pD1_rz(:, izD);

%% ============================================================
% Ultrasound xOz maps
%% ============================================================
xK_xoz = [-flipud(rhoK); rhoK];
pK_xoz = [flipud(pK_rz); pK_rz];

xD_xoz = [-flipud(rhoD); rhoD];
pD_xoz = [flipud(pD_rz); pD_rz];

AK_xoz  = abs(pK_xoz);
AD_xoz  = abs(pD_xoz);
PhK_xoz = angle(pK_xoz);
PhD_xoz = angle(pD_xoz);

if fig.mask_low_amp_phase
    thK_xoz = max(AK_xoz(:)) * fig.phase_amp_ratio;
    thD_xoz = max(AD_xoz(:)) * fig.phase_amp_ratio;
    PhK_xoz(AK_xoz < thK_xoz) = NaN;
    PhD_xoz(AD_xoz < thD_xoz) = NaN;
end

%% ============================================================
% Ultrasound radial phase for 1D plot
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

%% ============================================================
% Ultrasound plots
%% ============================================================
amp_lim_xoz = [0, max([AK_xoz(:); AD_xoz(:)])];
ph_lim      = [-pi, pi];

figure('Name','Ultrasound xOz |p| and phase: King vs DIM', 'Position',[80 60 1500 900]);

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

figure('Name','Ultrasound radial profile at z=1m', 'Position',[150 120 1200 480]);

subplot(1,2,1);
plot(rhoK, abs(pK_r), 'LineWidth', 1.8); hold on;
plot(rhoD, abs(pD_r), '--', 'LineWidth', 1.8);
grid on;
xlabel('\rho (m)');
ylabel('|p_{f1}|');
title(sprintf('Ultrasound amplitude at z = %.6f m', z_view));
legend('King','DIM','Location','best');
xlim([0 calc.fht.rho_max]);

subplot(1,2,2);
plot(rhoK, phK_r, 'LineWidth', 1.8); hold on;
plot(rhoD, phD_r, '--', 'LineWidth', 1.8);
grid on;
xlabel('\rho (m)');
ylabel('Phase (rad)');
title(sprintf('Ultrasound phase at z = %.6f m', z_view));
legend('King','DIM','Location','best');
xlim([0 calc.fht.rho_max]);

%% ============================================================
% Audio calculation
%% ============================================================
if calc.audio.enable
    fprintf('\n==== AUDIO calculation ====\n');
    fprintf('King audio Green spectrum: analytic form Gr00 = exp(1j*kza*abs(z))/kza\n');

    % --------------------------------------------------------
    % King audio pa_W on full (rho,z) grid using your m_FHT workflow
    % --------------------------------------------------------
    tic;
    [paK_W_rz, qK_rz, rhoAK, zAK] = local_audio_king_paW_using_user_mFHT( ...
        rhoK, zK, pK1_rz, pK2_rz, source, medium, calc);
    tAudioK = toc;

    % --------------------------------------------------------
    % DIM audio pa_W only on z = z_goal line by triple integral
    % --------------------------------------------------------
    tic;
    [paD_W_line, qD_rz] = local_audio_dim_paW_line_triple( ...
        rhoD, zD, pD1_rz, pD2_rz, source, medium, calc.audio);
    tAudioD = toc;

    % King line for comparison, interpolated to DIM rho points
    [~, izAK] = min(abs(zAK - calc.audio.z_goal));
    paK_W_line_onD = interp1(rhoAK, paK_W_rz(:,izAK), rhoD, 'linear', 0);

    % --------------------------------------------------------
    % King audio xOz plot
    % --------------------------------------------------------
    xAK_xoz   = [-flipud(rhoAK); rhoAK];
    paK_W_xoz = [flipud(paK_W_rz); paK_W_rz];

    SPL_aK_xoz = 20*log10(abs(paK_W_xoz) / medium.pref + eps);
    ph_aK_xoz  = angle(paK_W_xoz);

    figure('Name','Audio (King pa_W) on xOz', 'Position',[100 100 1200 500]);

    subplot(1,2,1);
    pcolor(zAK, xAK_xoz, SPL_aK_xoz); shading interp;
    colormap(gca, parula(256)); colorbar;
    xlabel('z (m)'); ylabel('x (m)');
    title('King audio pa_W (SPL) on xOz');
    xlim([0 zAK(end)]); ylim([-max(rhoAK) max(rhoAK)]);
    pbaspect([zAK(end), 2*max(rhoAK), 1]);

    subplot(1,2,2);
    pcolor(zAK, xAK_xoz, ph_aK_xoz); shading interp;
    colormap(gca, hsv(256)); clim([-pi pi]); colorbar;
    xlabel('z (m)'); ylabel('x (m)');
    title('King audio pa_W phase on xOz');
    xlim([0 zAK(end)]); ylim([-max(rhoAK) max(rhoAK)]);
    pbaspect([zAK(end), 2*max(rhoAK), 1]);

    % --------------------------------------------------------
    % Audio line comparison at z = z_goal
    % --------------------------------------------------------
    figure('Name','Audio line comparison at z=1m', 'Position',[160 120 1200 480]);

    subplot(1,2,1);
    plot(rhoD, abs(paK_W_line_onD), 'LineWidth', 1.8); hold on;
    plot(rhoD, abs(paD_W_line), '--', 'LineWidth', 1.8);
    grid on;
    xlabel('\rho (m)');
    ylabel('|p_a|');
    title(sprintf('Audio magnitude at z = %.6f m', calc.audio.z_goal));
    legend('King pa_W','DIM triple-integral','Location','best');
    xlim([0 max(rhoD)]);

    subplot(1,2,2);
    plot(rhoD, unwrap(angle(paK_W_line_onD)), 'LineWidth', 1.8); hold on;
    plot(rhoD, unwrap(angle(paD_W_line)), '--', 'LineWidth', 1.8);
    grid on;
    xlabel('\rho (m)');
    ylabel('Phase (rad)');
    title(sprintf('Audio phase at z = %.6f m', calc.audio.z_goal));
    legend('King pa_W','DIM triple-integral','Location','best');
    xlim([0 max(rhoD)]);

    % audio time stats
    N_audio_king = numel(paK_W_rz);
    N_audio_dim  = numel(paD_W_line);

    tAudioK_per_point = tAudioK / N_audio_king;
    tAudioD_per_point = tAudioD / N_audio_dim;

    fprintf('Audio King pa_W total time       : %.6f s\n', tAudioK);
    fprintf('Audio King audio points          : %d\n', N_audio_king);
    fprintf('Audio King avg time / point      : %.6e s/point\n', tAudioK_per_point);
    fprintf('Audio King avg time / point      : %.6f ms/point\n', tAudioK_per_point * 1e3);

    fprintf('Audio DIM line total time        : %.6f s\n', tAudioD);
    fprintf('Audio DIM audio points           : %d\n', N_audio_dim);
    fprintf('Audio DIM avg time / point       : %.6e s/point\n', tAudioD_per_point);
    fprintf('Audio DIM avg time / point       : %.6f ms/point\n', tAudioD_per_point * 1e3);
end

%% ============================================================
% Timing summary
%% ============================================================
fprintf('\nDone.\n');
fprintf('King z-view = %.6f m\n', z_view_K);
fprintf('DIM  z-view = %.6f m\n', z_view_D);

fprintf('\n==== Ultrasound timing summary ====\n');
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
% Local functions
%% ============================================================

function [pa_W, q, xh, z_audio] = local_audio_king_paW_using_user_mFHT(rho_in, z_in, p1_audio, p2_audio, source, medium, calc)
% Use your original solve_kappa0 + m_FHT workflow to compute pa_W
% Output pa_W is on (xh, z_audio)

c    = medium.c0;
rho0 = medium.rho0;
beta = medium.beta;

f1 = source.f1;
f2 = source.f2;
fa = source.fa;

wa = 2*pi*fa;
ka = wa/c + 1j*medium.atten_handle(fa);

m1 = source.m1;
m2 = source.m2;
ma = m2 - m1;

N_FHT   = calc.fht.N_FHT;
rho_max = calc.fht.rho_max;
za_max  = calc.fht.za_max;

% delta from current FHT setup
[~, fht_tmp] = make_source_velocity(source, medium, calc);
delta = fht_tmp.delta;

% ------------------------------------------------------------
% Your FHT-grid construction
% ------------------------------------------------------------
Nh = 1.2 * rho_max;
NH = 1.2 * (2*pi*f2) / c;

n_FHT = 0:N_FHT-1;
[a_solve, k0, x1, x0] = solve_kappa0(N_FHT, n_FHT);

xh = (x1 * Nh).';
xH = (x1 * NH).';

% crop rho to needed range
idx_r = find(xh <= rho_max, 1, 'last');
if isempty(idx_r), idx_r = numel(xh); end

z_audio = 0:delta:za_max;
Nza = numel(z_audio);

% ------------------------------------------------------------
% Interpolate provided ultrasound fields onto xh x z_audio
% ------------------------------------------------------------
[Zin, Rin]   = meshgrid(z_in(:).', rho_in(:));
[Zout, Rout] = meshgrid(z_audio(:).', xh(:));

p1a_full = interp2(Zin, Rin, p1_audio, Zout, Rout, 'linear', 0);
p2a_full = interp2(Zin, Rin, p2_audio, Zout, Rout, 'linear', 0);

% virtual source density, same as your reference code
q_full = conj(p1a_full) .* p2a_full * beta * wa / (1i * rho0^2 * c^4);

% ------------------------------------------------------------
% pa_W using your workflow, all transforms by m_FHT
% ------------------------------------------------------------
Qr00 = m_FHT(q_full, N_FHT, Nza, Nh, NH, a_solve, x0, x1, k0, ma);
Qr0  = [fliplr(Qr00(:,2:end)) Qr00];

Nz1    = size(Qr0,2);
N_conv = Nz1 + Nza - 1;

[M0, N0] = size(Qr0);
Qr = [Qr0, zeros(M0, N_conv - N0)];
Q  = (fft(Qr.')).';

% analytic Green spectrum
kza  = sqrt(ka^2 - xH.^2);
Gr00 = exp(1j * kza .* abs(z_audio)) ./ kza;
Gr0  = [fliplr(Gr00(:,2:end)) Gr00];
Gr   = [Gr0, zeros(size(Gr0,1), N_conv - size(Gr0,2))];
G    = (fft(Gr.')).';

Pa   = Q .* G;
par0 = (ifft(Pa.')).';
par  = par0(:, N_conv - Nza + 1 : N_conv);

phia_W  = -m_FHT(par, N_FHT, Nza, NH, Nh, a_solve, x0, x1, k0, ma) * delta * 1j / 2;
pa_W_full = 1j * rho0 * wa * phia_W;

% crop to rho <= rho_max
pa_W = pa_W_full(1:idx_r, :);
q    = q_full(1:idx_r, :);
xh   = xh(1:idx_r);
end

function [paW_line, q_rz] = local_audio_dim_paW_line_triple(rho, z, p1_rz, p2_rz, source, medium, audio_cfg)
% DIM audio only on the line (rho_obs, phi=0, z=z_goal)
% direct triple integral in cylindrical coordinates

rho = rho(:);
z   = z(:).';

Nr = numel(rho);
Nz = numel(z);

c0   = medium.c0;
rho0 = medium.rho0;
beta = medium.beta;

wa = 2*pi*source.fa;
ka = wa/c0 + 1i*medium.atten_handle(source.fa);

m1 = source.m1;
m2 = source.m2;
ma = m2 - m1;

z_goal = audio_cfg.z_goal;
Nphi   = audio_cfg.Nphi_dim;

q_rz = conj(p1_rz) .* p2_rz * beta * wa / (1i * rho0^2 * c0^4);

rho_obs = rho(:);

wrho = local_trapz_weights(rho);
wz   = local_trapz_weights(z);

phi_s = linspace(0, 2*pi, Nphi+1);
phi_s(end) = [];
dphi = 2*pi / Nphi;

eimphi = exp(1i * ma * phi_s);

phi_line = complex(zeros(size(rho_obs)));

use_par = isfield(audio_cfg,'use_parallel_dim') && audio_cfg.use_parallel_dim;

if use_par
    parfor io = 1:numel(rho_obs)
        xo = rho_obs(io);
        acc = 0;
        for iz = 1:Nz
            dz_w = wz(iz);
            zz = z(iz);
            q_col = q_rz(:,iz);
            dzeta = z_goal - zz;
            for ir = 1:Nr
                rs = rho(ir);
                base = q_col(ir) * rs * wrho(ir) * dz_w;
                R = sqrt(xo^2 + rs^2 - 2*xo*rs*cos(phi_s) + dzeta^2);
                G = exp(1i*ka*R) ./ (4*pi*R);
                acc = acc + base * sum(eimphi .* G) * dphi;
            end
        end
        phi_line(io) = acc;
    end
else
    for io = 1:numel(rho_obs)
        xo = rho_obs(io);
        acc = 0;
        for iz = 1:Nz
            dz_w = wz(iz);
            zz = z(iz);
            q_col = q_rz(:,iz);
            dzeta = z_goal - zz;
            for ir = 1:Nr
                rs = rho(ir);
                base = q_col(ir) * rs * wrho(ir) * dz_w;
                R = sqrt(xo^2 + rs^2 - 2*xo*rs*cos(phi_s) + dzeta^2);
                G = exp(1i*ka*R) ./ (4*pi*R);
                acc = acc + base * sum(eimphi .* G) * dphi;
            end
        end
        phi_line(io) = acc;
    end
end

paW_line = 1i * rho0 * wa * phi_line;
end

function w = local_trapz_weights(x)
x = x(:);
N = numel(x);
w = zeros(N,1);

if N == 1
    w(1) = 1;
    return;
end

w(1) = (x(2) - x(1)) / 2;
w(N) = (x(N) - x(N-1)) / 2;
for n = 2:N-1
    w(n) = (x(n+1) - x(n-1)) / 2;
end
end