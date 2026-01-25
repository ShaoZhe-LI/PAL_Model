%% ============================================================
% MAIN: King vs DIM (memory-lean) comparison + profiler (time & memory)
% + EXTRA FIGS (FHT only):
% (1) xOz (from King/FHT): |p| and SPL on (x,z)
% (2) xOy @ z≈1 m (from King/FHT): SPL and phase on (x,y)
% via planar phase extension: p(x,y)=p(rho)*exp(j*m*phi)
% ============================================================

clear; clc; close all;

%% -------------------- medium --------------------
medium.c0 = 343;
medium.rho0 = 1.21;
medium.beta = 1.2;
medium.pref = 2e-5;
medium.use_absorp = true;
medium.atten_handle = @(f) AbsorpAttenCoef(f);

%% -------------------- source --------------------
source.profile = 'Vortex-m'; % 'Uniform' | 'Focus' | 'Vortex-m' | 'Poly' | 'Custom'
source.a = 0.1;
source.v0 = 0.172;
source.v_ratio = 1;
source.m = 3;
source.F = 0.2;

source.f1 = 42e3;
source.fa = 2e3;

%% -------------------- calc parameters --------------------
calc = struct();

% --- FHT / King ---
calc.fht.N_FHT = 32768;
calc.fht.rho_max = 1.2;
calc.fht.Nh_scale = 1.2;
calc.fht.NH_scale = 4;
calc.fht.Nh_v_scale = 1.1;
calc.fht.zu_max = 4;
calc.fht.za_max = 2;

% --- DIM source discretization ---
calc.dim.use_freq = 'f2';
calc.dim.dis_coe = 16;
calc.dim.margin = 1;
calc.dim.method = 'rayleigh'; % 'rayleigh' or 'asm'

% --- King analytic spectrum stability ---
calc.king.method = 'analytic'; % 'analytic' or 'transform'
calc.king.eps_kzz = 1e-3; % keep your naming
% map to function-expected fields:
calc.king.eps_phase = calc.king.eps_kzz;
calc.king.kz_min = 1e-12;

% --- ASM settings (used only if calc.dim.method='asm') ---
calc.asm.pad_factor = 16;
calc.asm.kzz_eps = 1e-12;

%% -------------------- target plane --------------------
z_target = 1.0;
ds = 16;

%% ==================== PROFILER: helpers ====================
get_mem_mb = @() local_get_mem_mb(); % returns NaN if not available
fmt_mem = @(x) local_fmt_mem(x);

%% -------------------- King only (to get grids) --------------------
obs_grid_dummy.dim.x = 0;
obs_grid_dummy.dim.y = 0;
obs_grid_dummy.dim.z = 0;

fprintf('\n==================== KING (time & memory) ====================\n');
mem0 = get_mem_mb(); t0 = tic;

resK = calc_ultrasound_field( ...
source, medium, calc, obs_grid_dummy, ...
'king', calc.dim.method, calc.king.method);

tK = toc(t0); mem1 = get_mem_mb();

fprintf('King time: %.3f s\n', tK);
fprintf('King memory: start %s, end %s, delta %s\n', fmt_mem(mem0), fmt_mem(mem1), fmt_mem(mem1-mem0));

rho_full = resK.king.rho(:);
zK = resK.king.z(:);

[~, iz] = min(abs(zK - z_target));
z_use = zK(iz);
rho_ds = rho_full(1:ds:end);

%% -------------------- DIM --------------------
fprintf('\n==================== DIM (time & memory) =====================\n');
mem2 = get_mem_mb(); t1 = tic;

if strcmpi(calc.dim.method,'rayleigh')
% Rayleigh: compute directly on rho_ds
obs_grid.dim.x = rho_ds(:).';
obs_grid.dim.y = 0;
obs_grid.dim.z = z_use;
obs_grid.dim.block_size = 200000;

resD = calc_ultrasound_field( ...
source, medium, calc, obs_grid, ...
'dim', 'rayleigh', calc.king.method);

elseif strcmpi(calc.dim.method,'asm')
% ASM: compute on native uniform grid (function output), then interp in MAIN
obs_grid.dim.x = 0; % not used by ASM computation
obs_grid.dim.y = 0; % not used by ASM computation
obs_grid.dim.z = z_use;

resD = calc_ultrasound_field( ...
source, medium, calc, obs_grid, ...
'dim', 'asm', calc.king.method);

else
error('calc.dim.method must be ''rayleigh'' or ''asm''.');
end

tD = toc(t1); mem3 = get_mem_mb();

fprintf('DIM time: %.3f s\n', tD);
fprintf('DIM memory: start %s, end %s, delta %s\n', fmt_mem(mem2), fmt_mem(mem3), fmt_mem(mem3-mem2));

%% -------------------- extract fields (use f1 here to match your current code) --------------------
pK = resK.king.p_f1(:, iz);
pK = pK(1:ds:end); % compare on rho_ds

if strcmpi(calc.dim.method,'rayleigh')
pD = squeeze(resD.dim.p_f1(1,:,1)).'; % Nx x 1, already on rho_ds
else
% ASM output: p_f1 is Ny x Nx x Nz on native uniform grid
xA = resD.dim.x(:);
yA = resD.dim.y(:);
zA = resD.dim.z(:);

[~, iy0] = min(abs(yA - 0));
[~, iz0] = min(abs(zA - z_use));

pA_line = squeeze(resD.dim.p_f1(iy0,:,iz0)).'; % numel(xA) x 1

% interpolate to rho_ds ONLY in MAIN (as requested)
pD = interp1(xA, pA_line, rho_ds, 'linear', 0);
end

%% -------------------- metrics --------------------
magK = abs(pK);
magD = abs(pD);

phK = unwrap(angle(pK));
phD = unwrap(angle(pD));

eps0 = 1e-12;
rel_err_log = log10( abs(pD - pK) ./ (abs(pK) + eps0) );

%% -------------------- plots: 1D compare (rho line) --------------------
figure('Name',sprintf('King vs DIM-%s (z≈1 m)', upper(string(calc.dim.method))), ...
'position',[100 100 1200 1200]);

subplot(4,1,1);
plot(rho_ds, 20*log10(magK / medium.pref / sqrt(2)), 'LineWidth',1.5); hold on;
plot(rho_ds, 20*log10(magD / medium.pref / sqrt(2)), '--', 'LineWidth',1.5);
grid on;
xlabel('\rho (m)'); ylabel('SPL (dB)');
title(sprintf('Magnitude @ z = %.6f m', z_use));
legend('King','DIM');

subplot(4,1,2);
plot(rho_ds, magK, 'LineWidth',1.5); hold on;
plot(rho_ds, magD, '--', 'LineWidth',1.5);
grid on;
xlabel('\rho (m)'); ylabel('|p| (Pa)');
title(sprintf('Magnitude @ z = %.6f m', z_use));
legend('King','DIM');

subplot(4,1,3);
plot(rho_ds, phK, 'LineWidth',1.5); hold on;
plot(rho_ds, phD, '--', 'LineWidth',1.5);
grid on;
xlabel('\rho (m)'); ylabel('Phase (rad)');
title('Phase (unwrap)');
legend('King','DIM');

subplot(4,1,4);
plot(rho_ds, rel_err_log, 'LineWidth',1.5);
grid on;
xlabel('\rho (m)');
ylabel('log_{10} relative error');
title('log_{10}(|p_{DIM}-p_{King}| / |p_{King}|)');

fprintf('\n==================== SUMMARY ====================\n');
fprintf('dim_method = %s\n', string(calc.dim.method));
fprintf('z_target = %.3f m, z_use = %.6f m (nearest King grid)\n', z_target, z_use);
fprintf('King: %.3f s, DIM: %.3f s\n', tK, tD);
fprintf('NOTE: memory readings may be NaN if OS API is unavailable.\n');






%% ============================================================
%% ============================================================
%% ============================================================
% 2D FIG #1: xOz (FHT/King only)
% - show |p| (AMP) and SPL on xOz
% - each figure uses a single unified colorbar range (trivial, single axes)
% ============================================================

% choose which frequency field to plot
p_rz = resK.king.p_f1; % Nr x Nz
z_ultra = resK.king.z(:).'; % 1 x Nz
r_ultra = resK.king.rho(:); % Nr x 1

rho_max = calc.fht.rho_max;
idx_r = find(r_ultra <= rho_max, 1, 'last');
if isempty(idx_r); idx_r = numel(r_ultra); end

r_use = r_ultra(1:idx_r); % Nr_use x 1
p_use = p_rz(1:idx_r, :); % Nr_use x Nz

% mirror to xOz (no interpolation)
x_axis = [flipud(-r_use); r_use]; % (2*Nr_use) x 1
p_xz = [flipud(p_use); p_use]; % (2*Nr_use) x Nz

AMP_xz = abs(p_xz);
SPL_xz = 20*log10(AMP_xz / medium.pref / sqrt(2)); % RMS SPL

% ---- AMP (|p|) map ----
figure('Name','FHT: xOz AMP (|p|)');
pcolor(z_ultra, x_axis, AMP_xz);
colormap(MyColor('vik'));
shading interp;
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$|p|$ (Pa)';
set(clb,'Fontsize',20);

xlim([0 z_ultra(end)]); ylim([-rho_max rho_max]);
fontsize(gca,24,'points');
xlabel('$z$ (m)', 'Interpreter','latex','Fontsize',21);
ylabel('$x$ (m)', 'Interpreter','latex','Fontsize',21);
pbaspect([z_ultra(end), 2*rho_max, 1]);
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');

% ---- SPL map ----
figure('Name','FHT: xOz SPL');
pcolor(z_ultra, x_axis, SPL_xz);
colormap(MyColor('vik'));
shading interp;
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$\mathrm{SPL}$ (dB re $p_{\mathrm{ref}}$, RMS)';
set(clb,'Fontsize',20);

xlim([0 z_ultra(end)]); ylim([-rho_max rho_max]);
fontsize(gca,24,'points');
xlabel('$z$ (m)', 'Interpreter','latex','Fontsize',21);
ylabel('$x$ (m)', 'Interpreter','latex','Fontsize',21);
pbaspect([z_ultra(end), 2*rho_max, 1]);
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');

% optional:
% clim([45,75]);


%% ============================================================
% 2D FIG #2: xOy @ z≈1 m (FHT vs DIM) — SAME colorbar per figure
% - Construct polar grid (rho,theta) and do phase extension:
% p2D(r,theta) = p_line(r) * exp(1j*m*theta)
% - Plot SPL and AMP in separate figures
% - In EACH figure, the two subplots share the SAME clim
% ============================================================

% ---- view settings ----
r_boundary = 0.30; % radius to show (m)
theta_fig = 0:0.01:2*pi;

% ---- m ----
if isfield(resK,'source') && isfield(resK.source,'m_used')
m_use = resK.source.m_used;
else
m_use = source.m;
end

% ---- radial samples (NO interp) ----
rho_fig = rho_ds(:).'; % 1 x Nr_fig
pK_full = resK.king.p_f1(:, iz); % Nr_full x 1
pK_line = pK_full(1:ds:end).'; % 1 x Nr_fig
pD_line = pD(:).'; % 1 x Nr_fig (already on rho_ds)

% crop to r_boundary
idx_rb = find(rho_fig <= r_boundary, 1, 'last');
if isempty(idx_rb); idx_rb = 1; end
rho_fig = rho_fig(1:idx_rb);
pK_line = pK_line(1:idx_rb);
pD_line = pD_line(1:idx_rb);

% polar grid -> Cartesian
[TH, R] = meshgrid(theta_fig, rho_fig);
[X, Y] = pol2cart(TH, R);

% phase extension
pK_2D = (pK_line(:) * ones(1, numel(theta_fig))) .* exp(1i*m_use*TH);
pD_2D = (pD_line(:) * ones(1, numel(theta_fig))) .* exp(1i*m_use*TH);

% fields
AMP_K = abs(pK_2D);
AMP_D = abs(pD_2D);

SPL_K = 20*log10(AMP_K / medium.pref / sqrt(2)); % RMS SPL
SPL_D = 20*log10(AMP_D / medium.pref / sqrt(2));

PH_K = angle(pK_2D)/pi;
PH_D = angle(pD_2D)/pi;

% ---- unified color limits (important) ----
spl_lim = [min([SPL_K(:); SPL_D(:)]), max([SPL_K(:); SPL_D(:)])];
amp_lim = [min([AMP_K(:); AMP_D(:)]), max([AMP_K(:); AMP_D(:)])];
ph_lim = [-1 1];

%% ---- FIG: xOy SPL (two subplots, SAME colorbar range) ----
figure('Name',sprintf('xOy SPL @ z=%.6f m, m=%d (FHT vs DIM)', z_use, m_use), 'position',[100 100 1400 650]);

subplot(1,2,1);
pcolor(X, Y, SPL_K);
shading flat
xlim([-1.2*r_boundary 1.2*r_boundary]);
ylim([-1.2*r_boundary 1.2*r_boundary]);
pbaspect([1 1 1])
colormap(MyColor('vik'))
clim(spl_lim);
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$\mathrm{SPL}$ (dB re $p_{\mathrm{ref}}$, RMS)';
set(clb,'Fontsize',18);
set(gca,'position',[0.07 0.12 0.38 0.78]);
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');
xlabel('$x$ (m)','Interpreter','latex','Fontsize',18);
ylabel('$y$ (m)','Interpreter','latex','Fontsize',18);
title('FHT (from King)','Interpreter','latex','Fontsize',20);

subplot(1,2,2);
pcolor(X, Y, SPL_D);
shading flat
xlim([-1.2*r_boundary 1.2*r_boundary]);
ylim([-1.2*r_boundary 1.2*r_boundary]);
pbaspect([1 1 1])
colormap(MyColor('vik'))
clim(spl_lim);
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$\mathrm{SPL}$ (dB re $p_{\mathrm{ref}}$, RMS)';
set(clb,'Fontsize',18);
set(gca,'position',[0.55 0.12 0.38 0.78]);
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');
xlabel('$x$ (m)','Interpreter','latex','Fontsize',18);
ylabel('$y$ (m)','Interpreter','latex','Fontsize',18);
title(sprintf('DIM-%s', upper(string(calc.dim.method))), 'Interpreter','latex','Fontsize',20);

sgtitle(sprintf('$xOy$ SPL @ $z=%.6f$ m, $m=%d$, $r\\le%.2f$ m', z_use, m_use, r_boundary), ...
'Interpreter','latex','Fontsize',20);

%% ---- FIG: xOy AMP (|p|) (two subplots, SAME colorbar range) ----
figure('Name',sprintf('xOy AMP @ z=%.6f m, m=%d (FHT vs DIM)', z_use, m_use), 'position',[100 100 1400 650]);

subplot(1,2,1);
pcolor(X, Y, AMP_K);
shading flat
xlim([-1.2*r_boundary 1.2*r_boundary]);
ylim([-1.2*r_boundary 1.2*r_boundary]);
pbaspect([1 1 1])
colormap(MyColor('vik'))
clim(amp_lim);
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$|p|$ (Pa)';
set(clb,'Fontsize',18);
set(gca,'position',[0.07 0.12 0.38 0.78]);
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');
xlabel('$x$ (m)','Interpreter','latex','Fontsize',18);
ylabel('$y$ (m)','Interpreter','latex','Fontsize',18);
title('FHT (from King)','Interpreter','latex','Fontsize',20);

subplot(1,2,2);
pcolor(X, Y, AMP_D);
shading flat
xlim([-1.2*r_boundary 1.2*r_boundary]);
ylim([-1.2*r_boundary 1.2*r_boundary]);
pbaspect([1 1 1])
colormap(MyColor('vik'))
clim(amp_lim);
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$|p|$ (Pa)';
set(clb,'Fontsize',18);
set(gca,'position',[0.55 0.12 0.38 0.78]);
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');
xlabel('$x$ (m)','Interpreter','latex','Fontsize',18);
ylabel('$y$ (m)','Interpreter','latex','Fontsize',18);
title(sprintf('DIM-%s', upper(string(calc.dim.method))), 'Interpreter','latex','Fontsize',20);

sgtitle(sprintf('$xOy$ AMP ($|p|$) @ $z=%.6f$ m, $m=%d$, $r\\le%.2f$ m', z_use, m_use, r_boundary), ...
'Interpreter','latex','Fontsize',20);

%% ---- FIG: xOy Phase/pi (two subplots, SAME colorbar range) ----
figure('Name',sprintf('xOy Phase/pi @ z=%.6f m, m=%d (FHT vs DIM)', z_use, m_use), 'position',[100 100 1400 650]);

subplot(1,2,1);
pcolor(X, Y, PH_K);
shading flat
xlim([-1.2*r_boundary 1.2*r_boundary]);
ylim([-1.2*r_boundary 1.2*r_boundary]);
pbaspect([1 1 1])
colormap('hsv')
clim(ph_lim);
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$\angle p/\pi$';
set(clb,'Fontsize',18);
set(gca,'position',[0.07 0.12 0.38 0.78]);
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');
xlabel('$x$ (m)','Interpreter','latex','Fontsize',18);
ylabel('$y$ (m)','Interpreter','latex','Fontsize',18);
title('FHT (from King)','Interpreter','latex','Fontsize',20);

subplot(1,2,2);
pcolor(X, Y, PH_D);
shading flat
xlim([-1.2*r_boundary 1.2*r_boundary]);
ylim([-1.2*r_boundary 1.2*r_boundary]);
pbaspect([1 1 1])
colormap('hsv')
clim(ph_lim);
clb = colorbar;
clb.Title.Interpreter = 'latex';
clb.Title.String = '$\angle p/\pi$';
set(clb,'Fontsize',18);
set(gca,'position',[0.55 0.12 0.38 0.78]);
set(gca,'linewidth',2);
set(gca,'TickLabelInterpreter','latex');
xlabel('$x$ (m)','Interpreter','latex','Fontsize',18);
ylabel('$y$ (m)','Interpreter','latex','Fontsize',18);
title(sprintf('DIM-%s', upper(string(calc.dim.method))), 'Interpreter','latex','Fontsize',20);

sgtitle(sprintf('$xOy$ Phase$/\\pi$ @ $z=%.6f$ m, $m=%d$, $r\\le%.2f$ m', z_use, m_use, r_boundary), ...
'Interpreter','latex','Fontsize',20);








%% ==================== local functions ====================
function mem_mb = local_get_mem_mb()
% Returns MATLAB process memory usage in MB if possible; otherwise NaN.
mem_mb = NaN;

% Prefer Windows: use memory() for MATLAB process peak/usage
if ispc
try
m = memory();
mem_mb = double(m.MemUsedMATLAB) / (1024^2);
return;
catch
end
end

% macOS/Linux: Java heap usage (trend only, not full RSS)
try
rt = java.lang.Runtime.getRuntime();
used = double(rt.totalMemory() - rt.freeMemory()); % bytes
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

