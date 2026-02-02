%% ============================================================
%  MAIN: DIM-only Bessel beam (plane z≈1 m) velocity components
%  - Calls:
%       make_source_velocity()
%       calc_ultrasound_velocity_field()   (method = 'dim' only)
%  - Plots: v_x, v_y, v_z at f2 on z=1 m plane
%       each component: (magnitude, phase) as 1x2 subplots
%
%  NOTE:
%  This script generates a "Bessel-like" beam by applying an axicon-like
%  radial phase on a circular piston:
%       v_n(r) = v0 * exp(-j * k_ref * sin(theta) * r),  r<=a
%  which synthesizes a conical wave component -> truncated Bessel beam.
% ============================================================

clear; clc; close all;

%% -------------------- medium --------------------
medium.c0   = 343;
medium.rho0 = 1.21;
medium.beta = 1.2;
medium.pref = 2e-5;
medium.use_absorp = false;

%% -------------------- Bessel-like source (axicon phase) --------------------
source = struct();
source.profile = 'Custom';
source.a  = 0.10;      % aperture radius (m)
source.v0 = 0.172;     % normal velocity amplitude (m/s)

% PAL-style frequencies (you can change)
source.f1 = 42e3;      % carrier
source.fa = 1e3;       % audio diff
source.f2 = source.f1 + source.fa;

% Choose cone angle (controls kr = k*sin(theta))
theta_deg = 10;                 % adjust (e.g., 5~20 deg)
theta = deg2rad(theta_deg);

% Use f2 as reference for the axicon phase (consistent with your grids)
kref = 2*pi*source.f2 / medium.c0;

% Custom radial velocity at f1/f2 (here: same cone angle, different k)
k1 = 2*pi*source.f1/medium.c0;
k2 = 2*pi*source.f2/medium.c0;

source.custom_vrho_handle_1 = @(rho) source.v0 .* exp(-1i * (k1*sin(theta)) .* rho);
source.custom_vrho_handle_2 = @(rho) source.v0 .* exp(-1i * (k2*sin(theta)) .* rho);

%% -------------------- calc settings --------------------
calc = struct();

% ---- FHT settings are not used in DIM-only compute, but make_source_velocity needs them
calc.fht.N_FHT = 32768;
calc.fht.rho_max = 2;
calc.fht.Nh_scale = 1.2;
calc.fht.NH_scale = 4;
calc.fht.Nh_v_scale = 1.1;
calc.fht.zu_max = 15;
calc.fht.za_max = 4;

% ---- DIM discretization (source plane)
calc.dim.use_freq = 'f2';          % dx based on lambda(f2)
calc.dim.dis_coe  = 16;            % dx = lambda/16
calc.dim.src_discretization = 'polar';   % 'cart' or 'polar' (polar is efficient for disk)
calc.dim.margin   = 1;
calc.dim.center_x = 0;
calc.dim.center_y = 0;
calc.dim.z0       = 0;

% ---- DIM observation plane near z=1m
z_plot = 1.0;
calc.dim.z_use = z_plot;

% Observation grid (plane)
xmax = 0.30;          % view window half-width (m)
dx_obs = 0.003;       % observation spacing (m)
x_obs = -xmax:dx_obs:xmax;
y_obs = -xmax:dx_obs:xmax;

calc.dim.obs_grid = struct();
calc.dim.obs_grid.x = x_obs;
calc.dim.obs_grid.y = y_obs;
calc.dim.obs_grid.z = z_plot;      % must be scalar for this DIM velocity implementation

% Block sizes (tune if needed)
calc.dim.block_size     = 20000;   % observation block
calc.dim.src_block_size = 5000;    % source block

%% -------------------- compute (DIM ONLY) --------------------
result = calc_ultrasound_velocity_field(source, medium, calc, 'dim');

% Use f2 component (ultrasonic sideband)
X  = result.dim.X;
Y  = result.dim.Y;

VX = result.dim.v_x_f2;
VY = result.dim.v_y_f2;
VZ = result.dim.v_z_f2;

%% -------------------- plot helper --------------------
plot_comp = @(V, nameStr) local_plot_mag_phase(X, Y, V, nameStr, z_plot, source.f2, theta_deg);

plot_comp(VX, 'v_x');
plot_comp(VY, 'v_y');
plot_comp(VZ, 'v_z');

%% ===================== local function =====================
function local_plot_mag_phase(X, Y, V, compName, z_plot, f_hz, theta_deg)
    figure('Color','w'); 
    set(gcf,'Name',sprintf('%s @ z=%.3f m, f=%.1f kHz, theta=%.1f°', ...
        compName, z_plot, f_hz/1e3, theta_deg));

    % magnitude
    subplot(1,2,1);
    imagesc(X(1,:), Y(:,1), abs(V));
    axis image; set(gca,'YDir','normal');
    colorbar;
    title(sprintf('|%s| (m/s)', compName), 'Interpreter','none');
    xlabel('x (m)'); ylabel('y (m)');

    % phase
    subplot(1,2,2);
    imagesc(X(1,:), Y(:,1), angle(V));
    axis image; set(gca,'YDir','normal');
    colorbar;
    title(sprintf('phase(%s) (rad)', compName), 'Interpreter','none');
    xlabel('x (m)'); ylabel('y (m)');
end
