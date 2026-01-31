function result = calc_ultrasound_velocity_field(source, medium, calc, method)
% =========================================================================
% CALC_ULTRASOUND_VELOCITY_FIELD
% - Compute ultrasonic particle-velocity field v = (v_rho, v_phi, v_z)
%   at f1/f2 on the SAME (rho,z) grid as your King/FHT pressure solver.
%
% METHODS
%   method = 'king'  : King/FHT velocity (v_rho, v_phi, v_z)
%   method = 'dim'   : Direct-integration method (DIM) ONLY v_z
%   method = 'both'  : compute both King (all comps) + DIM (v_z)
%
% NOTES (DIM-vz)
%   - DIM implementation here is a "direct surface integral" on the source
%     plane z=0 using a polar grid (rho_s, phi_s).
%   - It computes ONLY v_z on (rho,z) (core field, without exp(i*m*phi)).
%   - Kernel used: d/dz [exp(i*k*R)/R]  where R is distance to source points.
%   - A scale factor "calc.dim.coef" is provided (default 1/(2*pi)) to help
%     align conventions with your King/FHT normalization if needed.
%
% OUTPUT (struct)
%   result.king (if computed): v_rho_f1/f2, v_phi_f1/f2, v_z_f1/f2
%   result.dim  (if computed): v_z_f1/f2
% =========================================================================

if nargin < 4 || isempty(method), method = 'king'; end
method = lower(strtrim(string(method)));

do_king = any(method == ["king","both"]);
do_dim  = any(method == ["dim","both"]);

if ~do_king && ~do_dim
    error('calc_ultrasound_velocity_field:BadMethod', ...
        'method must be ''king'', ''dim'', or ''both''.');
end

% -------------------- defaults from calc.king --------------------
if nargin < 3 || isempty(calc) || ~isstruct(calc), calc = struct(); end
if ~isfield(calc,'king') || isempty(calc.king), calc.king = struct(); end

if ~isfield(calc.king,'eps_phase') || isempty(calc.king.eps_phase)
    calc.king.eps_phase = 1e-3;
end
if ~isfield(calc.king,'kz_min') || isempty(calc.king.kz_min)
    calc.king.kz_min = 1e-12;
end
if ~isfield(calc.king,'gspec_method') || isempty(calc.king.gspec_method)
    calc.king.gspec_method = 'analytic';
end
king_gspec_method = lower(strtrim(string(calc.king.gspec_method)));

% -------------------- defaults from calc.dim (new) --------------------
if ~isfield(calc,'dim') || isempty(calc.dim), calc.dim = struct(); end

% NEW: source discretization for DIM (default cart)
if ~isfield(calc.dim,'src_discretization') || isempty(calc.dim.src_discretization)
    calc.dim.src_discretization = 'cart';    % 'cart' (default) | 'polar'
end
calc.dim.src_discretization = lower(strtrim(string(calc.dim.src_discretization)));
if ~any(calc.dim.src_discretization == ["cart","polar"])
    error('calc_ultrasound_velocity_field:BadSrcDiscretization', ...
        'calc.dim.src_discretization must be ''cart'' or ''polar''.');
end
if ~isfield(calc.dim,'Nphi') || isempty(calc.dim.Nphi)
    calc.dim.Nphi = 256;
end
if ~isfield(calc.dim,'block_size') || isempty(calc.dim.block_size)
    calc.dim.block_size = 20000;
end
if ~isfield(calc.dim,'coef') || isempty(calc.dim.coef)
    calc.dim.coef = 1/(2*pi);
end

% -------------------- build consistent source + grids --------------------
[source, fht, ~] = make_source_velocity(source, medium, calc);

result = struct();
result.source = source;
result.medium = source.medium;
result.calc   = calc;
result.fht    = fht;

% =========================================================================
% Common grids
% =========================================================================
rho = fht.xh(:);            % (Nrho x 1) physical rho grid
z   = fht.z_ultra(:).';     % (1 x Nz)  (typically >=0)
Nrho = numel(rho);
Nz   = numel(z);

m = source.m_used;

% =========================================================================
% KING / FHT velocity field (existing)
% =========================================================================
if do_king
    % ---- Hankel-domain grids ----
    kr = fht.xH(:);             % (Nkr x 1)
    Nkr = numel(kr);

    % ---- Vs at f1/f2 (H_m{ v_n(rho) } on xH) ----
    Vs1 = m_FHT(source.fht.vs1_on_xh_v, ...
        fht.N_FHT, 1, fht.Nh_v, fht.NH_v, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, m);

    Vs2 = m_FHT(source.fht.vs2_on_xh_v, ...
        fht.N_FHT, 1, fht.Nh_v, fht.NH_v, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, m);

    Vs1 = Vs1(:);
    Vs2 = Vs2(:);

    % ---- kz branch (imag(kz) >= 0), and |z| / sgn(z) ----
    Zabs = abs(z);
    sgnZ = sign(z); sgnZ(sgnZ==0) = 1;

    KZ1 = local_kz_branch(source.k1, kr, Nz);
    KZ2 = local_kz_branch(source.k2, kr, Nz);

    KZ1s = local_kz_floor(KZ1, calc.king.kz_min);
    KZ2s = local_kz_floor(KZ2, calc.king.kz_min);

    E1 = exp(1j * KZ1 .* (ones(Nkr,1)*Zabs));
    E2 = exp(1j * KZ2 .* (ones(Nkr,1)*Zabs));

    Scom1 = (Vs1 * ones(1,Nz)) .* E1 ./ KZ1s;
    Scom2 = (Vs2 * ones(1,Nz)) .* E2 ./ KZ2s;

    Sz1   = (Vs1 * ones(1,Nz)) .* E1;
    Sz2   = (Vs2 * ones(1,Nz)) .* E2;

    Sr1   = Scom1 .* (kr * ones(1,Nz));
    Sr2   = Scom2 .* (kr * ones(1,Nz));

    % v_z
    Vz1 = m_FHT(Sz1, fht.N_FHT, Nz, fht.NH, fht.Nh, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, m);
    Vz2 = m_FHT(Sz2, fht.N_FHT, Nz, fht.NH, fht.Nh, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, m);

    Vz1 = Vz1 .* (ones(Nrho,1) * sgnZ);
    Vz2 = Vz2 .* (ones(Nrho,1) * sgnZ);

    % v_phi core
    Vphi_core1 = m_FHT(Scom1, fht.N_FHT, Nz, fht.NH, fht.Nh, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, m);
    Vphi_core2 = m_FHT(Scom2, fht.N_FHT, Nz, fht.NH, fht.Nh, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, m);

    if m == 0
        Vphi1 = complex(zeros(Nrho, Nz));
        Vphi2 = complex(zeros(Nrho, Nz));
    else
        rho_safe = rho;
        rho_safe(rho_safe==0) = min(rho_safe(rho_safe>0)) * 1e-6 + 1e-12;
        Vphi1 = (m ./ rho_safe) * ones(1,Nz) .* Vphi_core1;
        Vphi2 = (m ./ rho_safe) * ones(1,Nz) .* Vphi_core2;
    end

    % v_rho core via (m±1)
    Inv_m1_f1 = m_FHT(Sr1, fht.N_FHT, Nz, fht.NH, fht.Nh, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, m-1);
    Inv_p1_f1 = m_FHT(Sr1, fht.N_FHT, Nz, fht.NH, fht.Nh, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, m+1);

    Inv_m1_f2 = m_FHT(Sr2, fht.N_FHT, Nz, fht.NH, fht.Nh, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, m-1);
    Inv_p1_f2 = m_FHT(Sr2, fht.N_FHT, Nz, fht.NH, fht.Nh, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, m+1);

    Vr1 = -(1i/2) * (Inv_m1_f1 - Inv_p1_f1);
    Vr2 = -(1i/2) * (Inv_m1_f2 - Inv_p1_f2);

    % pack
    result.king = struct();
    result.king.method            = "King-FHT velocity (Eq.13)";
    result.king.green_spec_method = king_gspec_method;
    result.king.m                 = m;
    result.king.rho               = rho;
    result.king.z                 = z;

    result.king.v_rho_f1 = Vr1;
    result.king.v_phi_f1 = Vphi1;
    result.king.v_z_f1   = Vz1;

    result.king.v_rho_f2 = Vr2;
    result.king.v_phi_f2 = Vphi2;
    result.king.v_z_f2   = Vz2;

    result.king.f1 = source.f1;
    result.king.f2 = source.f2;
end

% =========================================================================
% DIM (Rayleigh-style): velocity on ONE plane z=z_use
% Using the SAME convention as your DIM-pressure in calc_ultrasound_field:
%   p = -2j*rho0*w*phi,   phi(r)=∬ G(r,rs) q(rs)dS,  G=exp(ikR)/(4*pi*R), q=vn*dS
% Then:
%   v = (1/(j*w*rho0))∇p = -2∇phi
%
% We compute v_x, v_y, v_z by integrating ∂G/∂x, ∂G/∂y, ∂G/∂z, then convert to
% v_rho, v_phi at each (x,y).
% =========================================================================
if do_dim
    % ---------- forbid ASM with polar source discretization (keep consistent) ----------
    if isfield(calc,'dim') && isfield(calc.dim,'src_discretization')
        if strcmpi(calc.dim.src_discretization, "polar")
            % Your DIM velocity block is Rayleigh-style point summation;
            % it's OK for polar. Only ASM must be Cartesian.
            % (If you later add ASM here, forbid it.)
        end
    end

    % ---------- require z_use ----------
    if ~isfield(calc,'dim') || ~isfield(calc.dim,'z_use') || isempty(calc.dim.z_use)
        error('calc_ultrasound_velocity_field:DIMNeedZuse', ...
            'Set calc.dim.z_use in main script.');
    end
    z_pick = calc.dim.z_use;

    % ---------- require obs grid (x,y,z) ----------
    if ~isfield(calc.dim,'obs_grid') || ~isstruct(calc.dim.obs_grid)
        error('calc_ultrasound_velocity_field:DIMNeedObsGrid', ...
            'Set calc.dim.obs_grid = obs_grid.dim in main script (fields x,y,z).');
    end
    gd = calc.dim.obs_grid;

    if ~isfield(gd,'x') || isempty(gd.x), error('DIM obs_grid needs x'); end
    if ~isfield(gd,'y') || isempty(gd.y), error('DIM obs_grid needs y'); end
    if ~isfield(gd,'z') || isempty(gd.z), error('DIM obs_grid needs z'); end

    x_obs = gd.x(:).';     % 1 x Nx
    y_obs = gd.y(:).';     % 1 x Ny
    z_obs = gd.z(:).';     % must be scalar (one plane)

    if numel(z_obs) ~= 1
        error('calc_ultrasound_velocity_field:DIMPlaneOnly', ...
            'This DIM velocity block is plane-only: obs_grid.z must be scalar.');
    end

    % ---------- map z_use to nearest FHT z grid ----------
    [~, iz_dim] = min(abs(z(:) - z_pick));
    z_use_dim = z(iz_dim);

    % ---------- defaults: blocks ----------
    if ~isfield(calc.dim,'block_size') || isempty(calc.dim.block_size)
        calc.dim.block_size = 20000;
    end
    if ~isfield(calc.dim,'src_block_size') || isempty(calc.dim.src_block_size)
        calc.dim.src_block_size = 5000;
    end
    blk     = max(1000, round(calc.dim.block_size));
    src_blk = max(500,  round(calc.dim.src_block_size));

    % ---------- source points ----------
    if ~isfield(source,'dim') || ~isfield(source.dim,'Xs') || isempty(source.dim.Xs)
        error('calc_ultrasound_velocity_field:DIMNeedSourceDim', ...
            'source.dim.* not found. Ensure make_source_velocity builds DIM source points.');
    end
    Xs = source.dim.Xs(:);
    Ys = source.dim.Ys(:);
    Zs = source.dim.Zs(:);

    vn1 = source.dim.Vn_pts_f1(:);
    vn2 = source.dim.Vn_pts_f2(:);

    % q = vn * dS (MUST be consistent with pressure DIM)
    src_mode = "cart";
    if isfield(calc,'dim') && isfield(calc.dim,'src_discretization')
        src_mode = string(calc.dim.src_discretization);
    end

    if src_mode == "polar"
        % polar: require per-point area weights
        if ~isfield(source.dim,'dA_pts') || isempty(source.dim.dA_pts)
            error('calc_ultrasound_velocity_field:PolarNeedsDApts', ...
                'Polar source discretization requires source.dim.dA_pts (per-point area).');
        end
        dA_w = source.dim.dA_pts(:);
        if numel(dA_w) ~= numel(vn1)
            error('calc_ultrasound_velocity_field:BadDA', ...
                'Polar: numel(dA_pts) must equal numel(Vn_pts).');
        end
        q1 = vn1 .* dA_w;
        q2 = vn2 .* dA_w;

    else
        % cart: allow scalar dA or dA_pts (either works)
        if isfield(source.dim,'dA_pts') && ~isempty(source.dim.dA_pts)
            dA_w = source.dim.dA_pts(:);
            if numel(dA_w) ~= numel(vn1)
                error('calc_ultrasound_velocity_field:BadDA', ...
                    'Cart: numel(dA_pts) must equal numel(Vn_pts) if provided.');
            end
            q1 = vn1 .* dA_w;
            q2 = vn2 .* dA_w;
        else
            if ~isfield(source.dim,'dA') || isempty(source.dim.dA)
                error('calc_ultrasound_velocity_field:BadDA', ...
                    'Cart: need source.dim.dA (scalar) or dA_pts.');
            end
            q1 = vn1 * source.dim.dA;
            q2 = vn2 * source.dim.dA;
        end
    end

    Nsrc = numel(Xs);

    % ---------- observation plane mesh ----------
    [X2, Y2] = meshgrid(x_obs, y_obs);     % Ny x Nx
    Ny = size(X2,1); Nx = size(X2,2);
    xo = X2(:); yo = Y2(:);
    zo = z_use_dim;
    Nobs = numel(xo);

    % ---------- allocate vectors ----------
    vx1_vec = complex(zeros(Nobs,1));
    vy1_vec = complex(zeros(Nobs,1));
    vz1_vec = complex(zeros(Nobs,1));

    vx2_vec = complex(zeros(Nobs,1));
    vy2_vec = complex(zeros(Nobs,1));
    vz2_vec = complex(zeros(Nobs,1));

    % ---------- wavenumbers ----------
    k1 = source.k1;
    k2 = source.k2;

    % ---------- fixed coefficient (DO NOT expose) ----------
    % From v = -2 ∇phi given your pressure convention.
    C_DIM = -2;

    % ---------- integrate in blocks ----------
    for s = 1:blk:Nobs
        id = s:min(s+blk-1, Nobs);
        nb = numel(id);

        vx1_blk = complex(zeros(nb,1));
        vy1_blk = complex(zeros(nb,1));
        vz1_blk = complex(zeros(nb,1));

        vx2_blk = complex(zeros(nb,1));
        vy2_blk = complex(zeros(nb,1));
        vz2_blk = complex(zeros(nb,1));

        for js = 1:src_blk:Nsrc
            jd = js:min(js+src_blk-1, Nsrc);

            dX = xo(id) - Xs(jd).';
            dY = yo(id) - Ys(jd).';
            dZ = (zo)   - Zs(jd).';

            R  = sqrt(dX.^2 + dY.^2 + dZ.^2);
            R(R < 1e-9) = 1e-9;

            % ∂/∂coord [exp(ikR)/(4πR)] = (dCoord)*exp(ikR)/(4π)*(ik/R^2 - 1/R^3)
            common1 = exp(1j*k1.*R) ./ (4*pi) .* (1j*k1./(R.^2) - 1./(R.^3));
            common2 = exp(1j*k2.*R) ./ (4*pi) .* (1j*k2./(R.^2) - 1./(R.^3));

            Gx1 = dX .* common1;  Gy1 = dY .* common1;  Gz1 = dZ .* common1;
            Gx2 = dX .* common2;  Gy2 = dY .* common2;  Gz2 = dZ .* common2;

            vx1_blk = vx1_blk + Gx1 * q1(jd);
            vy1_blk = vy1_blk + Gy1 * q1(jd);
            vz1_blk = vz1_blk + Gz1 * q1(jd);

            vx2_blk = vx2_blk + Gx2 * q2(jd);
            vy2_blk = vy2_blk + Gy2 * q2(jd);
            vz2_blk = vz2_blk + Gz2 * q2(jd);
        end

        vx1_vec(id) = C_DIM * vx1_blk;
        vy1_vec(id) = C_DIM * vy1_blk;
        vz1_vec(id) = C_DIM * vz1_blk;

        vx2_vec(id) = C_DIM * vx2_blk;
        vy2_vec(id) = C_DIM * vy2_blk;
        vz2_vec(id) = C_DIM * vz2_blk;
    end

    % ---------- reshape back to plane ----------
    VX1 = reshape(vx1_vec, Ny, Nx);
    VY1 = reshape(vy1_vec, Ny, Nx);
    VZ1 = reshape(vz1_vec, Ny, Nx);

    VX2 = reshape(vx2_vec, Ny, Nx);
    VY2 = reshape(vy2_vec, Ny, Nx);
    VZ2 = reshape(vz2_vec, Ny, Nx);

    % ---------- cylindrical transform ----------
    PHI = atan2(Y2, X2);
    cP = cos(PHI); sP = sin(PHI);

    VR1   = VX1 .* cP + VY1 .* sP;
    VPHI1 = -VX1 .* sP + VY1 .* cP;

    VR2   = VX2 .* cP + VY2 .* sP;
    VPHI2 = -VX2 .* sP + VY2 .* cP;

    % ---------- pack result.dim ----------
    result.dim = struct();
    result.dim.method          = "DIM-Rayleigh velocity (plane-only)";
    result.dim.z_use           = z_use_dim;
    result.dim.iz_dim          = iz_dim;
    result.dim.x               = x_obs;
    result.dim.y               = y_obs;
    result.dim.X               = X2;
    result.dim.Y               = Y2;

    result.dim.v_x_f1          = VX1;
    result.dim.v_y_f1          = VY1;
    result.dim.v_z_f1          = VZ1;
    result.dim.v_rho_f1        = VR1;
    result.dim.v_phi_f1        = VPHI1;

    result.dim.v_x_f2          = VX2;
    result.dim.v_y_f2          = VY2;
    result.dim.v_z_f2          = VZ2;
    result.dim.v_rho_f2        = VR2;
    result.dim.v_phi_f2        = VPHI2;

    result.dim.block_size      = blk;
    result.dim.src_block_size  = src_blk;
    result.dim.f1              = source.f1;
    result.dim.f2              = source.f2;
end

end

% ======================================================================
% helpers
% ======================================================================
function KZ = local_kz_branch(k, kr_col, Nz)
% kz = sqrt(k^2 - kr^2), with branch chosen so imag(kz) >= 0
kz = sqrt(k.^2 - kr_col.^2);           % (Nkr x 1)
idx = imag(kz) < 0;
kz(idx) = -kz(idx);
KZ = kz * ones(1, Nz);                 % (Nkr x Nz)
end

function KZs = local_kz_floor(KZ, kz_min)
% floor small |kz| to avoid division blow-up
KZs = KZ;
mask = abs(KZs) < kz_min;
if any(mask(:))
    KZs(mask) = kz_min .* exp(1j * angle(KZs(mask)));
end
end

function w = local_trapz_weights_1d(x)
% Trapezoidal weights for possibly-nonuniform grid x (monotone increasing).
% such that sum(f .* w) ~ ∫ f(x) dx
x = x(:);
n = numel(x);
if n < 2
    w = 0;
    return;
end
dx = diff(x);
w = zeros(n,1);
w(1)   = dx(1)/2;
w(end) = dx(end)/2;
if n > 2
    w(2:end-1) = (dx(1:end-1) + dx(2:end))/2;
end
end
