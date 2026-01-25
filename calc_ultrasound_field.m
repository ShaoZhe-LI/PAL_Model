function result = calc_ultrasound_field(source, medium, calc, grid, method)
% =========================================================================
% CALC_ULTRASOUND_FIELD
% - ULTRASOUND fields at f1/f2 using:
%   (A) King/FHT axisymmetric
%       Green spectrum G(k_rho,z): 'analytic' (default) | 'transform'
%   (B) DIM:
%       'rayleigh' (default, memory-lean block integration)
%       'asm'
%
% DIM-Rayleigh is implemented in a "Cal_Ultra_Field_3D_block"-style:
%   - DO NOT build 3D meshgrid / observationPoints
%   - loop over z slices, build only 2D (x,y) plane per slice
%   - block over observation points in that slice
%   - distance matrix is (block_size x Nsrc) only
%
% INPUT
%   source, medium, calc : passed to make_source_velocity()
%   grid.dim : required if method includes DIM
%       grid.dim.x, grid.dim.y, grid.dim.z   (vectors)
%       grid.dim.block_size (optional, default 200000)
%
%   method: 'king' | 'dim' | 'both'          (default 'both')
%
%   NOTE (UPDATED):
%     dim_method and king_gspec_method are now read from calc:
%       - calc.dim.method        : 'rayleigh' | 'asm'       (default 'rayleigh')
%       - calc.king.gspec_method : 'analytic' | 'transform' (default 'analytic')
%
% OUTPUT
%   result.king.* and/or result.dim.*
% =========================================================================

% -------------------- defaults (method still as input) --------------------
if nargin < 5 || isempty(method)
    method = 'both';
end
method = lower(strtrim(string(method)));

% -------------------- defaults from calc --------------------
if nargin < 3 || isempty(calc) || ~isstruct(calc)
    calc = struct();
end

% dim method from calc.dim.method
if ~isfield(calc,'dim') || isempty(calc.dim); calc.dim = struct(); end
if ~isfield(calc.dim,'method') || isempty(calc.dim.method)
    calc.dim.method = 'rayleigh';
end
dim_method = lower(strtrim(string(calc.dim.method)));

% king gspec method from calc.king.gspec_method
if ~isfield(calc,'king') || isempty(calc.king); calc.king = struct(); end
if ~isfield(calc.king,'gspec_method') || isempty(calc.king.gspec_method)
    calc.king.gspec_method = 'analytic';
end
king_gspec_method = lower(strtrim(string(calc.king.gspec_method)));

do_king = any(method == ["king","both"]);
do_dim  = any(method == ["dim","both"]);

% -------------------- build consistent source + grids --------------------
[source, fht, dim_src] = make_source_velocity(source, medium, calc);

result = struct();
result.source  = source;
result.medium  = source.medium;
result.calc    = calc;
result.fht     = fht;
result.dim_src = dim_src;

% King analytic spectrum settings
if ~isfield(calc,'king') || isempty(calc.king); calc.king = struct(); end
% threshold on |kz*z| for Taylor patch (dimensionless)
if ~isfield(calc.king,'eps_phase') || isempty(calc.king.eps_phase)
    calc.king.eps_phase = 1e-3;
end
% floor for |kz| to avoid division blow-up
if ~isfield(calc.king,'kz_min') || isempty(calc.king.kz_min)
    calc.king.kz_min = 1e-12;
end

% ---- DIM parallel defaults ----
if ~isfield(calc,'dim') || isempty(calc.dim); calc.dim = struct(); end
if ~isfield(calc.dim,'use_parallel') || isempty(calc.dim.use_parallel)
    calc.dim.use_parallel = true;
end
if ~isfield(calc.dim,'num_workers') || isempty(calc.dim.num_workers)
    calc.dim.num_workers = 4;   % default 4 workers
end

% ASM padding settings
if ~isfield(calc,'asm') || isempty(calc.asm); calc.asm = struct(); end
if ~isfield(calc.asm,'pad_factor') || isempty(calc.asm.pad_factor)
    calc.asm.pad_factor = 8;        % 8 or 16
end
if ~isfield(calc.asm,'kzz_eps') || isempty(calc.asm.kzz_eps)
    calc.asm.kzz_eps = 1e-12;
end

% =========================================================================
% (A) KING / FHT
% =========================================================================
if do_king
    rho_k = fht.xh;          % (Nr x 1)
    z_k   = fht.z_ultra;     % (1 x Nz)
    Nz    = numel(z_k);

    Vs1 = m_FHT(source.fht.vs1_on_xh_v, ...
        fht.N_FHT, 1, fht.Nh_v, fht.NH_v, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, source.m_used);

    Vs2 = m_FHT(source.fht.vs2_on_xh_v, ...
        fht.N_FHT, 1, fht.Nh_v, fht.NH_v, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, source.m_used);

    kr   = fht.xH;                 % k_rho samples (Nr x 1)
    zrow = reshape(z_k, 1, []);    % (1 x Nz)

    switch king_gspec_method
        case "analytic"
            G1 = local_green_spec_analytic(source.k1, kr, zrow, calc.king.eps_phase, calc.king.kz_min);
            G2 = local_green_spec_analytic(source.k2, kr, zrow, calc.king.eps_phase, calc.king.kz_min);

        case "transform"
            Nr = numel(rho_k);

            G1 = complex(zeros(Nr, Nz));
            G2 = complex(zeros(Nr, Nz));

            for iz = 1:Nz
                z0 = z_k(iz);
                r  = hypot(rho_k, z0);
                r(r < 1e-9) = 1e-9;

                g1 = exp(1j * source.k1 .* r) ./ (4*pi*r);
                g2 = exp(1j * source.k2 .* r) ./ (4*pi*r);

                G1(:,iz) = m_FHT(g1, fht.N_FHT, 1, fht.Nh, fht.NH, ...
                    fht.a_solve, fht.x0, fht.x1, fht.k0, 0);

                G2(:,iz) = m_FHT(g2, fht.N_FHT, 1, fht.Nh, fht.NH, ...
                    fht.a_solve, fht.x0, fht.x1, fht.k0, 0);
            end

        otherwise
            error('calc_ultrasound_field:BadKingGspecMethod', ...
                'calc.king.gspec_method must be ''analytic'' or ''transform''.');
    end

    F1   = G1 .* Vs1;
    phi1 = -4*pi * m_FHT(F1, fht.N_FHT, Nz, fht.NH, fht.Nh, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, source.m_used);
    p1 = 1j * source.medium.rho0 * source.medium.c0 * real(source.k1) .* phi1;

    F2   = G2 .* Vs2;
    phi2 = -4*pi * m_FHT(F2, fht.N_FHT, Nz, fht.NH, fht.Nh, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, source.m_used);
    p2 = 1j * source.medium.rho0 * source.medium.c0 * real(source.k2) .* phi2;

    result.king = struct();
    result.king.method            = "King-FHT";
    result.king.green_spec_method = king_gspec_method;
    result.king.eps_phase         = calc.king.eps_phase;
    result.king.kz_min            = calc.king.kz_min;
    result.king.rho               = rho_k;
    result.king.z                 = z_k;
    result.king.p_f1              = p1;
    result.king.p_f2              = p2;
    result.king.f1                = source.f1;
    result.king.f2                = source.f2;
    result.king.G1                = G1;
    result.king.Vs1               = Vs1;
end

% =========================================================================
% (B) DIM
% =========================================================================
if do_dim
    if nargin < 4 || isempty(grid) || ~isstruct(grid) || ~isfield(grid,'dim')
        error('calc_ultrasound_field:MissingGrid', 'Using DIM requires grid.dim.');
    end

    gd = grid.dim;
    [x_obs, y_obs, z_obs] = local_parse_dim_grid(gd);
    Nx = numel(x_obs); Ny = numel(y_obs); Nz = numel(z_obs); %#ok<NASGU>

    % --- forbid ASM with polar source discretization ---
    if strcmpi(dim_method, "asm")
        if isfield(source,'dim') && isfield(source.dim,'src_discretization')
            if strcmpi(string(source.dim.src_discretization), "polar")
                error('calc_ultrasound_field:ASMNeedsCartesian', ...
                    'DIM-ASM requires Cartesian (rectangular) source sampling. Set calc.dim.src_discretization = ''cart''.');
            end
        end
    end

    switch dim_method
        case "rayleigh"
            % --- memory-lean block integration ---
            blk = 200000;
            if isfield(gd,'block_size') && ~isempty(gd.block_size)
                blk = gd.block_size;
            end

            % --- source-block size (like blk) ---
            src_blk = 5000;   % default
            if isfield(gd,'src_block_size') && ~isempty(gd.src_block_size)
                src_blk = gd.src_block_size;
            end

            % source points (vector)
            Xs = source.dim.Xs(:);
            Ys = source.dim.Ys(:);
            Zs = source.dim.Zs(:);

            vn1 = source.dim.Vn_pts_f1(:);
            vn2 = source.dim.Vn_pts_f2(:);

            % area weights: scalar (cart legacy) or per-point (polar)
            if isfield(source.dim,'dA_pts') && ~isempty(source.dim.dA_pts)
                dA_w = source.dim.dA_pts;
            else
                dA_w = source.dim.dA;   % legacy scalar
            end
            dA_w = dA_w(:);

            if isscalar(dA_w)
                q1 = vn1 * dA_w;
                q2 = vn2 * dA_w;
            else
                if numel(dA_w) ~= numel(vn1)
                    error('calc_ultrasound_field:BadDA', 'dA_pts size mismatch with Vn_pts.');
                end
                q1 = vn1 .* dA_w;
                q2 = vn2 .* dA_w;
            end

            Nsrc = numel(Xs);

            % Allocate output (Ny x Nx x Nz)
            P1 = complex(zeros(Ny, Nx, Nz));
            P2 = complex(zeros(Ny, Nx, Nz));

            % 2D plane coordinates (Ny x Nx)
            [X2, Y2] = meshgrid(x_obs, y_obs);
            xo = X2(:); yo = Y2(:);
            nP = numel(xo);

            % ===================== parallel policy (robust) =====================
            use_par = logical(calc.dim.use_parallel);
            nW_req  = calc.dim.num_workers;

            % Rule 1: if Nz < requested workers, do not parallelize (too little work)
            if Nz < nW_req
                use_par = false;
            end

            % Rule 2: if a pool already exists, use it (but still respect Rule 1)
            p = gcp('nocreate');
            if use_par
                if ~isempty(p)
                    nW_eff = p.NumWorkers;
                    if Nz < nW_eff
                        use_par = false;
                    end
                else
                    try
                        parpool('local', nW_req);
                    catch ME
                        warning('calc_ultrasound_field:ParpoolFailed', ...
                            'Failed to start parpool (%s). Falling back to serial.', ME.message);
                        use_par = false;
                    end
                end
            end
            % ====================================================================

            if use_par
                parfor iz = 1:Nz
                    z0 = z_obs(iz);
                    p1_vec = complex(zeros(nP,1));
                    p2_vec = complex(zeros(nP,1));

                    for s = 1:blk:nP
                        id = s:min(s+blk-1, nP);

                        phi1 = complex(zeros(numel(id),1));
                        phi2 = complex(zeros(numel(id),1));

                        for js = 1:src_blk:Nsrc
                            jd = js:min(js+src_blk-1, Nsrc);

                            dXj = xo(id) - Xs(jd).';
                            dYj = yo(id) - Ys(jd).';
                            dZj = (z0)   - Zs(jd).';
                            Rj  = sqrt(dXj.^2 + dYj.^2 + dZj.^2);
                            Rj(Rj < 1e-9) = 1e-9;

                            Gj = exp(1j*source.k1.*Rj) ./ (4*pi*Rj);
                            phi1 = phi1 + Gj * q1(jd);

                            Gj = exp(1j*source.k2.*Rj) ./ (4*pi*Rj);
                            phi2 = phi2 + Gj * q2(jd);
                        end

                        p1_vec(id) = -2j * source.medium.rho0 * source.w1 .* phi1;
                        p2_vec(id) = -2j * source.medium.rho0 * source.w2 .* phi2;
                    end

                    P1(:,:,iz) = reshape(p1_vec, Ny, Nx);
                    P2(:,:,iz) = reshape(p2_vec, Ny, Nx);
                end
            else
                for iz = 1:Nz
                    z0 = z_obs(iz);
                    p1_vec = complex(zeros(nP,1));
                    p2_vec = complex(zeros(nP,1));

                    for s = 1:blk:nP
                        id = s:min(s+blk-1, nP);

                        phi1 = complex(zeros(numel(id),1));
                        phi2 = complex(zeros(numel(id),1));

                        for js = 1:src_blk:Nsrc
                            jd = js:min(js+src_blk-1, Nsrc);

                            dXj = xo(id) - Xs(jd).';
                            dYj = yo(id) - Ys(jd).';
                            dZj = (z0)   - Zs(jd).';
                            Rj  = sqrt(dXj.^2 + dYj.^2 + dZj.^2);
                            Rj(Rj < 1e-9) = 1e-9;

                            Gj = exp(1j*source.k1.*Rj) ./ (4*pi*Rj);
                            phi1 = phi1 + Gj * q1(jd);

                            Gj = exp(1j*source.k2.*Rj) ./ (4*pi*Rj);
                            phi2 = phi2 + Gj * q2(jd);
                        end

                        p1_vec(id) = -2j * source.medium.rho0 * source.w1 .* phi1;
                        p2_vec(id) = -2j * source.medium.rho0 * source.w2 .* phi2;
                    end

                    P1(:,:,iz) = reshape(p1_vec, Ny, Nx);
                    P2(:,:,iz) = reshape(p2_vec, Ny, Nx);
                end
            end

            result.dim = struct();
            result.dim.method     = "DIM-Rayleigh";
            result.dim.dim_method = "rayleigh";
            result.dim.x          = x_obs;
            result.dim.y          = y_obs;
            result.dim.z          = z_obs;
            result.dim.p_f1       = P1;
            result.dim.p_f2       = P2;
            result.dim.f1         = source.f1;
            result.dim.f2         = source.f2;
            result.dim.block_size = blk;
            result.dim.src_block_size = src_blk;
            result.dim.use_parallel = use_par;
            if use_par
                result.dim.num_workers = calc.dim.num_workers;
            end

        case "asm"
            dx = source.dim.dx; dy = source.dim.dy;

            % --- uniformity check on 1D axes (correct) ---
            xg = source.dim.x_grid(:).';
            yg = source.dim.y_grid(:).';
            if numel(xg) < 2 || numel(yg) < 2
                error('ASM requires x_grid/y_grid length >= 2.');
            end
            if max(abs(diff(xg) - dx)) > 1e-12 || max(abs(diff(yg) - dy)) > 1e-12
                error('ASM requires uniform x_grid/y_grid consistent with dx/dy.');
            end

            % =============== ASM (with padding) ===============
            if ~isfield(source.dim,'dx') || ~isfield(source.dim,'dy')
                error('calc_ultrasound_field:ASMNeedDxDy', 'ASM needs source.dim.dx and source.dim.dy.');
            end
            if ~isfield(source.dim,'x_grid') || ~isfield(source.dim,'y_grid')
                error('calc_ultrasound_field:ASMNeedXYGrid', 'ASM needs source.dim.x_grid and source.dim.y_grid.');
            end

            % source plane (FFT grid)
            Vn0_f1 = source.dim.Vn_grid_f1;   % Ny0 x Nx0
            Vn0_f2 = source.dim.Vn_grid_f2;
            x0 = source.dim.x_grid(:).';      % 1 x Nx0
            y0 = source.dim.y_grid(:).';      % 1 x Ny0
            [Ny0, Nx0] = size(Vn0_f1);

            % ---- padding (center embed) ----
            pf = calc.asm.pad_factor;
            NyP = pf * Ny0;
            NxP = pf * Nx0;

            VnP_f1 = zeros(NyP, NxP);
            VnP_f2 = zeros(NyP, NxP);

            iy0 = floor((NyP - Ny0)/2) + (1:Ny0);
            ix0 = floor((NxP - Nx0)/2) + (1:Nx0);

            VnP_f1(iy0, ix0) = Vn0_f1;
            VnP_f2(iy0, ix0) = Vn0_f2;

            % padded coordinate axes
            x0c = x0( floor((Nx0+1)/2) );
            y0c = y0( floor((Ny0+1)/2) );

            x_pad = ((0:NxP-1) - floor((NxP+1)/2)) * dx + x0c;
            y_pad = ((0:NyP-1) - floor((NyP+1)/2)) * dy + y0c;

            % ---- propagate pressure directly from vn ----
            P1_full = local_asm_pressure_from_vn(VnP_f1, dx, dy, source.k1, source.medium.rho0, source.w1, z_obs, ...
                calc.asm.kzz_eps);
            P2_full = local_asm_pressure_from_vn(VnP_f2, dx, dy, source.k2, source.medium.rho0, source.w2, z_obs, ...
                calc.asm.kzz_eps);

            [x_out, y_out] = deal(x_pad, y_pad);
            [Xo, Yo] = meshgrid(x_out, y_out);

            result.dim = struct();
            result.dim.method     = "DIM-ASM";
            result.dim.dim_method = "asm";
            result.dim.x          = x_out;
            result.dim.y          = y_out;
            result.dim.z          = z_obs;
            result.dim.X          = Xo;
            result.dim.Y          = Yo;
            result.dim.p_f1       = P1_full;
            result.dim.p_f2       = P2_full;
            result.dim.f1         = source.f1;
            result.dim.f2         = source.f2;
            result.dim.dx         = dx;
            result.dim.dy         = dy;
            result.dim.pad_factor = pf;

        otherwise
            error('calc_ultrasound_field:BadDimMethod', ...
                'calc.dim.method must be ''rayleigh'' or ''asm''.');
    end
end

end

% ======================================================================
% helpers
% ======================================================================
function G = local_green_spec_analytic(k, kr_col, z_row, eps_phase, kz_min)
% Analytic Green spectrum: G = (1j/4pi) * exp(1j*kz*z) / kz
% with stable patch for small |kz*z| and floor for |kz|.

kz = sqrt(k.^2 - kr_col.^2);
idx = imag(kz) < 0;
kz(idx) = -kz(idx);

Nr = numel(kz);
Nz = numel(z_row);

KZ = kz * ones(1, Nz);
Z  = ones(Nr, 1) * reshape(z_row, 1, []);

% floor kz to avoid division blow-up
KZs = KZ;
mask_kz0 = abs(KZs) < kz_min;
if any(mask_kz0(:))
    KZs(mask_kz0) = kz_min .* exp(1j * angle(KZs(mask_kz0)));
end

% base expression
G = (1j/(4*pi)) * exp(1j * KZ .* Z) ./ KZs;

% small-phase patch: apply ONLY when kz is (almost) real
phase = abs(KZ .* Z);
rel_im = abs(imag(KZ)) ./ max(abs(KZ), kz_min);

mask = (phase < eps_phase) & (rel_im < 1e-6);
if any(mask(:))
    Gt = (1./KZs) + 1j*Z - (KZs .* (Z.^2))/2;
    G(mask) = (1j/(4*pi)) * Gt(mask);
end
end

function [x_obs, y_obs, z_obs] = local_parse_dim_grid(gd)
if isfield(gd,'x') && ~isempty(gd.x)
    x_obs = gd.x(:).';
else
    error('calc_ultrasound_field:BadGridX', 'grid.dim.x is required (vector).');
end

if isfield(gd,'y') && ~isempty(gd.y)
    y_obs = gd.y(:).';
else
    error('calc_ultrasound_field:BadGridY', 'grid.dim.y is required (vector).');
end

if isfield(gd,'z') && ~isempty(gd.z)
    z_obs = gd.z(:).';
else
    error('calc_ultrasound_field:BadGridZ', 'grid.dim.z is required (vector).');
end
end

function P = local_asm_pressure_from_vn(Vn0, dx, dy, k, rho0, w, z_list, kzz_eps)
% p_hat = (rho0 * w / kz) * Vn_hat * exp(1j*kz*z)

[Ny, Nx] = size(Vn0);

kx = (2*pi) * [0:floor(Nx/2), -ceil(Nx/2)+1:-1] / (Nx*dx);
ky = (2*pi) * [0:floor(Ny/2), -ceil(Ny/2)+1:-1] / (Ny*dy);
[KX, KY] = meshgrid(kx, ky);

kz = sqrt(k.^2 - KX.^2 - KY.^2);
kz(imag(kz)<0) = -kz(imag(kz)<0);

Vhat = fft2(Vn0);

coef = rho0 * w;

Nz = numel(z_list);
P = complex(zeros(Ny, Nx, Nz));

for iz = 1:Nz
    z = z_list(iz);
    H = exp(1j*kz*z) ./ (kz + kzz_eps);
    P(:,:,iz) = ifft2( coef * (Vhat .* H) );
end
end
