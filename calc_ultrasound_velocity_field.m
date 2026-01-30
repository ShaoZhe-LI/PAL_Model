function result = calc_ultrasound_velocity_field(source, medium, calc, method)
% =========================================================================
% CALC_ULTRASOUND_VELOCITY_FIELD
% - Compute ultrasonic particle-velocity field v = (v_rho, v_phi, v_z)
%   at f1/f2 on the SAME (rho,z) grid as your King/FHT pressure solver.
%
% - This follows the m-order King/Hankel velocity-field derivation:
%     v_rho, v_phi, v_z expressed by inverse Hankel transforms
%   (see Eq. (13) in your note). :contentReference[oaicite:0]{index=0}
%
% - Output is the "core" field on (rho,z) WITHOUT the exp(i*m*phi) factor.
%   If you need the full 3D cylindrical dependence:
%       v_total(rho,phi,z) = v_core(rho,z) * exp(1i*m*phi)
%
% INPUT
%   source, medium, calc : same as calc_ultrasound_field (passed to make_source_velocity)
%   method : 'king' | 'both' (default 'king')   (DIM velocity not implemented here)
%
% OUTPUT (struct)
%   result.source / result.medium / result.calc / result.fht
%   result.king:
%       .rho, .z, .m
%       .v_rho_f1, .v_phi_f1, .v_z_f1   (Nrho x Nz)
%       .v_rho_f2, .v_phi_f2, .v_z_f2
% =========================================================================

if nargin < 4 || isempty(method), method = 'king'; end
method = lower(strtrim(string(method)));
do_king = any(method == ["king","both"]);

if ~do_king
    error('calc_ultrasound_velocity_field:OnlyKing', ...
        'This function currently implements only the King/FHT velocity field.');
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

% -------------------- build consistent source + grids --------------------
[source, fht, ~] = make_source_velocity(source, medium, calc);

result = struct();
result.source = source;
result.medium = source.medium;
result.calc   = calc;
result.fht    = fht;

% =========================================================================
% KING / FHT velocity field
% =========================================================================
rho = fht.xh(:);            % (Nrho x 1) physical rho grid
z   = fht.z_ultra(:).';     % (1 x Nz)  (typically >=0)
Nrho = numel(rho);
Nz   = numel(z);

m = source.m_used;

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

% Make them (Nkr x 1) explicitly
Vs1 = Vs1(:);
Vs2 = Vs2(:);

% ---- kz branch (imag(kz) >= 0), and |z| / sgn(z) ----
Zabs = abs(z);                 % (1 x Nz)
sgnZ = sign(z); sgnZ(sgnZ==0) = 1;   % treat z=0 as + (consistent with your z>=0 grid)

% Build KZ matrices (Nkr x Nz)
KZ1 = local_kz_branch(source.k1, kr, Nz);
KZ2 = local_kz_branch(source.k2, kr, Nz);

% Floor kz to avoid blow-up
KZ1s = local_kz_floor(KZ1, calc.king.kz_min);
KZ2s = local_kz_floor(KZ2, calc.king.kz_min);

% exp(i*kz*|z|)
E1 = exp(1j * KZ1 .* (ones(Nkr,1)*Zabs));
E2 = exp(1j * KZ2 .* (ones(Nkr,1)*Zabs));

% -------------------------------------------------------------------------
% Core spectra per Eq.(13):
%   Common term (for v_phi and for v_rho):  Vs * exp(i*kz|z|) / kz
%   v_z uses:                              Vs * exp(i*kz|z|)
% -------------------------------------------------------------------------
Scom1 = (Vs1 * ones(1,Nz)) .* E1 ./ KZ1s;           % (Nkr x Nz)
Scom2 = (Vs2 * ones(1,Nz)) .* E2 ./ KZ2s;

Sz1   = (Vs1 * ones(1,Nz)) .* E1;                  % (Nkr x Nz)
Sz2   = (Vs2 * ones(1,Nz)) .* E2;

Sr1   = Scom1 .* (kr * ones(1,Nz));                % multiply by k_rho
Sr2   = Scom2 .* (kr * ones(1,Nz));

% -------------------------------------------------------------------------
% Inverse Hankel transforms:
%   - your pressure code uses: InvH_m(·) = m_FHT(·, N_FHT, Nz, NH, Nh, ..., m)
%   - v_z:  sgn(z) * InvH_m( Sz )
%   - v_phi: (m/rho) * InvH_m( Scom )
%   - v_rho: -(i/2) * [ InvH_{m-1}(Sr) - InvH_{m+1}(Sr) ]
% -------------------------------------------------------------------------

% v_z
Vz1 = m_FHT(Sz1, fht.N_FHT, Nz, fht.NH, fht.Nh, ...
    fht.a_solve, fht.x0, fht.x1, fht.k0, m);
Vz2 = m_FHT(Sz2, fht.N_FHT, Nz, fht.NH, fht.Nh, ...
    fht.a_solve, fht.x0, fht.x1, fht.k0, m);

Vz1 = Vz1 .* (ones(Nrho,1) * sgnZ);
Vz2 = Vz2 .* (ones(Nrho,1) * sgnZ);

% v_phi core (without exp(i m phi))
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

% v_rho core: use orders (m-1) and (m+1) directly (supports m<0)
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

% -------------------- pack --------------------
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
