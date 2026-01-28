%% demo_scheme2_narrowband_correction_with_G_V.m
% 方案2（简单验证版）：
% - 保留粗网格 FFT-FHT（m_FHT）全谱计算：phi_coarse
% - 只在窄带（k_rho≈Re(k) 临界带）用“直接积分”补偿：phi_corr = phi_coarse + Delta_phi
% - 用“全k的直接积分(较密)”做 reference：phi_ref
%
% 依赖（你已有）：
%   - m_FHT.m
%   - solve_kappa0.m
%
% 你要求的点：
%   1) 点数调高（默认 N=32768，可改）
%   2) G(k_rho,z) 定义按你 King analytic spectrum：
%        G = (1j/4pi) * exp(1j*kz*z) / kz   (带稳定补丁)
%   3) F = G .* Vs，其中 Vs 来自一个简单口径速度 vs(rho) 的 Hankel 变换（用 m_FHT）
%
% 注意：
%   - 为兼容你的 m_FHT：x1 必须是 1×N 行向量；h 必须是 N×1 列向量
%   - 本脚本严格在调用侧保证尺寸，不修改你的函数

clear; clc; close all;

%% ===================== 0) parameters you may edit =====================
N  = 32768;         % FHT points (high)
miu = 0;             % Hankel order (0 for axisymmetric)

% physical truncations (choose similar to your pipeline)
rho_max = 2.0;       % field rho range of interest (m)
Nh  = 1.2 * rho_max; % rho truncation length for field (m)

% choose an ultrasound wavenumber (can set absorption imag part too)
c0 = 343;
f  = 42000;                 % Hz
k  = 2*pi*f/c0 + 1j*0.15;    % 1/m

% choose z (m)
z0 = 0.30;

% k_rho truncation length for field (1/m)
NH = 4 * real(k);           % similar to your NH_scale*(w/c0)

% aperture model for Vs (simple demo)
a  = 0.10;                  % aperture radius (m)
v0 = 0.172;                 % velocity amplitude (m/s)
Nh_v = 1.1 * a;             % rho truncation for velocity transform (m)
NH_v = NH;                  % keep same k_rho truncation

% narrowband correction settings
dk = 20;                    % half bandwidth around k0=Re(k) (1/m)
Qband = 4000;               % dense samples in band only (for Delta_phi)

% reference integral settings (full [0,NH] dense)
Qref  = 60000;              % dense samples for reference (can reduce if too slow)

% Green spectrum stability
eps_phase = 1e-3;
kz_min    = 1e-12;

%% ===================== 1) FHT grid (x1 etc.) =====================
n_FHT = (0:N-1).';
[a_solve, k0, x1, x0] = solve_kappa0(N, n_FHT);

% IMPORTANT: enforce shapes for your m_FHT
x1 = x1(:).';       % 1×N (required by your m_FHT implementation)
% grids
rho = (Nh  * x1).'; % N×1
kr  = (NH  * x1).'; % N×1  (column for our use)

rho_v = (Nh_v * x1).';  % N×1 (velocity sampling rho)
kr_v  = (NH_v * x1).';  % N×1

%% ===================== 2) Build vs(rho) and compute Vs(k_rho) =====================
% simple uniform piston: v0 inside radius a, 0 outside
vs = v0 * (rho_v <= a);
vs = vs(:);  % N×1

% forward Hankel (rho->k) using your m_FHT
Vs = m_FHT(vs, N, 1, Nh_v, NH_v, a_solve, x0, x1, k0, miu);
Vs = Vs(:); % N×1

%% ===================== 3) Build G(k_rho,z0) and F = G.*Vs on coarse grid =====================
G = green_spec_analytic(k, kr_v, z0, eps_phase, kz_min);  % N×1
F_coarse = (G .* Vs);                                     % N×1

%% ===================== 4) Coarse full inverse by FHT =====================
% phi(rho) = -4*pi * inverseHankel{ F(k) }  (your convention)
phi_coarse = -4*pi * m_FHT(F_coarse, N, 1, NH, Nh, a_solve, x0, x1, k0, miu);
phi_coarse = phi_coarse(:);

%% ===================== 5) Narrowband direct correction Delta_phi =====================
k0c = real(k);
kL = max(0, k0c - dk);
kU = min(NH, k0c + dk);

k_band = linspace(kL, kU, Qband).';  % Qband×1
% interpolate coarse spectrum onto band and also build a "better" band model
% For validation: emulate "true" band with extra sharp variation near k0c
% (replace this block with your real higher-res band if you have it)
F_on_band_from_coarse = interp1(kr_v, F_coarse, k_band, 'pchip', 0);

% Add an artificial sharp feature in band to mimic your “剧烈变化”
sigma = max(dk/15, 1e-6);
sharp = 1 + 3.0 * exp(-((k_band-k0c)/sigma).^2) .* exp(1j*0.6);
F_true_band = F_on_band_from_coarse .* sharp;

dF_band = F_true_band - F_on_band_from_coarse;  % only band difference

% trapezoid weights in k
wtrap = ones(Qband,1); wtrap([1 end]) = 0.5;
dk_step = (kU-kL)/(Qband-1);
wtrap = wtrap * dk_step;

integrand = (dF_band .* k_band) .* wtrap;  % Qband×1 (contains k dk)

Delta_phi = complex(zeros(N,1));
blk = 2048;  % rho blocking
for s = 1:blk:N
    ii = s:min(s+blk-1, N);
    rr = rho(ii);                 % nb×1
    Jm = besselj(miu, rr * (k_band.'));  % nb×Qband
    Delta_phi(ii) = -4*pi * (Jm * integrand);
end

phi_corr = phi_coarse + Delta_phi;

%% ===================== 6) Reference: full dense direct integral =====================
k_ref = linspace(0, NH, Qref).';
% build reference spectrum as: coarse spectrum + same sharp feature (band only)
F_ref = interp1(kr_v, F_coarse, k_ref, 'pchip', 0);
idb = (k_ref >= kL) & (k_ref <= kU);
% apply same sharp feature inside band (so reference contains the "true" band)
F_ref(idb) = F_ref(idb) .* (1 + 3.0 * exp(-((k_ref(idb)-k0c)/sigma).^2) .* exp(1j*0.6));

wtrap_ref = ones(Qref,1); wtrap_ref([1 end]) = 0.5;
dk_ref = NH/(Qref-1);
wtrap_ref = wtrap_ref * dk_ref;

integrand_ref = (F_ref .* k_ref) .* wtrap_ref;  % Qref×1

phi_ref = complex(zeros(N,1));
blk = 1024;
for s = 1:blk:N
    ii = s:min(s+blk-1, N);
    rr = rho(ii);
    Jm = besselj(miu, rr * (k_ref.'));      % nb×Qref
    phi_ref(ii) = -4*pi * (Jm * integrand_ref);
end

%% ===================== 7) plots + error =====================
figure;
plot(rho, 20*log10(abs(phi_ref)+1e-30), 'LineWidth', 1.2); hold on;
plot(rho, 20*log10(abs(phi_coarse)+1e-30), '--', 'LineWidth', 1.2);
plot(rho, 20*log10(abs(phi_corr)+1e-30), ':', 'LineWidth', 1.8);
grid on; xlabel('\rho (m)'); ylabel('|\phi| (dB)');
legend('ref (dense integral)','coarse FHT','coarse + narrowband correction','Location','best');
title(sprintf('Magnitude  N=%d, k0=%.2f, dk=%.2f', N, k0c, dk));

figure;
plot(rho, unwrap(angle(phi_ref)), 'LineWidth', 1.2); hold on;
plot(rho, unwrap(angle(phi_coarse)), '--', 'LineWidth', 1.2);
plot(rho, unwrap(angle(phi_corr)), ':', 'LineWidth', 1.8);
grid on; xlabel('\rho (m)'); ylabel('phase(\phi) (rad)');
legend('ref (dense integral)','coarse FHT','coarse + narrowband correction','Location','best');
title('Phase');

fprintf('RelErr coarse = %.3e, corrected = %.3e\n', ...
    norm(phi_coarse-phi_ref)/norm(phi_ref), norm(phi_corr-phi_ref)/norm(phi_ref));

%% ===================== local functions =====================
function G = green_spec_analytic(k, kr_col, z, eps_phase, kz_min)
% Your King analytic Green spectrum in k_rho domain:
%   G(k_rho,z) = (i/4pi) * exp(i*kz*z) / kz
% with:
%   - branch choice Im(kz)>=0
%   - floor on kz to avoid division blow-up
%   - small-phase Taylor patch when kz ~ real and |kz*z| small

kr_col = kr_col(:);
kz = sqrt(k.^2 - kr_col.^2);
kz(imag(kz)<0) = -kz(imag(kz)<0);

kzs = kz;
mask0 = abs(kzs) < kz_min;
if any(mask0)
    kzs(mask0) = kz_min .* exp(1j*angle(kzs(mask0)));
end

G = (1j/(4*pi)) .* exp(1j*kz*z) ./ kzs;

phase = abs(kz*z);
rel_im = abs(imag(kz)) ./ max(abs(kz), kz_min);
mask = (phase < eps_phase) & (rel_im < 1e-6);
if any(mask)
    % exp(i kz z)/kz ≈ 1/kz + i z - kz z^2/2
    Gt = (1./kzs) + 1j*z - (kzs*(z^2))/2;
    G(mask) = (1j/(4*pi)) .* Gt(mask);
end
end
