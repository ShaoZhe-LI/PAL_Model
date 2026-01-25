

clear;clc;close all;


%%
isradial = 1;       %是否径向（否则就是轴向）,1使用z_idx0

r_idx0 = 11000;      % 8000是0.02  11000是0.116
z_idx0 = 30;       % 250是1m


m1_v = [5];
name_file = 'result\';

a = 0.15;
v0 = 0.172;
c = 343;
rho = 1.21;
beta = 1.2;

fu = 42e3;

% fu = 10e3;
% a = 18 / (2*pi*fu/343);

fa = 2e3;
f1 = fu;
f2 = fu + fa;

N_FHT = 16384 * 2;
delta = c/f2*4;
rho_max = 5;
zu_max = 4;
za_max = 2 + delta;
r_boundary = 0.2;

isprofile = 'm-king';
isdraw = 0;
isSPL = 1;
z_goal_phi = 1;

results = cell(1,length(m1_v));
max_results = zeros(1,length(m1_v));
min_results = zeros(1,length(m1_v));

%% ===================== m 扫描主循环 =====================

clc; close all

tic

m1 = m1_v(1);
m2 = m1;
ma = m2 - m1;
iscomp = 0;

%% ===================== 基本频率参数 =====================
w1 = 2*pi*f1;
w2 = 2*pi*f2;
wa = 2*pi*fa;

k1 = w1/c + 1j*AbsorpAttenCoef(f1);
k2 = w2/c + 1j*AbsorpAttenCoef(f2);
ka = wa/c + 1j*AbsorpAttenCoef(fa);

%% ===================== FHT 参数 =====================
Nh = 1.2*rho_max;
NH = 1.2*w2/c;
NH = 5*w2/c;
n_FHT = 0:N_FHT-1;

[a_solve,k0,x1,x0] = solve_kappa0(N_FHT,n_FHT);
xh = (x1*Nh).';
xH = (x1*NH).';

%% 画图参数
r_interval = 8 * 4;
z_interval = 4;

if isradial
    r_idx = 1:r_interval:N_FHT;
    z_idx = z_idx0;
else
    z_idx = 1:z_interval:length(z);
    r_idx = r_idx0;
end

%% ===================== 速度分布 =====================
syms rho_v
switch isprofile
    case 'uniform'
        m1 = 0; m2 = 0;
        vs_sym1 = v0*(heaviside(rho_v)-heaviside(rho_v-a));
        vs_sym2 = vs_sym1;
    case 'focus'
        m1 = 0; m2 = 0;
        F = 0.2;
        vs_sym1 = v0*exp(-1j*real(k1)*sqrt(rho_v.^2+F^2)) ...
            .* (heaviside(rho_v)-heaviside(rho_v-a));
        vs_sym2 = v0*exp(-1j*real(k2)*sqrt(rho_v.^2+F^2)) ...
            .* (heaviside(rho_v)-heaviside(rho_v-a));
    case 'm-king'
        vs_sym1 = v0*(heaviside(rho_v)-heaviside(rho_v-a));
        vs_sym2 = vs_sym1;
        % vs_sym1 = v0*(rho_v^5)*(heaviside(rho_v)-heaviside(rho_v-a));
end

Nh_v = 1.1*a;
NH_v = NH;
xh_v = (x1*Nh_v).';
xH_v = (x1*NH_v).';

vs_f1 = matlabFunction(vs_sym1);
vs1 = vs_f1(xh_v);
vs_f2 = matlabFunction(vs_sym2);
vs2 = vs_f2(xh_v);

Vs1 = m_FHT(vs1,N_FHT,1,Nh_v,NH_v,a_solve,x0,x1,k0,m1);
Vs2 = m_FHT(vs2,N_FHT,1,Nh_v,NH_v,a_solve,x0,x1,k0,m2);

%% ===================== 超声场计算 =====================
z = 0:delta:zu_max; Nz = length(z);
z_audio = 0:delta:za_max; Nza = length(z_audio);

rs=sqrt(xh.^2+z.^2);
g1=exp(1j*k1*rs)./(4*pi*rs);
G1_FHT = m_FHT(g1,N_FHT,Nz,Nh,NH,a_solve,x0,x1,k0,0);
% 这种方法更好


kr = xH;                      % 你的 k_rho 采样点
kz = sqrt(k1.^2 - kr.^2);      % 主值
idx = imag(kz) < 0;
kz(idx) = -kz(idx);
G1_ana = 1j/(4*pi) * exp(1i*(kz*z)) ./ (kz*ones(1,numel(z)));   % Nk×Nz

eps0 = 1e-3;
mask = abs(kz*z) < eps0;                       % Nk×Nz
KZmat = kz*ones(1,numel(z));                   % Nk×Nz
Zmat  = ones(numel(kz),1)*z;                   % Nk×Nz
G_taylor = (1./KZmat) + 1j*Zmat - (KZmat.*(Zmat.^2))/2;  % Nk×Nz
G1_ana(mask) = 1j/(4*pi) * G_taylor(mask);
% z=0 列（若存在）保持为 1/kz
iz0 = (z == 0);
if any(iz0)
    G1_ana(:,iz0) = 1j/(4*pi) ./ (kz * ones(1,sum(iz0)));
end

% G1_ana_without = 1j/(4*pi) * exp(1j*kz.*z)./kz;

G1_ana_plot = G1_ana(r_idx,z_idx);
G1_FHT_plot = G1_FHT(r_idx,z_idx);
F1 = G1_ana.*Vs1;
clear G1_ana
phi1 = -1j*m_FHT(F1,N_FHT,Nz,NH,Nh,a_solve,x0,x1,k0,m1);
clear F1
p1_ana = 1j*rho*c*real(k1)*phi1 * 4*pi;
clear phi1

F1 = G1_FHT.*Vs1;
clear G1_FHT
phi1 = -1j*m_FHT(F1,N_FHT,Nz,NH,Nh,a_solve,x0,x1,k0,m1);
clear F1
p1_FHT = 1j*rho*c*real(k1)*phi1 * 4*pi;
clear phi1
% p1_audio = p1(:,1:Nza);


%% 直接积分法
x_obs0 = xh(r_idx);
y_obs0 = 0;
z_obs0 = z(z_idx);
[X, Y, Z] = meshgrid(x_obs0, y_obs0, z_obs0);
r_points = [X(:), Y(:), Z(:)];  % Nx3 矩阵




% --- 2. 定义圆形声源的几何参数 ---
source_radius = a; % 圆形声源的半径 (m)
source_center_x = 0; % 圆心 x 坐标 (m)
source_center_y = 0; % 圆心 y 坐标 (m)
source_z_plane = 0;  % 声源所在的 z 平面 (m)


% --- 3. 定义离散间隔 ---
dis_coe = 8;
lambda = c / f1;        % 当前频率的波长
dx = lambda / dis_coe; % 离散间隔，等于半波长
dy = lambda / dis_coe;
dz = 0; % 如果声源在平面上，这个值可以为0或者保持dz不变

% --- 4. 优化后的离散圆形声源 (使用向量化操作) ---
disp('正在离散圆形声源 (优化中)...');

% 定义覆盖圆形声源的矩形区域
% 稍微扩大范围以确保包含边界，防止浮点数误差导致边界点遗漏
x_min_rect = source_center_x - source_radius - dx;
x_max_rect = source_center_x + source_radius + dx;
y_min_rect = source_center_y - source_radius - dy;
y_max_rect = source_center_y + source_radius + dy;

% 生成密集的直角坐标网格点
x_grid = x_min_rect:dx:x_max_rect;
y_grid = y_min_rect:dy:y_max_rect;

% 使用 meshgrid 生成所有可能的 (X, Y) 坐标对
[X_s_all, Y_s_all] = meshgrid(x_grid, y_grid);

% 将 Z 坐标设置为声源所在平面 (对于所有点都相同)
Z_s_all = source_z_plane * ones(size(X_s_all));

% 计算每个网格点到圆心的二维距离 (向量化操作)
distance_to_center = sqrt((X_s_all - source_center_x).^2 + ...
    (Y_s_all - source_center_y).^2);

% 使用逻辑索引筛选出在圆形区域内的点
% `in_circle_idx` 是一个逻辑矩阵，其中为 true 的位置表示该点在圆内
in_circle_idx = (distance_to_center <= source_radius);

% 提取圆形区域内的 X, Y, Z 坐标
x_s_circular = X_s_all(in_circle_idx);
y_s_circular = Y_s_all(in_circle_idx);
z_s_circular = Z_s_all(in_circle_idx); % Z 坐标也会被筛选，但值都相同

% 将筛选出的点合并到 circular_source_points 矩阵中
circular_source_points = [x_s_circular(:), y_s_circular(:), z_s_circular(:)];

num_source_points = size(circular_source_points, 1);


if num_source_points == 0
    error('未生成任何声源离散点。请检查声源参数和离散步长。');
end
disp(['圆形声源已离散为 ', num2str(num_source_points), ' 个点。']);

% --- 5. 定义声源面振动函数 u(rs_coord) ---
% 替换您原有的 U = ones(size(Rho_s));
% 对于每个离散声源点，定义其振动值。
% 示例：均匀振动面，所有声源点的 u_val 都为 1
% 如果需要复杂分布，可以在这里定义 u_vals = f(circular_source_points);
Xs = circular_source_points(:,1);
Ys = circular_source_points(:,2);
phi_s = atan2(Ys, Xs);   % [-pi, pi]
rho_s = sqrt(Xs.^2 + Ys.^2);
m = m1;
U_vals_at_sources = exp(1i * m * phi_s);
% U_vals_at_sources = (rho_s.^5)

% --- 6. 定义每个离散声源点代表的面积元 (dA) ---
% 在直角坐标系下，每个离散点代表的面积元是 dx * dy
dA = dx * dy;

% --- 7. 计算声压 ---
Nr = length(x_obs0)
Nz = length(z_obs0)
p_DIM = zeros(Nr,Nz); % 存储每个观测点的声压
% pause(2);

for idr = 1:Nr
    idr
    for idz = 1:Nz
        x_obs = x_obs0(idr);
        y_obs = y_obs0;
        z_obs = z_obs0(idz);

        current_sum_integrand = 0; % 用于累加所有声源点的贡献

        for s_idx = 1:num_source_points
            % 当前声源点坐标
            x_s = circular_source_points(s_idx, 1);
            y_s = circular_source_points(s_idx, 2);
            z_s = circular_source_points(s_idx, 3);

            % 声源面振动值 u(r_s)
            u_val = U_vals_at_sources(s_idx); % 从预定义的 U_vals_at_sources 获取

            % 距离计算 R = |r - r_s|
            R = sqrt( (x_obs - x_s)^2 + (y_obs - y_s)^2 + (z_obs - z_s)^2 );

            % Green函数 g(r; r_s; ω_i)
            G = exp(1i * k1 * R) / (4 * pi * R);

            % 处理 R 接近零的情况 (避免除以零)
            if R < 1e-9 % R 非常小，可能表示观测点和声源点重合
                G = 0; % 或者一个更合适的处理方式，根据物理模型而定
            end

            % 累加被积函数 u * g * dA
            current_sum_integrand = current_sum_integrand + u_val * G * dA;
        end

        % 计算声压
        p_DIM(idr, idz) = -2i * rho * v0 * w1 * current_sum_integrand;
    end
end

p_DIM = 1j * p_DIM;
% p_all = [p_comp;p];


% %% 角谱法
% %% ===================== 角谱法（ASM）计算：z = z(z_idx) 平面相位 =====================
% % 约定：与上文一致使用 exp(+i*k*R) 的外辐射解；频域默认 e^{-i*w*t}
% % 输入：与 DIM/FHT 同一声源：v_n(rho,phi)=vs_f1(rho)*exp(i*m1*phi)，disk 外为 0
% % 输出：p_ASM(x,y,z0) 以及相位图
% 
% % ---------- 0) 取目标 z 平面 ----------
% z0 = z(z_idx);        % 你当前选的 z 平面（isradial=1 时 z_idx 是标量）
% 
% % ---------- 1) 设定源面离散网格（必须是二维均匀网格，供 FFT） ----------
% % 建议：网格覆盖至少 ±r_boundary（最好更大一点以降低边界截断影响）
% dx_ASM = dx;          % 直接用你 DIM 的 dx（但注意：dx 必须足够小，建议 <= lambda/8）
% dy_ASM = dx_ASM;
% 
% Lx = 2.5 * r_boundary * 2;   % 覆盖范围（可按需增大）
% Ly = Lx;
% 
% Nx = 2^nextpow2( round(Lx/dx_ASM) );   % FFT 友好
% Ny = Nx;
% 
% x_ASM = ( (-Nx/2):(Nx/2-1) ) * dx_ASM;
% y_ASM = ( (-Ny/2):(Ny/2-1) ) * dy_ASM;
% [X_ASM, Y_ASM] = meshgrid(x_ASM, y_ASM);
% 
% RHO_ASM = sqrt(X_ASM.^2 + Y_ASM.^2);
% PHI_ASM = atan2(Y_ASM, X_ASM);
% 
% % ---------- 2) 源面法向速度 v_n(x,y,0)（与 DIM/FHT 一致） ----------
% % 关键：vs_f1(rho) 的定义必须与你 FHT 那边一致（是否包含 v0）
% v_radial_ASM = vs_f1(RHO_ASM);                 % disk 内分布（若 vs_f1 不含 v0，则在这里乘 v0）
% v_n_xy = v_radial_ASM .* exp(1i*m1*PHI_ASM);   % m 阶角向项
% v_n_xy(RHO_ASM > a) = 0;                       % disk 外为 0（刚性挡板）
% 
% % ---------- 3) 构造 (kx,ky) 以及 kz（选择辐射分支：imag(kz) >= 0） ----------
% k = k1;     % 用与你超声场一致的复波数（含吸收时可直接用 k1）
% % 频域采样（单位：rad/m）
% kx = 2*pi * ( (-Nx/2):(Nx/2-1) ) / (Nx*dx_ASM);
% ky = 2*pi * ( (-Ny/2):(Ny/2-1) ) / (Ny*dy_ASM);
% [KX, KY] = meshgrid(kx, ky);
% 
% KZ = sqrt(k.^2 - KX.^2 - KY.^2);     % 主值
% idx = imag(KZ) < 0;                  % 选取使场不发散的分支
% KZ(idx) = -KZ(idx);
% 
% % ---------- 4) 角谱关系：由 v_z(kx,ky,0) 得到 p(kx,ky,0) ----------
% % 由动量方程频域形式：v_z = (1/(i*rho*w)) * ∂p/∂z
% % 对每个平面波分量：∂/∂z -> i*KZ  => v_z = (KZ/(rho*w)) * p  => p = rho*w * v_z / KZ
% % 注意：KZ->0 附近要做数值保护
% V_spec = fftshift( fft2( ifftshift(v_n_xy) ) ) * (dx_ASM*dy_ASM);  % 连续近似的尺度（含面积元）
% epsK = 1e-8;
% P0_spec = rho*w1 * V_spec ./ (KZ + epsK);    % z=0 平面压力谱
% 
% % ---------- 5) 传播到 z0：P(kx,ky,z0)=P0*exp(i*KZ*z0)，再 IFFT ----------
% Pz_spec = P0_spec .* exp(1i*KZ*z0);
% 
% p_ASM = fftshift( ifft2( ifftshift(Pz_spec) ) ) / (dx_ASM*dy_ASM); % 回到 x-y 平面复声压
% 
% % （可选）如果你发现整体幅值与 DIM/FHT 有常数比例差，通常来自 FFT 尺度/定义差；
% % 相位一般不会受这个常数影响。


%% 画图比较
% xh(r_idx)
close all
if isradial
    z_goal_phi = z(z_idx);
    theta_fig = 0:0.01:2*pi;
    rho_fig = xh(1:r_interval:end);
    [~,index_z]=min(abs(z_audio-z_goal_phi));
    p1_ana_2d = p1_ana(1:r_interval:end,z_idx).*exp(1j*m1*theta_fig);
    p1_FHT_2d = p1_FHT(1:r_interval:end,z_idx).*exp(1j*m1*theta_fig);
    p_DIM_2d = p_DIM(:).*exp(1j*m1*theta_fig);
    [TH, R] = meshgrid(theta_fig,rho_fig);
    [X,Y] = pol2cart(TH, R);
    XY = X.^2+Y.^2;
    [XY_r,XY_c]=size(XY);
end

if isradial
    figure('Name','Radial');
    plot(xH(r_idx), abs(G1_FHT_plot), '--', 'LineWidth', 2);
    hold on
    plot(xH(r_idx), abs(G1_ana_plot), 'LineWidth', 1.2);
    yyaxis right
    plot(xH(r_idx), abs(Vs1(r_idx)));
    hold off
    legend('FHT', 'ana', 'V');
else
    figure('Name','Axial');
    plot(z(z_idx), abs(G1_FHT_plot), '--', 'LineWidth', 2);
    hold on
    plot(z(z_idx), abs(G1_ana(r_idx,z_idx)), 'LineWidth', 1.2);
    yyaxis right
    plot(xH(r_idx), abs(Vs1(r_idx)));
    hold off
    legend('FHT', 'ana', 'V');
end

% err = log10(abs(p_DIM - p1(r_idx,z_idx)) ./ abs(p_DIM));     % 相对误差写作对数
% err_mean = mean(err);

if isradial
    figure('Name','p_Radial');
    plot(xh(r_idx), abs(p1_ana(r_idx,z_idx)), '--', 'LineWidth', 3);
    hold on
    plot(xh(r_idx), abs(p1_FHT(r_idx,z_idx)), '-.', 'LineWidth', 2);
    plot(xh(r_idx), abs(p_DIM), 'LineWidth', 1.2);
    hold off
    legend('FHT-G-ana', 'FHT-G-FHT', 'DIM');
else
    figure('Name','p_Axial');
    plot(z(z_idx), abs(p1_ana(r_idx,z_idx)), '--', 'LineWidth', 3);
    plot(xh(r_idx), abs(p1_FHT(r_idx,z_idx)), '-.', 'LineWidth', 2);
    hold on
    plot(z(z_idx), abs(p_DIM), 'LineWidth', 1.2);
    hold off
    legend('FHT-G-ana', 'FHT-G-FHT', 'DIM');
end


% p1_ana_phs = unwrap(angle(p1_ana(r_idx,z_idx)));
% p1_FHT_phs = unwrap(angle(p1_FHT(r_idx,z_idx)));
% p_DIM_phs = unwrap(angle(p_DIM));
p1_ana_phs = (angle(p1_ana(r_idx,z_idx)));
p1_FHT_phs = (angle(p1_FHT(r_idx,z_idx)));
p_DIM_phs = (angle(p_DIM));
% p1_ana_phs = p1_ana_phs - (p1_ana_phs(round(0.5*length(p1_ana_phs))) - p_DIM_phs(round(0.5*length(p_DIM_phs))));
% p1_FHT_phs = p1_FHT_phs - (p1_FHT_phs(round(0.5*length(p1_ana_phs))) - p_DIM_phs(round(0.5*length(p_DIM_phs))));

if isradial
    figure('Name','p_phs_Radial');
    plot(xh(r_idx), p1_ana_phs/pi, '--', 'LineWidth', 3);
    hold on
    plot(xh(r_idx), p1_FHT_phs/pi, '-.', 'LineWidth', 2);
    plot(xh(r_idx), p_DIM_phs/pi, 'LineWidth', 1.2);
    hold off
    legend('FHT-G-ana', 'FHT-G-FHT', 'DIM');
else
    figure('Name','p_phs_Axial');
    plot(z(z_idx), p1_ana_phs/pi, '--', 'LineWidth', 3);
    plot(z(z_idx), p1_FHT_phs/pi, '-.', 'LineWidth', 2);
    hold on
    plot(z(z_idx), p_DIM_phs/pi, 'LineWidth', 1.2);
    hold off
    legend('FHT-G-ana', 'FHT-G-FHT', 'DIM');
end


if isradial
    fig = figure('Name','p_ana_phs','position',[200 200 1000 1000]);
    pcolor(X,Y,angle(p1_ana_2d)/pi);
    % pcolor(X,Y,abs(pa_K_fig));
    shading flat
    xlim([-1.2*r_boundary 1.2*r_boundary]);
    ylim([-1.2*r_boundary 1.2*r_boundary]);
    pbaspect([1 1 1])
    colormap('hsv')
    % colormap(MyColor('vik'))
    clim([-1 1])
    colorbar
    set(gca,'position', [0.15 0.15 0.7 0.7]);

    fig = figure('Name','p_FHT_phs','position',[200 200 1000 1000]);
    pcolor(X,Y,angle(p1_FHT_2d)/pi);
    % pcolor(X,Y,abs(pa_K_fig));
    shading flat
    xlim([-1.2*r_boundary 1.2*r_boundary]);
    ylim([-1.2*r_boundary 1.2*r_boundary]);
    pbaspect([1 1 1])
    colormap('hsv')
    % colormap(MyColor('vik'))
    clim([-1 1])
    colorbar
    set(gca,'position', [0.15 0.15 0.7 0.7]);

    fig = figure('Name','p_DIM_phs','position',[200 200 1000 1000]);
    pcolor(X,Y,angle(p_DIM_2d)/pi);
    % pcolor(X,Y,abs(pa_K_fig));
    shading flat
    xlim([-1.2*r_boundary 1.2*r_boundary]);
    ylim([-1.2*r_boundary 1.2*r_boundary]);
    pbaspect([1 1 1])
    colormap('hsv')
    % colormap(MyColor('vik'))
    clim([-1 1])
    colorbar
    set(gca,'position', [0.15 0.15 0.7 0.7]);
end


% %% ===================== 画图：特定 z0 平面相位（与前面风格一致） =====================
% % 只看 r <= r_boundary 的区域
% mask_rb = (RHO_ASM <= r_boundary);
% 
% phase_ASM = angle(p_ASM)/pi;
% phase_ASM(~mask_rb) = NaN;
% 
% figure('Name','ASM phase at z0','position',[250 150 900 800]);
% pcolor(X_ASM, Y_ASM, phase_ASM);
% shading flat;
% axis equal;
% xlim([-1.2*r_boundary, 1.2*r_boundary]);
% ylim([-1.2*r_boundary, 1.2*r_boundary]);
% colormap('hsv');
% clim([-1 1]);
% colorbar;
% title(sprintf('ASM phase (angle/\\pi) at z = %.4f m, m = %d', z0, m1));
% 
% % （可选）与 DIM / FHT 同一条径向切线比较相位（y=0）
% [~, iy0] = min(abs(y_ASM-0));
% figure('Name','Phase cut y=0','position',[250 150 900 450]);
% plot(x_ASM, angle(p_ASM(iy0,:))/pi, 'LineWidth', 1.2);
% grid on;
% xlabel('x (m)'); ylabel('phase/\pi');
% title(sprintf('ASM phase cut (y=0) at z = %.4f m', z0));
% xlim([-r_boundary, r_boundary]);
