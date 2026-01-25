%% 完整示例程序：根据给定参数计算单点或多点声压

clear; clc; close all;


%% 导入King积分算法的相关参数及计算结果
load('data\parameters.mat');
load('data\pa_save.mat');   % pa_save = pa_W(1:r_interval:end,1:z_interval:end);
load('data\p1_save.mat');   % p1_save = p1(1:r_interval:end,1:z_interval:end);

xh = xh(1:r_interval:end);
z = z(1:z_interval:end);
z_audio = z_audio(1:z_interval:end);

a = 0.15;

m=5;

%% 全局变量
isultra = true;
dis_coe = 16;           % 波长与离散间隔之比(16够用, 32更好)
isall = false;           % 是否计算整个声场


%% ==== 参数设置 ====
% 用于比较的King结果


% 观测点坐标 (Nx3)
if isall
    if isultra
        p_King = p1_save;
        z_idx = 1 : 8 : length(z);
    else
        p_King = pa_save;
        z_idx = 1 : 8 : length(z_audio);
    end
    r_idx = 1 : 16 : length(xh);

else
    if isultra
        p_King = p1_save;
        z_idx = [1 2 5 10 15 20 30 50 75 100 150 200 300 500 700];
        z_idx = 97;
    else
        p_King = pa_save;
        z_idx = [1 2 5 10 15 20 30 50 75 100 150];
    end
    r_idx = 1:2:length(xh);
end

x_obs0 = xh(r_idx);
y_obs0 = 0;
z_obs0 = z(z_idx);
p_comp = p_King(r_idx,z_idx);


[X, Y, Z] = meshgrid(x_obs0, y_obs0, z_obs0);
r_points = [X(:), Y(:), Z(:)];  % Nx3 矩阵

p_comp = p_King(r_idx,z_idx);

% 声学参数
f1 = fu;
f2 = fu + fa;
w1 = 2 * pi * f1; k1 = w1 / c + 1j*AbsorpAttenCoef(f1);
w2 = 2 * pi * f2; k2 = w2 / c + 1j*AbsorpAttenCoef(f2);
wa = 2 * pi * fa; ka = wa / c + 1j*AbsorpAttenCoef(fa);


%% 直接积分法计算
% --- 2. 定义圆形声源的几何参数 ---
source_radius = a; % 圆形声源的半径 (m)
source_center_x = 0; % 圆心 x 坐标 (m)
source_center_y = 0; % 圆心 y 坐标 (m)
source_z_plane = 0;  % 声源所在的 z 平面 (m)


% --- 3. 定义离散间隔 ---
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
U_vals_at_sources = exp(1i * m * phi_s);

% --- 6. 定义每个离散声源点代表的面积元 (dA) ---
% 在直角坐标系下，每个离散点代表的面积元是 dx * dy
dA = dx * dy;

% --- 7. 计算声压 ---
Nr = length(x_obs0)
Nz = length(z_obs0)
p = zeros(Nr,Nz); % 存储每个观测点的声压
pause(2);

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
        p(idr, idz) = -2i * rho * v0 * w1 * current_sum_integrand;
    end
end

% p_all = [p_comp;p];

err = log10(abs(p - p_comp) ./ abs(p));     % 相对误差写作对数
err_mean = mean(err);


%% 绘制误差图像
% figure
% pcolor(z_obs0, x_obs0, err)
% shading interp
% colorbar

%%
r_boundary = 0.25;
[~,idr_boundary]=min(abs(x_obs0 - r_boundary));
idx_pump = 1508;
% p = besselh(m, 1, x_obs0);

close all
figure('Name','abs')
plot(x_obs0(1:idr_boundary), 20*log10(abs(p(1:idr_boundary))),'LineWidth',1.2);
hold on
plot(x_obs0(1:idr_boundary), abs(p_comp(1:idr_boundary)),'LineWidth',2);
yyaxis right
plot(x_obs0(1:idr_boundary), err(1:idr_boundary),'LineWidth',1.2);
legend('DIM','King','err');

figure('Name','phs')
plot(x_obs0(1:idr_boundary), angle(p(1:idr_boundary)),'LineWidth',1.2);
hold on
plot(x_obs0, angle(p_comp),'LineWidth',2);
legend('DIM','King');

figure('Name','phs_unwrap')
p_unwrap = unwrap(angle(p));
p_comp_unwrap = unwrap(angle(p_comp));
plot(x_obs0(1:idr_boundary), p_unwrap(1:idr_boundary),'LineWidth',1.2);
hold on
plot(x_obs0(1:idr_boundary), p_comp_unwrap(1:idr_boundary),'LineWidth',2);
legend('DIM','King');
% plot(x_obs0(idx_pump), p_unwrap(idx_pump), 'o');
% label_text = sprintf('(%.4f, %.4f)', x_obs0(idx_pump), p_unwrap(idx_pump));
% text(x_obs0(idx_pump) + 0.005, p_unwrap(idx_pump) + 0.1, label_text, ...
%     'FontSize', 10, ...
%     'Color', 'red', ...
%     'HorizontalAlignment', 'left', ...
%     'VerticalAlignment', 'bottom');


theta_fig = 0:0.01:2*pi;
rho_fig = xh(1:2:end);
[TH, R] = meshgrid(theta_fig,rho_fig);
[X,Y] = pol2cart(TH, R);
XY = X.^2+Y.^2;
[XY_r,XY_c]=size(XY);

pa_K_fig0 = p(:,1).*exp(1j*m*theta_fig);
p_DIM_fig = pa_K_fig0;
% pa_K_fig(find(XY>r_boundary^2)) = NaN;
pa_K_fig0 = p_comp(:,1).*exp(1j*m*theta_fig);
p_King_fig = pa_K_fig0;


fig = figure('Name','p_DIM','position',[200 200 1000 1000]);
pcolor(X,Y,angle(p_DIM_fig)/pi);
% pcolor(X,Y,abs(pa_K_fig));
shading interp
xlim([-1.2*r_boundary 1.2*r_boundary]);
ylim([-1.2*r_boundary 1.2*r_boundary]);
pbaspect([1 1 1])
colormap('hsv')
% colormap(MyColor('vik'))
clim([-1 1])
colorbar
set(gca,'position', [0.15 0.15 0.7 0.7]);
set(fig,'PaperSize',[16,16]);


fig = figure('Name','p_King','position',[200 200 1000 1000]);
pcolor(X,Y,angle(p_King_fig)/pi);
% pcolor(X,Y,abs(pa_K_fig));
shading interp
xlim([-1.2*r_boundary 1.2*r_boundary]);
ylim([-1.2*r_boundary 1.2*r_boundary]);
pbaspect([1 1 1])
colormap('hsv')
% colormap(MyColor('vik'))
clim([-1 1])
colorbar
set(gca,'position', [0.15 0.15 0.7 0.7]);
set(fig,'PaperSize',[16,16]);


fig = figure('Name','p_DIM_amp','position',[200 200 1000 1000]);
pcolor(X,Y,abs(p_DIM_fig));
% pcolor(X,Y,abs(pa_K_fig));
shading interp
xlim([-1.2*r_boundary 1.2*r_boundary]);
ylim([-1.2*r_boundary 1.2*r_boundary]);
pbaspect([1 1 1])
% colormap('hsv')
colormap(MyColor('vik'))
% clim([-1 1])
colorbar
set(gca,'position', [0.15 0.15 0.7 0.7]);
set(fig,'PaperSize',[16,16]);

fig = figure('Name','p_King_amp','position',[200 200 1000 1000]);
pcolor(X,Y,abs(p_King_fig));
% pcolor(X,Y,abs(pa_K_fig));
shading interp
xlim([-1.2*r_boundary 1.2*r_boundary]);
ylim([-1.2*r_boundary 1.2*r_boundary]);
pbaspect([1 1 1])
% colormap('hsv')
colormap(MyColor('vik'))
% clim([-1 1])
colorbar
set(gca,'position', [0.15 0.15 0.7 0.7]);
set(fig,'PaperSize',[16,16]);