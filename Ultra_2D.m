clear; clc; close all; % 清理工作区、关闭图形窗口、清空命令行

try
    if strcmp(profile('status').ProfilerStatus, 'on')
        profile off;
    end
catch
    % 旧版本兼容处理：直接尝试 profile off
    profile off;
end



% --- 1. 参数设置 ---
% 声学参数 (根据您的要求更新)
v0 = 0.12;              % 振幅
rho0 = 1.21;            % 空气密度 kg/m^3
f = 40e3;             % 频率 39.5 kHz
omega = 2 * pi * f;     % 角频率
c = 343;                % 声速 m/s
k = omega / c + 1j * AbsorpAttenCoef(f); % 复数波数
k = omega / c;
a = 0.05; % m


% --- 2. 二维声场计算区域设置 (xOz 平面) ---
% 定义声场计算平面的 X 和 Z 范围 (单位：米)
x_range = [0 0];
z_range = [0 3];

% 定义声场计算的采样间隔 (单位：米)
field_dx = 0.005; % X方向间隔 0.0025
field_dz = 0.005; % Z方向间隔 0.0025

field_dx = c/f/3; % X方向间隔 0.0025
field_dz = c/f/3; % Z方向间隔 0.0025

y_goal = 0;       % 计算平面的y坐标 (m)





% 计算，并输出计算时间和内存占用
profile on -memory
tic

[P_field, x_field, z_field] = Cal_Ultra_Field_2D(v0, rho0, f, c, k, a, x_range, z_range, field_dx, field_dz, y_goal);

toc
profile off
printProfileMemorySummary();

vars = whos;
total_bytes = sum([vars.bytes]);
total_MB = total_bytes / 1024^2;
total_GB = total_bytes / 1024^3;

fprintf('当前工作区所有变量共占用内存：%.2f MB (%.2f GB)\n', total_MB, total_GB);


% --- 3. 可视化 ---
figure;
pcolor(z_field,[fliplr(-x_field) x_field],...
    20*log10(abs([flipud(P_field);P_field])/2e-5));
colormap(MyColor('vik'));
shading interp;
clb = colorbar;
clb.Title.Interpreter = 'latex';
set(clb,'Fontsize',20);
% xlim([0 (z_audio(end))]);ylim([-rho_max rho_max]);
fontsize(gca,24,'points');
xlabel('$z$ (m)', 'Interpreter','latex','Fontsize',21);
ylabel('$x$ (m)', 'Interpreter','latex','Fontsize',21);
% xticks([0 2 4 6 8]);yticks([-2 -1 0 1 2]);
% pbaspect([z_audio(end), 2*rho_max, 1]);
set(gca, 'linewidth', 2);
set(gca, 'TickLabelInterpreter', 'latex');
clim([90,135]);