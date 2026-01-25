%% 完整示例程序：根据给定参数计算单点或多点声压

clear; clc; close all;


%% 导入King积分算法的相关参数及计算结果
load('data\parameters.mat');
load('data\pa_save.mat');   % pa_save = pa_W(1:r_interval:end,1:z_interval:end);
load('data\p1_save.mat');   % p1_save = p1(1:r_interval:end,1:z_interval:end);

xh = xh(1:r_interval:end);
z = z(1:z_interval:end);
z_audio = z_audio(1:z_interval:end);
P_King = p1_save(1 : 4*8 : length(xh), 1 : 4*4 : length(z));
P_field = squeeze(P_field);

err = log10(abs(P_King - P_field) ./ abs(P_field));     % 相对误差写作对数

%% 打开检测器
try
    if strcmp(profile('status').ProfilerStatus, 'on')
        profile off;
    end
catch
    % 旧版本兼容处理：直接尝试 profile off
    profile off;
end

%% 全局变量
isultra = true;
dis_coe = 16;           % 波长与离散间隔之比(16够用, 32更好)
isall = true;           % 是否计算整个声场
iscomp = true;          % 用于比较表示此处选点与King相同；否则是均匀取点
% block_size = 5e4;       % 块大小，默认10e4

%% 参数设置
% --- 1. 声学参数 ---
f1 = fu;
f2 = fu + fa;
w1 = 2 * pi * f1; k1 = w1 / c + 1j*AbsorpAttenCoef(f1);
w2 = 2 * pi * f2; k2 = w2 / c + 1j*AbsorpAttenCoef(f2);
wa = 2 * pi * fa; ka = wa / c + 1j*AbsorpAttenCoef(fa);

% --- 2. 三维声场计算区域设置 ---
if ~iscomp
    x_range = [-1.2, 1.2];     % x 范围 (m)
    y_range = [-1.2, 1.2];       % y 范围 (m)
    z_range = [0, 6];     % z 范围 (m)

    delta = c/f2/4;     % 1/4间隔
    dx = delta;             % x 采样间隔 (m) 0.0025
    dy = delta;             % y 采样间隔 (m) 0.0025
    dz = delta;             % z 采样间隔 (m) 0.0025

    x_field = x_range(1):dx:x_range(2);
    y_field = y_range(1):dy:y_range(2);
    z_field = z_range(1):dz:z_range(2);
end

if iscomp
    z_idx = 1 : 4*4 : length(z);
    r_idx = 1 : 4*8 : length(xh);
    x_field = xh(r_idx);
    y_field = 0;
    z_field = z(z_idx);
end

% --- 3. 计算声场，并输出时间和内存 ---
profile on -memory
tic

[P_field] = Cal_Ultra_Field_3D_block(v0, rho, f1, c, k1, a, x_field, y_field, z_field);

toc
profile off
printProfileMemorySummary();

vars = whos;
total_bytes = sum([vars.bytes]);
fprintf('当前工作区所有变量共占用内存：%.2f MB (%.2f GB)\n', total_bytes/1024^2, total_bytes/1024^3);

save('data\p_DIM.mat','P_field','x_field','y_field','z_field')

%%
% p = zeros(4096,709);
% for i =1:4096
%     for ii =1:709
%         p(i,ii) = P_field(i,1,ii);
%     end
% end
% save('p.mat','p')


%%
% err = log10(abs(p1_save - p)./abs(p));


% %% 二维音频声场
% issave = 1;
% 
% 
% x_show = 1;
% z_show = 4;
% x_grid = x_field;
% y_grid = y_field;
% z_grid = z_field;
% 
% x_idx = abs(x_grid) < x_show+0.1;
% x_idx = vec2ndim(x_idx, 1);
% y_idx = y_grid == 0;
% y_idx = vec2ndim(y_idx, 2);
% z_idx = abs(z_grid) < z_show+0.1;
% z_idx = vec2ndim(z_idx, 3);
% % z_idx = (-0.1 < z_grid) & (z_grid < z_show+0.1);
% 
% x_grid_show = vec2ndim(x_grid(x_idx), 1);
% y_grid_show = (y_grid(y_idx));
% z_grid_show = (vec2ndim(z_grid(z_idx), 2));
% prs_a_show = permute(P_field(x_idx, y_idx, z_idx), [1 3 2]);
% 
% 
% spl_show = prs2spl(prs_a_show);
% % spl_show = spl_show - max(spl_show,[],"all");
% % spl_show = abs(prs_a_show)./max(abs(prs_a_show),[],"all");
% fig11 = figure;
% % pcolor(z_grid_show, x_grid_show, spl_show - max(spl_show,[],"all")); 
% pcolor(z_grid_show, x_grid_show, spl_show);
% shading interp;
% colormap(slanCM('vik'));
% colorbar;
% c = colorbar;
% c.Label.Interpreter = 'latex';
% c.TickLabelInterpreter = 'latex';
% c.Title.String = 'SPL (dB)';
% c.Title.Interpreter = 'latex';
% xlim([0 z_show]);
% ylim([-x_show x_show]);
% colormap(slanCM('vik'));
% 
% xlabel('$z$ (m)', 'Interpreter', 'latex');
% ylabel('$x$ (m)', 'Interpreter', 'latex');
% % title('f1 超声场 $x$-$z$ 平面声压级分布', 'Interpreter', 'latex');
% 
% set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex');
% clim([60 100])
% c.Ticks = [60:10:100];
% 
% fig_paper = gcf;
% fig_paper.PaperPositionMode = 'auto';
% fig_pos = fig_paper.PaperPosition;
% fig_paper.PaperSize = [fig_pos(3) fig_pos(4)];
% pbaspect([4, 2, 1])
% 
% 
% if issave == 1
%     exportgraphics(fig11, fullfile(['Audio_xOz', '.pdf']), 'ContentType', 'vector');
% end
% 
% 
% %% xOy
% z_goal = 2;
% [~, idx_z2] = min(abs(z_grid - z_goal));
% pa_xy_show = (P_field(:,:,idx_z2));
% y_show1 = y_grid;
% x_show1 = (x_grid.');
% 
% 
% fig12 = figure;
% pcolor(y_show1, x_show1, prs2spl(pa_xy_show));
% shading interp;
% colormap(slanCM('vik'));
% colorbar;
% c = colorbar;
% c.Label.Interpreter = 'latex';
% c.TickLabelInterpreter = 'latex';
% c.Title.String = 'SPL (dB)';
% c.Title.Interpreter = 'latex';
% xlim([-x_show x_show]);
% ylim([-x_show x_show]);
% 
% xlabel('$x$ (m)', 'Interpreter', 'latex');
% ylabel('$y$ (m)', 'Interpreter', 'latex');
% 
% set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex');
% clim([50 80])
% % c.Ticks = [45:10:75];
% 
% fig_paper = gcf;
% fig_paper.PaperPositionMode = 'auto';
% fig_pos = fig_paper.PaperPosition;
% fig_paper.PaperSize = [fig_pos(3) fig_pos(4)];
% pbaspect([1, 1, 1]);
% 
% if issave == 1
%     exportgraphics(fig12, fullfile(['Audio_xOy_Abs', '.pdf']), 'ContentType', 'vector');
% end
% 
% 
% 
% p_draw = prs2spl(pa_xy_show(:,y_show1==0));
% fig121 = figure;
% plot(x_show1,p_draw,'LineWidth', 2);
% xlabel('$x$ (m)', 'Interpreter', 'latex');
% ylabel('SPL (dB)', 'Interpreter', 'latex');
% set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex');
% max_lim = max(p_draw) - mod(max(p_draw),10) + 10;
% ylim([max_lim-10 max_lim])
% ylim([64 70]);yticks([64:2:70]);
% xlim([-1 1])
% pbaspect([2.5, 2, 1]);
% 
% if issave == 1
%     exportgraphics(fig121, fullfile( ['Audio_xOy_Abs_y0', '.pdf']), 'ContentType', 'vector');
% end
